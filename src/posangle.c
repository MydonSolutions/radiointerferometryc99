#include <string.h>

#include "radiointerferometryc99.h"

// all transcribed and surmised from pyuvdata@4fba712cz:
// https://github.com/RadioAstronomySoftwareGroup/pyuvdata/blob/4fba712c51d638ed12b6669b157c2926b152b2a0/pyuvdata/utils.py#L3191C3-L3201C3

int _itrs_transform_app_to_icrs(
  double* time_jd,
  double* app_ra_radians,
  double* app_dec_radians,
  double* pm_x_arcsec,
  double* pm_y_arcsec,
  double* ut1_utc_sec,
  size_t count,
  double longitude_rad,
  double latitude_rad,
  double altitude,
  double* icrs_ra,
  double* icrs_dec
) {
  /*
  Transform a set of coordinates in topocentric/apparent to ICRS coordinates.

  This utility uses either astropy or erfa to calculate the ICRS  coordinates of
  a given set of apparent source coordinates. These coordinates are most typically
  used for defining the celestial/catalog position of a source. Note that at present,
  this is only implemented in astropy and pyERFA, although it could hypothetically
  be extended to NOVAS at some point.

  Parameters
  ----------
  time_jd :
    Julian dates to calculate coordinate positions for.
  app_ra_radians :
    ICRS RA of the celestial target, expressed in units of radians.
  app_dec_radians :
    ICRS Dec of the celestial target, expressed in units of radians.
  count :
    The number of elements behind each of the 3 pointers.

  Returns
  -------
  icrs_ra :
    ICRS right ascension coordinates, in units of radians. Taken to be allocated.
  icrs_dec :
    ICRS declination coordinates, in units of radians. Taken to be allocated.

  : int
    Zero if success, otherwise `(index+1)*10+errcode` encoding the index of the
    erroneous element and the errorcodes:
    - 0 being dubious year
    - 1 being unacceptable date.
    - [2, 8] being iers_get() errcode + 3
  */
  double rbpn_matrix[3][3];
  double cip_x, cip_y;
  double cio_s;
  double eqn_org;
  int rv;
  for (size_t i = 0; i < count; i++) {
    eraPnm06a(time_jd[i], 0, rbpn_matrix);
    eraBpn2xy(rbpn_matrix, &cip_x, &cip_y);
    cio_s = eraS06(time_jd[i], 0, cip_x, cip_y);
    eqn_org = eraEors(rbpn_matrix, cio_s);

    // Observed to ICRS via ERFA
    // status: +1 = dubious year (Note 4), 0 = OK, -1 = unacceptable date
    rv = eraAtoc13(
      "R",
      app_ra_radians[i] + eqn_org,
      app_dec_radians[i],
      time_jd[i], 0,
      ut1_utc_sec[i],
      longitude_rad,
      latitude_rad,
      altitude,
      pm_x_arcsec[i] * (RADIOINTERFEROMETERY_PI/(180 * 3600)), // convert arcsec to radian
      pm_y_arcsec[i] * (RADIOINTERFEROMETERY_PI/(180 * 3600)), // convert arcsec to radian
      0, // atm pressure, used for refraction (ignored)
      0, // amb temperature, used for refraction (ignored)
      0, // rel humidity, used for refraction (ignored)
      0, // wavelength, used for refraction (ignored)
      icrs_ra + i,
      icrs_dec + i
    );
    if (rv != 0) {
      // {-1, +1} -> {1, 0}
      return (i+1)*10+((rv+2)%3);
    }
  }
  
  return 0;
}

int calc_itrs_icrs_frame_pos_angle(
  double* time_jd,
  double* app_ra_radians,
  double* app_dec_radians,
  size_t count,
  double longitude_rad,
  double latitude_rad,
  double altitude,
  double offset_pos,
  const char* iers_filepath,
  double* pos_angle
) {
  /*
  Calculate an position angle given apparent position and reference frame.

  Accesses IERS data then calls `calc_itrs_icrs_frame_pos_angle_with_pm_and_ut1_utc`
  */

  // Get IERS data, which is needed for highest precision
  radiointerferometry_iers_record_t iers_rec = {0};
  double* pm_x_arcsec = malloc(count*sizeof(double));
  double* pm_y_arcsec = malloc(count*sizeof(double));
  double* ut1_utc_sec = malloc(count*sizeof(double));
  int rv;
  
  for (size_t i = 0; i < count; i++) {
    iers_rec.mjd = time_jd[i] - 2400000.5;
    rv = radiointerferometry_iers_get(
      iers_filepath,
      &iers_rec
    );
    if (rv != 0) {
      return (i+1)*10+(rv+3);
    }
    pm_x_arcsec[i] = iers_rec.pm_x_a;
    pm_y_arcsec[i] = iers_rec.pm_y_a;
    ut1_utc_sec[i] = iers_rec.ut1_utc_a;
  }

  rv = calc_itrs_icrs_frame_pos_angle_with_pm_and_ut1_utc(
    time_jd,
    app_ra_radians,
    app_dec_radians,
    pm_x_arcsec,
    pm_y_arcsec,
    ut1_utc_sec,
    count,
    longitude_rad,
    latitude_rad,
    altitude,
    offset_pos,
    pos_angle
  );

  free(pm_x_arcsec);
  free(pm_y_arcsec);
  free(ut1_utc_sec);
  return rv;
}

int calc_itrs_icrs_frame_pos_angle_with_pm_and_ut1_utc(
  double* time_jd,
  double* app_ra_radians,
  double* app_dec_radians,
  double* pm_x_arcsec,
  double* pm_y_arcsec,
  double* ut1_utc_sec,
  size_t count,
  double longitude_rad,
  double latitude_rad,
  double altitude,
  double offset_pos,
  double* pos_angle
) {
  /*
  Calculate an position angle given apparent position and reference frame.

  This function is used to determine the position angle between the great
  circle of declination in apparent coordinates, versus that in a given
  reference frame. Note that this is slightly different than parallactic
  angle, which is the difference between apparent declination and elevation.

  The telescope frame is taken to be 'itrs'.
  The reference frame is taken to be 'icrs'.

  Note that this computation is intensive. As such the provided arrays should
  ideally only express unique combinations, with results being repeated
  appropriately by the caller.

  Parameters
  ----------
  time_jd :
    Julian dates to calculate coordinate positions for.
  app_ra_radians :
    ICRS RA of the celestial target, expressed in units of radians.
  app_dec_radians :
    ICRS Dec of the celestial target, expressed in units of radians.
  count :
    The number of elements behind each of the 3 pointers.
  offset_pos : 
    Distance of the offset position used to calculate the frame PA. Recommendation is
    `PI/360` such that the PA is determined over a 1 deg arc. Must be > 0.

  Returns
  -------
  frame_pa :
    Array of position angles, in units of radians. Taken to be allocated.

  : int
    Zero if success, otherwise `(index+1)*10+errcode` encoding the index of the
    erroneous element and the errorcodes:
    - -1 being early date
    - 0 being dubious year
    - 1 being unacceptable date.
    - >=2 being iers_get() errcode + 1
  */

  double *_time_jd, *_app_ra_radians, *_app_dec_radians, *_pm_x_arcsec, *_pm_y_arcsec, *_ut1_utc_sec, *icrs_ra, *icrs_dec;
  size_t array_size = sizeof(double)*count;

  _time_jd = malloc(2*array_size);
  _app_ra_radians = malloc(2*array_size);
  _app_dec_radians = malloc(2*array_size);
  _pm_x_arcsec = malloc(2*array_size);
  _pm_y_arcsec = malloc(2*array_size);
  _ut1_utc_sec = malloc(2*array_size);
  icrs_ra = malloc(2*array_size);
  icrs_dec = malloc(2*array_size);
  
  memcpy(_time_jd, time_jd, array_size);
  memcpy(_time_jd+count, time_jd, array_size);
  memcpy(_app_ra_radians, app_ra_radians, array_size);
  memcpy(_app_ra_radians+count, app_ra_radians, array_size);
  memcpy(_pm_x_arcsec, pm_x_arcsec, array_size);
  memcpy(_pm_x_arcsec+count, pm_x_arcsec, array_size);
  memcpy(_pm_y_arcsec, pm_y_arcsec, array_size);
  memcpy(_pm_y_arcsec+count, pm_y_arcsec, array_size);
  memcpy(_ut1_utc_sec, ut1_utc_sec, array_size);
  memcpy(_ut1_utc_sec+count, ut1_utc_sec, array_size);
  for (size_t i = 0; i < count; i++) {
    _app_dec_radians[i] = app_dec_radians[i] - offset_pos;
    // Wrap the positions if they happen to go over the poles
    if (-_app_dec_radians[i] > (RADIOINTERFEROMETERY_PI/2.0)) {
      _app_dec_radians[i] = (RADIOINTERFEROMETERY_PI/2.0) - _app_dec_radians[i];
    }
    
    _app_dec_radians[count+i] = app_dec_radians[i] + offset_pos;
    // Wrap the positions if they happen to go over the poles
    if (_app_dec_radians[count+i] > (RADIOINTERFEROMETERY_PI/2.0)) {
      _app_dec_radians[count+i] = (RADIOINTERFEROMETERY_PI/2.0) - _app_dec_radians[count+i];
    }
  }

  // Run the set of offset coordinates through the "reverse" transform. The two offset
  // positions are concat'd together to help reduce overheads
  int rv = _itrs_transform_app_to_icrs(
    _time_jd,
    _app_ra_radians,
    _app_dec_radians,
    _pm_x_arcsec,
    _pm_y_arcsec,
    _ut1_utc_sec,
    2*count,
    longitude_rad,
    latitude_rad,
    altitude,
    icrs_ra,
    icrs_dec
  );
  if (rv != 0) {
    count = ((rv/10)-1)/2;
  }

  // Use the pas function from ERFA to calculate the position angle. The negative sign
  // is here because we're measuring PA of app -> frame, but we want frame -> app.
  for (size_t i = 0; i < count; i++) {
    pos_angle[i] = -eraPas(
      icrs_ra[i], icrs_dec[i], icrs_ra[count+i], icrs_dec[count+i]
    );
  }

  free(_time_jd);
  free(_app_ra_radians);
  free(_app_dec_radians);
  free(_pm_x_arcsec);
  free(_pm_y_arcsec);
  free(_ut1_utc_sec);
  free(icrs_ra);
  free(icrs_dec);
  return rv;
}