#include "uvh5/uvh5_calc.h"

float UVH5calc_julian_date_from_unix(float unix_sec) {
	return (unix_sec / DAYSEC) / 2440587.5;
}

float UVH5calc_julian_date_from_guppi_param(
	float tbin,
	size_t sampleperblk,
	size_t piperblk,
	size_t synctime,
	size_t pktidx
) {
	 float tperblk = sampleperblk * tbin;
	 float tperpktidx = tperblk / piperblk;
	 return UVH5calc_julian_date_from_unix(synctime + tperpktidx*pktidx);
}

void UVH5calc_ha_dec_rad(
	double ra_rad,
	double dec_rad,
	double longitude_rad,
	double latitude_rad,
	double altitude,
	double timemjd,
	double dut1,
	double* hour_angle_rad,
	double* declination_rad
) {
	double aob, zob, rob, eo;
	eraAtco13(
		ra_rad, dec_rad,
		0, 0, 0, 0,
		timemjd, 0,
		dut1,
		longitude_rad, latitude_rad, altitude,
		0, 0,
		0, 0, 0, 0,
		&aob, &zob, hour_angle_rad,
    declination_rad, &rob, &eo
	);
}

/*
 * https://github.com/liberfa/erfa/blob/master/src/gst06a.c#L44-L47
 * This uses UT1 for both UT1 and TT, which results
 * in an error on the order of 100 microarcseconds or approximately
 * 7 microseconds.
 */
double UVH5calc_lst(
	double timemjd,
	double dut1
) {
	return eraGst06a(timemjd, dut1, timemjd, dut1);
}

float UVH5calc_hypotenuse_f(float* position, int dims) {
	double sum = 0.0;
	while(--dims > 0) {
		sum += position[dims]*position[dims];
	}
	return (float) sqrt(sum);
}

double UVH5calc_hypotenuse(double* position, int dims) {
	double sum = 0.0;
	while(--dims > 0) {
		sum += position[dims]*position[dims];
	}
	return sqrt(sum);
}

void UVH5calc_frame_translate(double* positions, int position_count, double translation[3]) {
	while(--position_count >= 0)
	{
		positions[position_count*3+0] += translation[0];
		positions[position_count*3+1] += translation[1];
		positions[position_count*3+2] += translation[2];
	}
}

/*
 * https://github.com/JuliaGeo/Geodesy.jl/blob/dc2b3bd4d73a5fb4ed6f2f9c5462763ac54e5196/src/transformations.jl#L175-L188
 */
void UVH5calc_ecef_from_lla(
	double ecef[3],
	const double longitude_rad,
	const double latitude_rad,
	const double altitude,
	const geodesy_t* geo
) {
    double sin_phi = sin(latitude_rad);
		double cos_phi = cos(latitude_rad);
    double sin_lambda = sin(longitude_rad);
		double cos_lambda = cos(longitude_rad);

    double N = geo->a / sqrt(1 - geo->e2 * sin_phi * sin_phi);  // Radius of curvature (meters)

    ecef[0] = (N + altitude) * cos_phi * cos_lambda;
    ecef[1] = (N + altitude) * cos_phi * sin_lambda;
    ecef[2] = (N * (1 - geo->e2) + altitude) * sin_phi;
}

/*
 * Clockwise (right-hand curl) rotations
 * y' = cos*vec.y + sin*vec.z
 * z' = -sin*vec.y + cos*vec.z
 */
static inline void _rotate_around_x_cached_trig(
	double vec[3],
	double sin_val,
	double cos_val
) {
	double y, z;
	y = vec[1];
	z = vec[2];
	vec[1] = cos_val*y - sin_val*z;
	vec[2] = sin_val*y + cos_val*z;
}

/*
 * Clockwise (right-hand curl) rotations
 * x' = cos*vec.x - sin*vec.z
 * z' = sin*vec.x + cos*vec.z
 */
static inline void _rotate_around_y_cached_trig(
	double vec[3],
	double sin_val,
	double cos_val
) {
	double x, z;
	x = vec[0];
	z = vec[2];
	vec[0] = cos_val*x + sin_val*z;
	vec[2] = -sin_val*x + cos_val*z;
}

/*
 * Clockwise (right-hand curl) rotations
 * x' = cos*vec.x - sin*vec.y
 * y' = sin*vec.x + cos*vec.y
 */
static inline void _rotate_around_z_cached_trig(
	double vec[3],
	double sin_val,
	double cos_val
) {
	double x, y;
	x = vec[0];
	y = vec[1];
	vec[0] = cos_val*x - sin_val*y;
	vec[1] = sin_val*x + cos_val*y;
}

/*
 * Subtracts ECEF(LLA, WGS84) from positions.
 */
void UVH5calc_position_to_xyz_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
) {
	double ecef[3];
	geodesy_t wgs84 = {0};
	geodesy_from_af_inv(&wgs84, WGS84_A_METERS, WGS84_F_INV);

	UVH5calc_ecef_from_lla(
		ecef,
		longitude_rad,
		latitude_rad,
		altitude,
		&wgs84
	);
	ecef[0] *= -1.0;
	ecef[1] *= -1.0;
	ecef[2] *= -1.0;
	UVH5calc_frame_translate(positions, position_count, ecef);
}

/*
 * Adds ECEF(LLA, WGS84) to positions.
 */
void UVH5calc_position_to_ecef_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
) {
	double ecef[3];
	geodesy_t wgs84 = {0};
	geodesy_from_af_inv(&wgs84, WGS84_A_METERS, WGS84_F_INV);

	UVH5calc_ecef_from_lla(
		ecef,
		longitude_rad,
		latitude_rad,
		altitude,
		&wgs84
	);
	UVH5calc_frame_translate(positions, position_count, ecef);
}

/*
 * https://github.com/david-macmahon/RadioInterferometry.jl/blob/3a084d47919ddb422be109c41faade9d365c2a35/src/RadioInterferometry.jl#L659-L662
 * Rotates the ENU frame anticlockwise about the East (i.e. first)
 * axis by `lat_rad`, producing a (East,Z,X') frame, then rotates that frame
 * anticlockwise about the Z (i.e. second) axis by `-lon_rad`, producing a
 * (Y,Z,X) frame which is then permuted to (X,Y,Z).
 */
void UVH5calc_position_to_xyz_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
) {
	double sin_longitude = sin(longitude_rad);
	double cos_longitude = cos(longitude_rad);
	double sin_latitude = sin(latitude_rad);
	double cos_latitude = cos(latitude_rad);
	double tmp;
	while(--position_count >= 0)
	{
		// RotX(longitude) anti-clockwise
		_rotate_around_x_cached_trig(
			positions + position_count*3,
			-sin_latitude,
			cos_latitude
		);
		// RotY(longitude) clockwise
		_rotate_around_y_cached_trig(
			positions + position_count*3,
			sin_longitude,
			cos_longitude
		);
		// Permute (YZX) to (XYZ)
		tmp = positions[position_count*3 + 2]; // save X
		positions[position_count*3 + 2] = positions[position_count*3 + 1]; // move Z
		positions[position_count*3 + 1] = positions[position_count*3 + 0]; // move Y
		positions[position_count*3 + 0] = tmp; // move X
	}
}

/*
 * https://github.com/david-macmahon/RadioInterferometry.jl/blob/3a084d47919ddb422be109c41faade9d365c2a35/src/RadioInterferometry.jl#L487-490
 * Rotates the XYZ frame anticlockwise about the Z (i.e. third)
 * axis by `lon_rad`, producing a (X',East,Z) frame, then rotates that frame
 * anticlockwise about the E (i.e. second) axis by `-lat_rad`, producing a
 * (U,E,N) frame which is then permuted to (E,N,U).
 */
void UVH5calc_position_to_enu_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
) {
	double sin_longitude = sin(longitude_rad);
	double cos_longitude = cos(longitude_rad);
	double sin_latitude = sin(latitude_rad);
	double cos_latitude = cos(latitude_rad);
	double tmp;
	while(--position_count >= 0)
	{
		// RotZ(longitude) anti-clockwise
		_rotate_around_z_cached_trig(
			positions + position_count*3,
			-sin_longitude,
			cos_longitude
		);
		// RotY(longitude) clockwise
		_rotate_around_y_cached_trig(
			positions + position_count*3,
			sin_latitude,
			cos_latitude
		);
		// Permute (UEN) to (ENU)
		tmp = positions[position_count*3 + 0];
		positions[position_count*3 + 0] = positions[position_count*3 + 1];
		positions[position_count*3 + 1] = positions[position_count*3 + 2];
		positions[position_count*3 + 2] = tmp;
	}
}

/*
 * Effects `ecef -> xyz -> enu`.
 */
void UVH5calc_position_to_enu_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
) {
  UVH5calc_position_to_xyz_frame_from_ecef(
		positions,
		position_count,
		longitude_rad,
		latitude_rad,
		altitude
	);
  UVH5calc_position_to_enu_frame_from_xyz(
		positions,
		position_count,
		longitude_rad,
		latitude_rad,
		altitude
	);
}

/*
 * Effects `enu -> xyz -> ecef`.
 */
void UVH5calc_position_to_ecef_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
) {
  UVH5calc_position_to_xyz_frame_from_enu(
		positions,
		position_count,
		longitude_rad,
		latitude_rad,
		altitude
	);
  UVH5calc_position_to_ecef_frame_from_xyz(
		positions,
		position_count,
		longitude_rad,
		latitude_rad,
		altitude
	);
}

/*
 * https://github.com/david-macmahon/RadioInterferometry.jl/blob/3a084d47919ddb422be109c41faade9d365c2a35/src/RadioInterferometry.jl#L584-L589
 * Rotates the ENU frame anticlockwise about the East (i.e. first)
 * axis by `lat_rad`, producing a (East,Z,X') frame, then rotates that frame
 * anticlockwise about the Z (i.e. second) axis by `-ha_rad`, producing a
 * (U,Z,X") frame, then rotates anticlockwise about the U (i.e. first) axis by
 * `-dec_rad`, producing the (U,V,W) frame where U is east, V is north, and W is
 * in the direction of projection.
 */
void UVH5calc_position_to_uvw_frame_from_enu(
	double* positions,
	int position_count,
	double hour_angle_rad,
	double declination_rad,
	double latitude_rad
) {	
	double sin_hour_angle = sin(hour_angle_rad);
	double cos_hour_angle = cos(hour_angle_rad);
	double sin_declination = sin(declination_rad);
	double cos_declination = cos(declination_rad);
	double sin_latitude = sin(latitude_rad);
	double cos_latitude = cos(latitude_rad);

	while(--position_count >= 0) {
		 // anti-clockwise
		_rotate_around_x_cached_trig(
			positions + position_count*3,
			-sin_latitude,
			cos_latitude
		);
		 // clockwise
		_rotate_around_y_cached_trig(
			positions + position_count*3,
			sin_hour_angle,
			cos_hour_angle
		);
		 // clockwise
		_rotate_around_x_cached_trig(
			positions + position_count*3,
			sin_declination,
			cos_declination
		);
	}
}

/*
 * https://github.com/david-macmahon/RadioInterferometry.jl/blob/3a084d47919ddb422be109c41faade9d365c2a35/src/RadioInterferometry.jl#L422-427
 * Rotates the XYZ frame anticlockwise about the Z (i.e. third)
 * axis by `lon_rad-ha_rad`, producing a (X',U,Z) frame, then rotates that frame
 * anticlockwise about the U (i.e. second) axis by `-dec_rad`, producing an
 * (W,U,V) frame which is then permuted to (U,V,W) where U is east, V is north,
 * and W is in the direction of the given hour angle and declination as seen from
 * the given longitude.
 */
void UVH5calc_position_to_uvw_frame_from_xyz(
	double* positions,
	int position_count,
	double hour_angle_rad,
	double declination_rad,
	double longitude_rad
) {
	double sin_long_minus_hangle = sin(longitude_rad-hour_angle_rad);
	double cos_long_minus_hangle = cos(longitude_rad-hour_angle_rad);
	double sin_declination = sin(declination_rad);
	double cos_declination = cos(declination_rad);
	double tmp;

	while (--position_count >= 0) {
		// RotZ(long-ha) anti-clockwise
		_rotate_around_z_cached_trig(
			positions + 3*position_count,
			-sin_long_minus_hangle,
			cos_long_minus_hangle
		);
		// RotY(declination) clockwise
		_rotate_around_y_cached_trig(
			positions + 3*position_count,
			sin_declination,
			cos_declination
		);
		// Permute (WUV) to (UVW)
		tmp = positions[3*position_count + 0]; // save W
		positions[3*position_count + 0] = positions[3*position_count + 1]; // move U
		positions[3*position_count + 1] = positions[3*position_count + 2]; // move U
		positions[3*position_count + 2] = tmp; // move U
	}
}