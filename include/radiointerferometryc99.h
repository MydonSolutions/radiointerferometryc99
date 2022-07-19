#ifndef _UVH5_CALC_H
#define _UVH5_CALC_H

#include <stddef.h>
#include <math.h>
#include "_geodesy.h"
#include "erfa.h"
#include "erfam.h"

#define RADIOINTERFEROMETERY_DAYSEC ERFA_DAYSEC
#define RADIOINTERFEROMETERY_PI 3.14159265358979323846
#define RADIOINTERFEROMETERY_C 299792458.0

static inline double calc_deg2rad(double deg) {return (deg/180)*RADIOINTERFEROMETERY_PI;};

enum position_frames {
	FRAME_ENU,
	FRAME_XYZ,
	FRAME_ECEF,
	FRAME_UVW
};

static inline double calc_modified_from_julian_date(double julian_date) {return  julian_date + 2400000.5;}
static inline double calc_julian_date_from_modified(double modified_jd) {return  modified_jd - 2400000.5;}

double calc_julian_date_from_unix(double unix_sec);

double calc_epoch_seconds_from_guppi_param(
	const double tbin,
	const size_t sampleperblk,
	const size_t piperblk,
	const size_t synctime,
	const size_t pktidx
);

double calc_julian_date_from_guppi_param(
	const double tbin,
	const size_t sampleperblk,
	const size_t piperblk,
	const size_t synctime,
	const size_t pktidx
);

void calc_ha_dec_rad(
	double ra_rad,
	double dec_rad,
	double longitude_rad,
	double latitude_rad,
	double altitude,
	double timemjd,
	double dut1,
	double* hour_angle_rad,
	double* declination_rad
);

double calc_lst(double timemjd, double dut1);

float calc_hypotenuse_f(float* position, int dims);
double calc_hypotenuse(double* position, int dims);

void calc_frame_translate(double* positions, int position_count, double translation[3]);

void calc_independent_astrom(
	double longitude_rad,
	double latitude_rad,
	double altitude,
	double timemjd,
	double dut1,
    eraASTROM* astrom
);

void calc_ha_dec_rad_with_independent_astrom(
	double ra_rad,
	double dec_rad,
    eraASTROM* astrom,
	double* hour_angle_rad,
	double* declination_rad
);

void calc_ecef_from_lla(
	double ecef[3],
	const double longitude_rad,
	const double latitude_rad,
	const double altitude,
	const geodesy_t* geo
);

void calc_position_to_xyz_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

void calc_position_to_ecef_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

void calc_position_to_xyz_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
);

void calc_position_to_enu_frame_from_xyz(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude // Not used
);

void calc_position_to_enu_frame_from_ecef(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

void calc_position_to_ecef_frame_from_enu(
	double* positions,
	int position_count,
	double longitude_rad,
	double latitude_rad,
	double altitude
);

void calc_position_to_uvw_frame_from_enu(
	double* positions,
	int position_count,
	double hour_angle_rad,
	double declination_rad,
	double latitude_rad
);

void calc_position_to_uvw_frame_from_xyz(
	double* positions,
	int position_count,
	double hour_angle_rad,
	double declination_rad,
	double longitude_rad
);

void calc_position_delays(
	double* positions_xyz_in_uvw_out,
	int position_count,
	int reference_position_index,
	double hour_angle_rad,
	double declination_rad,
	double longitude_rad,
	double* delays
);

#endif
