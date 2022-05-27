#ifndef _UVH5_CALC_H
#define _UVH5_CALC_H

#include <stddef.h>
#include <math.h>
#include "_geodesy.h"
#include "erfa.h"
#include "erfam.h"

#define DAYSEC ERFA_DAYSEC

static inline double calc_deg2rad(double deg) {return (deg/180)*M_PI;};

enum position_frames {
	FRAME_ENU,
	FRAME_XYZ,
	FRAME_ECEF,
	FRAME_UVW
};

float calc_julian_date_from_unix(float unix_sec);

float calc_julian_date_from_guppi_param(
	float tbin,
	size_t sampleperblk,
	size_t piperblk,
	size_t synctime,
	size_t pktidx
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

#endif
