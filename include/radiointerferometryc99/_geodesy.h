// Currently this is only for internal use!
#ifndef __GEODESY_H
#define __GEODESY_H

#include <stddef.h>

#define WGS84_A_METERS 6378137
#define WGS84_F (((double)1.0) / ( 298257223563LL / 1000000000 ))
#define WGS84_F_INV 298.257223563

typedef struct {
	double a;
	double b;
	double f;
	double e2;
} geodesy_t;

static inline void geodesy_from_ab(geodesy_t* geo, double a, double b) {
	geo->a = a;
	geo->b = b;
	geo->f = 1 - geo->b/geo->a;
	geo->e2 = geo->f*(2-geo->f);
}

static inline void geodesy_from_af(geodesy_t* geo, double a, double f) {
	double b = a*(1 - f);
	geodesy_from_ab(geo, a, b);
}

static inline void geodesy_from_af_inv(geodesy_t* geo, double a, double f_inv) {
	double b = a*(1 - 1.0/f_inv);
	geodesy_from_ab(geo, a, b);
}

#endif