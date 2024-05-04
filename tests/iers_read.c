#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "radiointerferometryc99.h"

void print(int rv, radiointerferometry_iers_record_t *iers_rec) {
  printf("return code:  %d\n", rv);
  printf("--------------\n");

  printf("year:         %d\n", iers_rec->year); // int
  printf("month:        %d\n", iers_rec->month); // int
  printf("day:          %d\n", iers_rec->day); // int
  printf("mjd:          %f\n", iers_rec->mjd); // double
  printf("polpmflag_a:  %d\n", iers_rec->polpmflag_a); // bool
  printf("pm_x_a:       %f\n", iers_rec->pm_x_a); // double
  printf("e_pm_x_a:     %f\n", iers_rec->e_pm_x_a); // double
  printf("pm_y_a:       %f\n", iers_rec->pm_y_a); // double
  printf("e_pm_y_a:     %f\n", iers_rec->e_pm_y_a); // double
  printf("ut1flag_a:    %d\n", iers_rec->ut1flag_a); // bool
  printf("ut1_utc_a:    %f\n", iers_rec->ut1_utc_a); // double
  printf("e_ut1_utc_a:  %f\n", iers_rec->e_ut1_utc_a); // double
  printf("lod_a:        %f\n", iers_rec->lod_a); // double
  printf("e_lod_a:      %f\n", iers_rec->e_lod_a); // double
  printf("nutflag_a:    %d\n", iers_rec->nutflag_a); // bool
  printf("dx_2000a_a:   %f\n", iers_rec->dx_2000a_a); // double
  printf("e_dx_2000a_a: %f\n", iers_rec->e_dx_2000a_a); // double
  printf("dy_2000a_a:   %f\n", iers_rec->dy_2000a_a); // double
  printf("e_dy_2000a_a: %f\n", iers_rec->e_dy_2000a_a); // double
  printf("pm_x_b:       %f\n", iers_rec->pm_x_b); // double
  printf("pm_y_b:       %f\n", iers_rec->pm_y_b); // double
  printf("ut1_utc_b:    %f\n", iers_rec->ut1_utc_b); // double
  printf("dx_2000a_b:   %f\n", iers_rec->dx_2000a_b); // double
  printf("dy_2000a_b:   %f\n", iers_rec->dy_2000a_b); // double
}

int main(int argc, const char * argv[]) {

  radiointerferometry_iers_record_t iers_rec = {0};
  iers_rec.mjd = 60709.5;
  int rv = radiointerferometry_iers_get(argv[1], &iers_rec);

  print(rv, &iers_rec);

  return rv;
}