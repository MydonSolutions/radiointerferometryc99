#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "radiointerferometryc99.h"

int main(int argc, const char * argv[]) {
  
  int test_value = calc_modified_julian_date_from_ymd(1999, 07, 10);
  printf("1999-07-10: %d MJD (expecting 51369)\n", test_value);
  if (test_value != 51369) {
    return 1;
  }
  test_value = calc_modified_julian_date_from_ymd(1899, 12, 31);
  printf("1899-12-31: %d MJD (expecting 15019)\n", test_value);
  if (test_value != 15019) {
    return 1;
  }
  test_value = calc_modified_julian_date_from_ymd(1900, 1, 1);
  printf("1900-01-01: %d MJD (expecting 15020)\n", test_value);
  if (test_value != 15020) {
    return 1;
  }
  test_value = calc_modified_julian_date_from_ymd(2000, 1, 1);
  printf("2000-01-01: %d MJD (expecting 51544)\n", test_value);
  if (test_value != 51544) {
    return 1;
  }
  test_value = calc_modified_julian_date_from_ymd(2026, 2, 16);
  printf("2026-02-16: %d MJD (expecting 61087)\n", test_value);
  if (test_value != 61087) {
    return 1;
  }
  
  return 0;
}