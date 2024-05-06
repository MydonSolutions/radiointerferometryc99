#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "radiointerferometryc99.h"

int main(int argc, const char * argv[]) {
  size_t count = 1;
  double time_jd[] = {2400000.5+41691.5};
  double app_ra_radians[] = {8.3*RADIOINTERFEROMETERY_PI/180};
  double app_dec_radians[] = {16.3*RADIOINTERFEROMETERY_PI/180};
  double pos_angle[] = {1};
  double latitude = 33.97391383157283*RADIOINTERFEROMETERY_PI/180.0;
  double longitude = -116.5833461618117*RADIOINTERFEROMETERY_PI/180.0;
  double altitude = 1073.4610445341686;
  double offset_pos = RADIOINTERFEROMETERY_PI/360.0;
  
  int rv = calc_itrs_icrs_frame_pos_angle(
    time_jd,
    app_ra_radians,
    app_dec_radians,
    count,
    longitude,
    latitude,
    altitude,
    offset_pos,
    argv[1],
    pos_angle
  );

  printf("rv: %d, count: %d\n------------------\n", rv, ((rv/10)-1)/2);
  for (size_t i = 0; i < count; i++)
  {
    printf("posangle %ld: %f\n", i, pos_angle[i]);
  }
  
  return rv;
}