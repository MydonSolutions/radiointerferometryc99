/*
================================================================================
Byte-by-byte Description of each record (187 characters) in the file:
--------------------------------------------------------------------------------
   Bytes Format Units  Label  Explanations
--------------------------------------------------------------------------------
  1-   2   I2    ---     year         To get true calendar year, add 1900 for 
                                      MJD<=51543 or add 2000 for MJD>=51544)
  3-   4   I2    ---     month
  5-   6   I2    ---     day          of month
  8-  15   F8.2    d     MJD          fractional Modified Julian Date (MJD UTC)
      17   A1    ---     PolPMFlag_A  IERS (I) or Prediction (P) flag for 
                                      Bull. A polar motion values
 19-  27   F9.6  arcsec  PM_x_A       Bull. A PM-x
 28-  36   F9.6  arcsec  e_PM_x_A     error in PM-x (sec. of arc)
 38-  46   F9.6  arcsec  PM_y_A       Bull. A PM-y (sec. of arc)
 47-  55   F9.6  arcsec  e_PM_y_A     error in PM-y (sec. of arc)
      58   A1    ---     UT1Flag_A    IERS (I) or Prediction (P) flag for 
                                      Bull. A UT1-UTC values
 59-  68   F10.7 s       UT1_UTC_A    Bull. A UT1-UTC (sec. of time)
 69-  78   F10.7 s       e_UT1_UTC_A  error in UT1-UTC (sec. of time)
 80-  86   F7.4  ms      LOD_A        Bull. A LOD (msec. of time)
                                      -- NOT ALWAYS FILLED
 87-  93   F7.4  ms      e_LOD_A      error in LOD (msec. of time) 
                                      -- NOT ALWAYS FILLED
      96   A1    ---     NutFlag_A    IERS (I) or Prediction (P) flag for 
                                      Bull. A nutation values
 98- 106   F9.3  marcsec dX_2000A_A   Bull. A dX wrt IAU2000A Nutation 
                                      Free Core Nutation NOT Removed
107- 115   F9.3  marcsec e_dX_2000A_A error in dX (msec. of arc)
117- 125   F9.3  marcsec dY_2000A_A   Bull. A dY wrt IAU2000A Nutation
                                      Free Core Nutation NOT Removed
126- 134   F9.3  marcsec e_dY_2000A_A error in dY (msec. of arc)
135- 144   F10.6 arcsec  PM_X_B       Bull. B PM-x (sec. of arc)
145- 154   F10.6 arcsec  PM_Y_B       Bull. B PM-y (sec. of arc)
155- 165   F11.7 s       UT1_UTC_B    Bull. B UT1-UTC (sec. of time)
166- 175   F10.3 marcsec dX_2000A_B   Bull. B dX wrt IAU2000A Nutation
176- 185   F10.3 marcsec dY_2000A_B   Bull. B dY wrt IAU2000A Nutation
--------------------------------------------------------------------------------

********************************************************************************
An	for a character column made of n characters;
In	for a column containing an integer number of n digits;
Fn.d	for a column containing a number of width n digits and up to d digits in the fractional part;
En.d    for a number using the exponential notation
Dn.d	
********************************************************************************
*/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

typedef struct
{
  int year;            // Bytes   1-   2 ---     
  int month;           // Bytes   3-   4 ---     
  int day;             // Bytes   5-   6 ---     of month
  double mjd;          // Bytes   8-  15   d     fractional Modified Julian Date (MJD UTC)
  bool polpmflag_a;    // Bytes       17 ---     IERS (false) or Prediction (true) flag for Bull. A polar motion values
  double pm_x_a;       // Bytes  19-  27 arcsec  Bull. A PM-x
  double e_pm_x_a;     // Bytes  28-  36 arcsec  error in PM-x (sec. of arc)
  double pm_y_a;       // Bytes  38-  46 arcsec  Bull. A PM-y (sec. of arc)
  double e_pm_y_a;     // Bytes  47-  55 arcsec  error in PM-y (sec. of arc)
  bool ut1flag_a;      // Bytes       58 ---     IERS (false) or Prediction (true) flag for Bull. A UT1-UTC values
  double ut1_utc_a;    // Bytes  59-  68 s       Bull. A UT1-UTC (sec. of time)
  double e_ut1_utc_a;  // Bytes  69-  78 s       error in UT1-UTC (sec. of time)
  double lod_a;        // Bytes  80-  86 ms      Bull. A LOD (msec. of time) -- NOT ALWAYS FILLED
  double e_lod_a;      // Bytes  87-  93 ms      error in LOD (msec. of time) -- NOT ALWAYS FILLED
  bool nutflag_a;      // Bytes       96 ---     IERS (false) or Prediction (true) flag for Bull. A nutation values
  double dx_2000a_a;   // Bytes  98- 106 marcsec Bull. A dX wrt IAU2000A Nutation Free Core Nutation NOT Removed
  double e_dx_2000a_a; // Bytes 107- 115 marcsec error in dX (msec. of arc)
  double dy_2000a_a;   // Bytes 117- 125 marcsec Bull. A dY wrt IAU2000A Nutation Free Core Nutation NOT Removed
  double e_dy_2000a_a; // Bytes 126- 134 marcsec error in dY (msec. of arc)
  double pm_x_b;       // Bytes 135- 144 arcsec  Bull. B PM-x (sec. of arc)
  double pm_y_b;       // Bytes 145- 154 arcsec  Bull. B PM-y (sec. of arc)
  double ut1_utc_b;    // Bytes 155- 165 s       Bull. B UT1-UTC (sec. of time)
  double dx_2000a_b;   // Bytes 166- 175 marcsec Bull. B dX wrt IAU2000A Nutation
  double dy_2000a_b;   // Bytes 176- 185 marcsec Bull. B dY wrt IAU2000A Nutation
} radiointerferometry_iers_record_t;

/*
* The return argument `record` is expected to have its `mjd` field populated with the 
* date of the record to be read.
*
* Returns:
*  -1: error record->mjd < 41684.0, predating IERS records.
*  0: success
*  1: error opening filepath
*  2: error seeking to record, probably exceeding IERS records.
*  3: error reading record, should not occur.
*  4: error reading subsequent record, could not interpolate, MJD clipped, probably exceeding IERS records.
*  5: error interpolating, consequent record is not next day, sparse IERS records not expected.
*/
int radiointerferometry_iers_get(
  char* filepath,
  radiointerferometry_iers_record_t* record
);