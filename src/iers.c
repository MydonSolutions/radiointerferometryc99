#include <radiointerferometryc99/iers.h>

int _iers_record_parse(
  char* char_record,
  radiointerferometry_iers_record_t* record
) {
  char* char_record_end;
  
  char_record_end = char_record+2;
  record->year = strtol(
    char_record+1-1,
    &char_record_end,
    10
  );
  
  char_record_end = char_record+4;
  record->month = strtol(
    char_record+3-1,
    &char_record_end,
    10
  );
  
  char_record_end = char_record+6;
  record->day = strtol(
    char_record+5-1,
    &char_record_end,
    10
  );

  char_record_end = char_record+15;
  record->mjd = strtod(
    char_record+8-1,
    &char_record_end
  );
  record->year += record->mjd < 51544.0 ? 1900 : 2000;

  record->polpmflag_a = char_record[17-1] != 'I';

  char_record_end = char_record+27;
  record->pm_x_a = strtod(
    char_record+19-1,
    &char_record_end
  );

  char_record_end = char_record+36;
  record->e_pm_x_a = strtod(
    char_record+28-1,
    &char_record_end
  );

  char_record_end = char_record+46;
  record->pm_y_a = strtod(
    char_record+38-1,
    &char_record_end
  );

  char_record_end = char_record+55;
  record->e_pm_y_a = strtod(
    char_record+47-1,
    &char_record_end
  );
  record->ut1flag_a = char_record[58-1] != 'I';

  char_record_end = char_record+68;
  record->ut1_utc_a = strtod(
    char_record+59-1,
    &char_record_end
  );

  char_record_end = char_record+78;
  record->e_ut1_utc_a = strtod(
    char_record+69-1,
    &char_record_end
  );

  char_record_end = char_record+86;
  record->lod_a = strtod(
    char_record+80-1,
    &char_record_end
  );

  char_record_end = char_record+93;
  record->e_lod_a = strtod(
    char_record+87-1,
    &char_record_end
  );
  record->nutflag_a = char_record[96-1] != 'I';

  char_record_end = char_record+106;
  record->dx_2000a_a = strtod(
    char_record+98-1,
    &char_record_end
  );

  char_record_end = char_record+115;
  record->e_dx_2000a_a = strtod(
    char_record+107-1,
    &char_record_end
  );

  char_record_end = char_record+125;
  record->dy_2000a_a = strtod(
    char_record+117-1,
    &char_record_end
  );

  char_record_end = char_record+134;
  record->e_dy_2000a_a = strtod(
    char_record+126-1,
    &char_record_end
  );

  char_record_end = char_record+144;
  record->pm_x_b = strtod(
    char_record+135-1,
    &char_record_end
  );

  char_record_end = char_record+154;
  record->pm_y_b = strtod(
    char_record+145-1,
    &char_record_end
  );

  char_record_end = char_record+165;
  record->ut1_utc_b = strtod(
    char_record+155-1,
    &char_record_end
  );

  char_record_end = char_record+175;
  record->dx_2000a_b = strtod(
    char_record+166-1,
    &char_record_end
  );

  char_record_end = char_record+185;
  record->dy_2000a_b = strtod(
    char_record+176-1,
    &char_record_end
  );
  return 0;
}

int radiointerferometry_iers_get(
  const char* filepath,
  radiointerferometry_iers_record_t* record
) {
  if (record->mjd < 41684.0) {
    return -1;
  }

  int fd = open(filepath, O_RDONLY);
  if(fd < 0) {
    return 1;
  }

  off_t file_size = lseek(fd, 0, SEEK_END);
  int record_index = record->mjd - 41684.0;

  // seek as close to record as naively can
  if ((record_index+2)*188 > file_size) {
    lseek(fd, ((file_size/188)-1)*188, SEEK_SET);
  }
  else {
    lseek(fd, record_index*188, SEEK_SET);
  }

  // search for record
  double mjd = record->mjd;
  char char_record[187+2+188] = {0};
  char* char_record_mjd_end = char_record+15;
  read(
    fd,
    char_record,
    15
  );
  record->mjd = strtod(
    char_record+8-1,
    &char_record_mjd_end
  );
  int search_direction = record->mjd == mjd ? 0 : (record->mjd < mjd ? 1 : -1);

  // search in one direction
  while (search_direction != 0) {
    if (search_direction == -1) {
      if (record->mjd <= mjd) {
        break;
      }
      lseek(fd, -15-188, SEEK_CUR);
    }
    else {
      if (record->mjd >= mjd) {
        break;
      }
      lseek(fd, (188-15), SEEK_CUR);
    }
    
    if (15 > read(fd, char_record, 15)) {
      close(fd);
      return 2;
    }
    record->mjd = strtod(
      char_record+8-1,
      &char_record_mjd_end
    );
  }

  // read the rest of the record and the following (for interpolation)
  int bytes_read = read(
    fd,
    char_record+15,
    (188*2)-15
  );
  close(fd);
  if (188 > bytes_read) {
    return 3;
  }

  _iers_record_parse(
    char_record,
    record
  );
  if (record->mjd == mjd) {
    // no need to interpolate
    return 0;
  }
  if ((2*188)-15 > bytes_read) {
    // could not read the subsequent record,
    // cannot interperolate
    return 4;
  }
  
  radiointerferometry_iers_record_t next_record;
  _iers_record_parse(
    char_record+188,
    &next_record
  );

  if ((int)(next_record.mjd - record->mjd) > 1) {
    return 5;
  }

  // linearly interpolate between records
  double fraction = mjd - record->mjd;
  record->pm_x_a        += fraction*(next_record.pm_x_a       - record->pm_x_a);
  record->e_pm_x_a      += fraction*(next_record.e_pm_x_a     - record->e_pm_x_a);
  record->pm_y_a        += fraction*(next_record.pm_y_a       - record->pm_y_a);
  record->e_pm_y_a      += fraction*(next_record.e_pm_y_a     - record->e_pm_y_a);
  record->ut1_utc_a     += fraction*(next_record.ut1_utc_a    - record->ut1_utc_a);
  record->e_ut1_utc_a   += fraction*(next_record.e_ut1_utc_a  - record->e_ut1_utc_a);
  record->lod_a         += fraction*(next_record.lod_a        - record->lod_a);
  record->e_lod_a       += fraction*(next_record.e_lod_a      - record->e_lod_a);
  record->dx_2000a_a    += fraction*(next_record.dx_2000a_a   - record->dx_2000a_a);
  record->e_dx_2000a_a  += fraction*(next_record.e_dx_2000a_a - record->e_dx_2000a_a);
  record->dy_2000a_a    += fraction*(next_record.dy_2000a_a   - record->dy_2000a_a);
  record->e_dy_2000a_a  += fraction*(next_record.e_dy_2000a_a - record->e_dy_2000a_a);
  record->pm_x_b        += fraction*(next_record.pm_x_b       - record->pm_x_b);
  record->pm_y_b        += fraction*(next_record.pm_y_b       - record->pm_y_b);
  record->ut1_utc_b     += fraction*(next_record.ut1_utc_b    - record->ut1_utc_b);
  record->dx_2000a_b    += fraction*(next_record.dx_2000a_b   - record->dx_2000a_b);
  record->dy_2000a_b    += fraction*(next_record.dy_2000a_b   - record->dy_2000a_b);
  record->mjd += fraction;

  return 0;
}