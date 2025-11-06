import ctypes
import sys
assert len(sys.argv) == 3, "Provide the RadioInterferometryC99 library .so filepath and the IERS filepath."

lib_so_path, iers_path = sys.argv[1:3]
print(lib_so_path, iers_path)

from astropy.utils import iers
iers.IERS_A_FILE = iers_path
from pyuvdata.utils.phasing import calc_frame_pos_angle as pyuvcalc_frame_pos_angle, np

libri = ctypes.CDLL(lib_so_path)
libri.calc_itrs_icrs_frame_pos_angle.argtypes = (
  ctypes.POINTER(ctypes.c_double), #time_jd,
  ctypes.POINTER(ctypes.c_double), #app_ra_radians,
  ctypes.POINTER(ctypes.c_double), #app_dec_radians,
  ctypes.c_size_t, #count,
  ctypes.c_double, #longitude_rad,
  ctypes.c_double, #latitude_rad,
  ctypes.c_double, #altitude,
  ctypes.c_double, #offset_pos,
  ctypes.c_char_p, #iers_filepath,
  ctypes.POINTER(ctypes.c_double), #pos_angle
)
libri.calc_itrs_icrs_frame_pos_angle.restypes = ctypes.c_int

def calc_itrs_icrs_frame_pos_angle(
    time_jd: list[float],
    app_ra_radians: list[float],
    app_dec_radians: list[float],
    longitude_rad: float,
    latitude_rad: float,
    altitude: float,
    iers_filepath: str,
    offset_pos: float = np.pi/360.0
) -> list[float]:
    assert len(time_jd) == len(app_ra_radians)
    assert len(app_ra_radians) == len(app_dec_radians)
    count = len(time_jd)

    ret = (ctypes.c_double*count)(*([0.0]*count))
    time_jd_c = (ctypes.c_double*count)(*time_jd)
    app_ra_radians_c = (ctypes.c_double*count)(*app_ra_radians)
    app_dec_radians_c = (ctypes.c_double*count)(*app_dec_radians)
    rcode = libri.calc_itrs_icrs_frame_pos_angle(
        time_jd_c,
        app_ra_radians_c,
        app_dec_radians_c,
        count,
        longitude_rad,
        latitude_rad,
        altitude,
        offset_pos,
        iers_filepath.encode(),
        ret
    )
    assert rcode == 0, f"Non-zero Return Code: {rcode}"
    return [r for r in ret]

def pyuvcalc_itrs_icrs_frame_pos_angle(
    time_jd: list[float],
    app_ra_radians: list[float],
    app_dec_radians: list[float],
    longitude_rad: float,
    latitude_rad: float,
    altitude: float,
    iers_filepath: str,
    offset_pos: float = np.pi/360.0
):  
    return pyuvcalc_frame_pos_angle(
        time_array = np.array(time_jd), # time_jd
        app_ra = np.array(app_ra_radians), # app_ra_radians
        app_dec = np.array(app_dec_radians), # app_dec_radians
        telescope_loc = ( # telescope_loc
            longitude_rad,
            latitude_rad,
            altitude
        ),
        ref_frame = "icrs", # ref_frame
        offset_pos = offset_pos #
    )

static_test_cases = [
    (
        [2400000.5+41691.5], # time_jd
        [8.3*np.pi/180], # app_ra_radians
        [16.3*np.pi/180], # app_dec_radians
        33.97391383157283*np.pi/180, # latitude
        -116.5833461618117*np.pi/180, # longitude
        1073.4610445341686, # altitude
    ),
    (
        [2400000.5+41700.00], # time_jd
        [60.3*np.pi/180], # app_ra_radians
        [37.5555*np.pi/180], # app_dec_radians
        33.97391383157283*np.pi/180, # latitude
        -116.5833461618117*np.pi/180, # longitude
        1073.4610445341686, # altitude
    ),
    (
        [2400000.5+41708.00], # time_jd
        [60.3*np.pi/180], # app_ra_radians
        [37.5555*np.pi/180], # app_dec_radians
        33.97391383157283*np.pi/180, # latitude
        -116.5833461618117*np.pi/180, # longitude
        1073.4610445341686, # altitude
    ),
    (
        [2400000.5+60706.00], # time_jd
        [60.3*np.pi/180], # app_ra_radians
        [37.5555*np.pi/180], # app_dec_radians
        33.97391383157283*np.pi/180, # latitude
        -116.5833461618117*np.pi/180, # longitude
        1073.4610445341686, # altitude
    )
]

for test_index, test_case in enumerate(static_test_cases):
    time_jd, app_ra_radians, app_dec_radians, latitude, longitude, altitude = test_case
    ri_ret = calc_itrs_icrs_frame_pos_angle(
        time_jd,
        app_ra_radians,
        app_dec_radians,
        latitude,
        longitude,
        altitude,
        iers_path
    )
    uv_ret = pyuvcalc_itrs_icrs_frame_pos_angle(
        time_jd,
        app_ra_radians,
        app_dec_radians,
        latitude,
        longitude,
        altitude,
        iers_path
    )
    print(f"ri: {ri_ret}")
    print(f"uv: {uv_ret}")
    assert np.isclose(ri_ret, uv_ret, rtol=1e-3), f"Test #{test_index}"
