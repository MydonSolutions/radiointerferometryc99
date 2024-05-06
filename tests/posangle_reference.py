from pyuvdata.utils import calc_frame_pos_angle as pyuvcalc_frame_pos_angle, np

print(
    pyuvcalc_frame_pos_angle(
        np.array([2400000.5+41691.5]), # time_jd
        np.array([8.3*np.pi/180]), # app_ra_radians
        np.array([16.3*np.pi/180]), # app_dec_radians
        ( # telescope_loc
            33.97391383157283*np.pi/180, # latitude
            -116.5833461618117*np.pi/180, # longitude
            1073.4610445341686, # altitude
        ),
        "icrs", # ref_frame
        offset_pos = np.pi/360.0, #
    )
)
