import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from obspy.core import read, UTCDateTime
from scipy.interpolate import interp1d


stations = ["OBS01","OBS02","OBS03","OBS04",
            "OBS05","OBS06","OBS07","OBS08",
            "OBS09","OBS10","OBS11","OBS12",
            "NORTH","SOUTH"]

# all data of a given station must be in one folder
# e.g., /home/valeroe/data/network/station/*
data_path = "/data/valeroe/red_sea_obs/data_clock_corrected/ZF"
opath = "/data/valeroe/red_sea_obs/data_clock_corrected/ZF"
data_fmt = "mseed"

skew_file = "./obs_skew.txt"
estimated_drift_file = "./obs_estimated_clock_drift.txt"
show_figures = False

# EDIT WITH CAUTION
# -------------------
for station in stations:
    output_path = os.path.join(opath, station)

    if not os.path.isdir(output_path):
        print(f"Output directory {output_path} does not exist")
        raise Exception

    # load clock drift information
    if station in ["OBS04","NORTH","SOUTH"]:
        skew_info = np.loadtxt(estimated_drift_file, dtype="str", skiprows=(1))

        staidx = np.argwhere(skew_info[:, 0] == station)[0][0]

        begin_time = skew_info[staidx, 1]
        dps = skew_info[staidx, 2]

        begin_time = UTCDateTime(begin_time)
        dps = float(dps)

        print(station, begin_time, dps)
    else:
        skew_info = np.loadtxt(skew_file, dtype="str", skiprows=(1))

        staidx = np.argwhere(skew_info[:, 0] == station)[0][0]

        begin_time = skew_info[staidx, 1]
        end_time = skew_info[staidx, 2]
        skew = skew_info[staidx, 3]

        begin_time = UTCDateTime(begin_time)
        end_time = UTCDateTime(end_time)
        skew = float(skew)

        recording_time = end_time - begin_time  # in seconds
        dps = skew / recording_time  # clock drift per second

        print(station, begin_time, end_time, skew, dps)

    # fix clock drift
    sta_path = os.path.join(data_path, station)
    waveform_files = glob(os.path.join(sta_path, "*"))
    waveform_files.sort()

    for wfile in waveform_files:
        st = read(wfile, format=data_fmt)

        for tr in st:
            # step 1: get the GPS time at which the signal was actually recorded
            t_obs = tr.times(reftime=begin_time)  # seconds after begin time
            t_gps = t_obs/(1+dps)  # remove clock drift

            # step 2: interpolate the signal into the time axis of interest (OBS)
            ip1d_obj = interp1d(t_gps, tr.data, kind="linear", bounds_error=False,
                                fill_value="extrapolate")

            data2 = ip1d_obj(t_obs)

            if np.isnan(data2).any():
                raise Exception

            if show_figures:
                plt.plot(t_obs, tr.data, "r", alpha=0.7,
                         label="recorded wrong time")

                plt.plot(t_gps, tr.data, "k", alpha=0.7,
                         label="recorded correct time")

                plt.plot(t_obs, data2, "g--",alpha=0.7,
                         label="interpolated to times of interest")

                plt.legend()
                plt.show()
                plt.close()

            # overwrite data
            tr.data = data2.copy()

        # save file
        file_name = os.path.basename(wfile)
        st.write(os.path.join(output_path, file_name), format=data_fmt)
        print(f"{file_name} done")
