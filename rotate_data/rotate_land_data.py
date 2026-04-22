import matplotlib.pyplot as plt
import numpy as np
import os
from obspy.core import read, UTCDateTime, Stream
from glob import glob


def rotate_RT(cN, cE, az):
    """
        counter clock-wise axes rotation
        alpha in degrees
    """
    az = np.deg2rad(az)
    cR = np.cos(az) * cN + np.sin(az) * cE
    cT = -np.sin(az) * cN + np.cos(az) * cE
    return cR,cT


data_path = "/data/valeroe/red_sea_obs/data_clock_corrected"
output_path = "/data/valeroe/red_sea_obs/data_clock_corrected"

orientations_file = "./stations_orientations.txt"
dfmt = "mseed"
net = "ZF"

obs_stations = [f"OBS{x:02}" for x in range(1,13)]
obs_stations.extend(["NORTH", "SOUTH"])

# -------------------------------
# EDIT WITH CAUTION
# -------------------------------
stations_info = np.loadtxt(orientations_file, comments="#",
                           dtype=[("sta","U10"),("az",float),("unct",float)]
                          )

for i in range(stations_info.shape[0]):
    sta = stations_info[i]["sta"]

    if sta in obs_stations:
        continue

    sta_az = stations_info[i]["az"]
    print(sta, sta_az)

    stadir = os.path.join(data_path, net, sta)
    filesN = glob(os.path.join(data_path, net, sta, "*.?HN.*.mseed"))

    staoutdir = os.path.join(output_path, net, sta)
    if not os.path.isdir(staoutdir):
        os.makedirs(staoutdir)

    for fileN in filesN:
        fileE = fileN.replace("HN","HE")
        if not os.path.isfile(fileN) or not os.path.isfile(fileE):
            print("missing one component\n")
            continue
        else:
            stN = read(fileN, format="mseed")
            stE = read(fileE, format="mseed")

        # check both streams have same number of traces
        if len(stN) != len(stE):
            print("streams must have same number of traces\n")
            continue

        for j in range(0, len(stN)):
            # check they span the same time range
            if stN[j].stats.starttime != stE[j].stats.starttime:
                newstartime = max(stN[j].stats.starttime,
                                  stE[j].stats.starttime)

                print(stN[j])
                print(stE[j])
                print("components N and E have a different starttime")
                print("cutting traces to new starttime: ", newstartime,"\n")

                stN[j] = stN[j].slice(starttime=newstartime,
                                      endtime=stN[j].stats.endtime,
                                      nearest_sample=True)

                stE[j] = stE[j].slice(starttime=newstartime,
                                      endtime=stE[j].stats.endtime,
                                      nearest_sample=True)

            if stN[j].stats.endtime != stE[j].stats.endtime:
                newendtime = min(stN[j].stats.endtime,
                                 stE[j].stats.endtime)

                print(stN[j])
                print(stE[j])
                print("components N and E have a different endtime")
                print("cutting traces to new endtime: ", newendtime,"\n")

                stN[j] = stN[j].slice(starttime=stN[j].stats.starttime,
                                      endtime=newendtime,
                                      nearest_sample=True)

                stE[j] = stE[j].slice(starttime=stE[j].stats.starttime,
                                      endtime=newendtime,
                                      nearest_sample=True)

            cN, cE = rotate_RT(stN[j].data, stE[j].data, -sta_az)

            stN[j].data = cN.copy()
            stN[j].data = stN[j].data.astype("float32")
            stE[j].data = cE.copy()
            stE[j].data = stE[j].data.astype("float32")

        # write north data
        ofnameN = os.path.basename(fileN)
        stN.write(os.path.join(staoutdir, ofnameN), format=dfmt)

        # write east data
        ofnameE = os.path.basename(fileE)
        stE.write(os.path.join(staoutdir, ofnameE), format=dfmt)

    print(f'{sta} done')
