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


data_path = "/data/valeroe/red_sea_obs/data_tilt_clock_corrected"
output_path = "/data/valeroe/red_sea_obs/data_tilt_clock_corrected"

orientations_file = "./stations_orientations.txt"
dfmt = "mseed"
net = "ZF"

# -------------------------------
# EDIT WITH CAUTION
# -------------------------------
stations_info = np.loadtxt(orientations_file, comments="#",
                           dtype=[("sta","U10"),("az",float),("unct",float)]
                          )

for i in range(stations_info.shape[0]):
    sta = stations_info[i]["sta"]
    sta_az = stations_info[i]["az"]
    print(sta, sta_az)

    stadir = os.path.join(data_path, net, sta)
    files1 = glob(os.path.join(data_path, net, sta, "*.?H1.*.mseed"))

    staoutdir = os.path.join(output_path, net, sta)
    if not os.path.isdir(staoutdir):
        os.makedirs(staoutdir)

    for file1 in files1:
        file2 = file1.replace("H1","H2")
        if not os.path.isfile(file1) or not os.path.isfile(file2):
            print("missing one component")
            continue
        else:
            st1 = read(file1, format="mseed")
            st2 = read(file2, format="mseed")

        # check there is only one trace per stream
        if len(st1) != 1 or len(st2) != 1:
            print("streams must contain onlye one trace")
            raise Exception
        else:
            tr1 = st1[0]
            tr2 = st2[0]

        # check they span the same time range
        if tr1.stats.starttime != tr2.stats.starttime:
            print("components 1 and 2 have a different starttime")
            raise Exception

        if tr1.stats.endtime != tr2.stats.endtime:
            print("components 1 and 2 have a different endtime")
            raise Exception

        # For FUGRO OBS
        # This script rotates h1 to the north and h2 to the east
        # assuming stations follow a left-hand convention.
        # However, FUGRO obs follow right-hand convention, meaning that
        # h2 is -90 degrees from h1.
        # Here, we change polarization so h2 is 90 degrees from h1 as in
        # the left-hand convention.
        if sta == "NORTH" or sta == "SOUTH":
            tr2.data = -tr2.data

        # rotate to north-east given a coordinate system with
        # h1 and h2 as main axes
        cN, cE = rotate_RT(tr1.data, tr2.data, -sta_az)

        # write north data
        trN = tr1.copy()
        trN.data = cN.copy()
        trN.data = trN.data.astype("float32")
        trN.stats.channel = trN.stats.channel.replace("1","N")

        ofnameN = os.path.basename(file1)
        ofnameN = ofnameN.replace("H1","HN")

        trN.write(os.path.join(staoutdir, ofnameN), format=dfmt)

        # write east data
        trE = tr2.copy()
        trE.data = cE.copy()
        trE.data = trE.data.astype("float32")
        trE.stats.channel = trE.stats.channel.replace("2","E")

        ofnameE = os.path.basename(file2)
        ofnameE = ofnameE.replace("H2","HE")

        trE.write(os.path.join(staoutdir, ofnameE), format=dfmt)

    print(f'{sta} done')
