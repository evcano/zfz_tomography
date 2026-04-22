import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from scipy.stats import zscore

#net = "ZF"
#stations_list = [f"OBS{x:02}" for x in range(1, 13)]
#stations_list.extend(['BREEM','KHUF','LAVA','NORTH','QUMAN','SOUTH'])

#net = "KL"
#stations_list = ["DEEP","EAST","WEST","SOUTH"]

net = "SA"
stations_list = ["EWJHS","WJHS"]

cmp = "N"
indir = "./bad_windows"
thr_zscore  = 3.

# --------------------------------------------------
# station loop
for sta in stations_list:
    sta = f"{net}.{sta}"

    infiles = glob(os.path.join(indir,f"{sta}.{cmp}_WAMP_*.dat"))
    sta_outliers = []

    # frequency band loop
    for file_ in infiles:
        band_outliers = []

        dates = np.loadtxt(file_,usecols=0, dtype= str)
        values = np.loadtxt(file_, usecols=1, dtype=float)

        # detect outliers based on percentile
        p = np.percentile(values, 95)
        oidx = np.where(values > p)[0].tolist()
        band_outliers.extend(dates[oidx])

#        plt.plot(values, "ok", alpha=0.5, markersize=1.0)
#        plt.plot(oidx, values[oidx], "or", alpha=0.5, markersize=1.0)
#        plt.show()
#        plt.close()

        # keep inliers 
        iidx = np.where(values <= p)[0].tolist()
        dates2 = dates[iidx]
        values2 = values[iidx]

        # detect outliers based on zscore
        z = zscore(values2)
        oidx = np.where(abs(z) > thr_zscore)[0]
        band_outliers.extend(dates2[oidx])

#        plt.plot(values2, "ok", alpha=0.5, markersize=1.0)
#        plt.plot(oidx, values2[oidx], "or", alpha=0.5, markersize=1.0)
#        plt.show()
#        plt.close()

        # pass band outliers to sta outliers
        sta_outliers.extend(band_outliers)

    # eliminate duplicated outliers
    sta_outliers = list(set(sta_outliers))
    sta_outliers.sort()

    # print info
    nwin = len(dates)
    nout = len(sta_outliers)
    print(f"station : {sta}")
    print(f"nbands: {len(infiles)}")
    print(f"outliers: {nout}/{nwin} ({nout/nwin}%)")

    outfile = f"{sta}.{cmp}_OUTLIERS.dat"
    with open(os.path.join(indir,outfile), "w") as file_:
        for out in sta_outliers:
            file_.write(f"{out}\n")
