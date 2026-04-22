import matplotlib.pyplot as plt
import numpy as np
import obspy as obs
import os
import yaml
from glob import glob
from mpi4py import MPI


comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nproc = comm.Get_size()


#data_path = "/data/valeroe/red_sea_obs/data_tilt_clock_corrected"
#net = "ZF"
#cha = "?HZ"
#stations_list = [f"OBS{x:02}" for x in range(1, 13)]
#stations_list.extend(['BREEM','KHUF','LAVA','NORTH','QUMAN','SOUTH'])

#data_path = "/data/valeroe/harrat_lunayir/data_processed"
#net = "KL"
#cha = "H?E"
#stations_list = ["DEEP", "EAST", "WEST", "SOUTH"]

data_path = "/data/valeroe/sgs_for_zafran/data_processed"
net = "SA"
cha = "?HN"
stations_list = ["EWJHS","WJHS"]

# -------------------
fmt = "mseed"
outdir = "./bad_windows"

windur = 3600.0
overlap = 1.0
fqbands = [(1/50.,1/25.), (1/25.,1/12.5), (1/12.5,1/6.25),  (1/6.25,1.0)]


# EDIT WITH CAUTION
# --------------------
sta = stations_list[myrank]
print(myrank,sta)

waveform_files = glob(os.path.join(data_path,net,sta,f"{net}.{sta}.{cha}*{fmt}"))
waveform_files.sort()

win_start = {}
win_amp = {}

for band in fqbands:
    bcode = f"{band[0]}-{band[1]}"
    win_start[bcode] = []
    win_amp[bcode] = []

for file_ in waveform_files:
    st = obs.read(file_, format=fmt)

    st.detrend("linear")
    st.detrend("demean")
    st.taper(0.05,type="hann")

    for band in fqbands:
        bcode = f"{band[0]}-{band[1]}"

        st2 = st.copy()
        if band[0] == band[1]:
            st2.filter("highpass",freq=band[0],
                       corners=2, zerophase=True)
        else:
            st2.filter("bandpass",freqmin=band[0],freqmax=band[1],
                       corners=2, zerophase=True)

        for tr in st2:
            for tr_win in tr.slide(windur, windur-overlap):
                starttime = tr_win.stats.starttime.format_iris_web_service()
                amp = np.max(np.abs(tr_win.data))

                win_start[bcode].append(starttime)
                win_amp[bcode].append(amp)

for band in fqbands:
    bcode = f"{band[0]}-{band[1]}"
    outfile = f"{net}.{sta}.{cha[-1]}_WAMP_{bcode}.dat"

    with open(os.path.join(outdir, outfile), "w") as file_:
        for i in range(0, len(win_start[bcode])):
            file_.write(f"{win_start[bcode][i]} {win_amp[bcode][i]}\n")

print(myrank,"done")
