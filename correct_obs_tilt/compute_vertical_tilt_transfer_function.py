import matplotlib.pyplot as plt
import numpy as np
import obspy as obs
import os
import yaml
from mpi4py import MPI
from glob import glob
from obstools.atacr import DayNoise, StaNoise
from obstools.atacr.utils import coherence, phase


# NOTES
# Z1: vertical minus horizontal 1
# Z2-1: vertical minus horizontal 1 and 2
# ZP-21: vertical minus horizontal 1 and 2 and pressure
# ZH: vertical minus rotated horizontal
# ZP-H: vertical minus rotated horizontal and pressure
# ZP: vertical minus pressure

# Bowden computes the transfer functions in periods of
# 12 hours and tapers them to zero outside the 5 to 15 band

# Tian sets a cutoff frequenty of 0.1 Hz to the transfer
# functions

def read_transient_win(file_):
    win_dates = np.loadtxt(file_, dtype=str)
    dates_utc = {}
    for w in win_dates:
        utc = obs.UTCDateTime(w)
        if utc.date.isoformat() not in dates_utc.keys():
            dates_utc[utc.date.isoformat()] = []
        dates_utc[utc.date.isoformat()].append(utc)
    return dates_utc

# -------------------------------------------------
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nproc = comm.Get_size()

# PARAMETERS
# ------------------
data_path = "/data/valeroe/red_sea_obs/data_processed/ZF"
psd_win_dur = 1800.0  # duration of window to compute PSD
wdur = 3600.0  # duration of windows with transient signals

# frequency band at which QC of daily spectra is performed
# default is 0.004 to 0.2
fqmin = 0.004
fqmax = 0.2

save = True
show_figure = False

stations_list = [f"OBS{x:02}" for x in range(1,13)]
stations_list.extend(["NORTH","SOUTH"])

# EDIT WITH CAUTION
# ------------------------------------------------
sta = stations_list[myrank]
print(sta)

outdir = f"./transfer_functions/{sta}"
if not os.path.isdir(outdir):
    os.mkdir(outdir)

cut_win_file = f"./bad_windows/ZF.{sta}_OUTLIERS.dat"
cut_win_dates = read_transient_win(cut_win_file)

wfiles = glob(os.path.join(data_path,sta,f"*.*HZ.*.mseed"))
wfiles.sort()

for fileZ in wfiles:
    # read data
    file1 = fileZ.replace("HZ", "H1")
    file2 = fileZ.replace("HZ", "H2")

    st = obs.Stream()
    try:
        st += obs.read(fileZ, format="mseed")
        st += obs.read(file1, format="mseed")
        st += obs.read(file2, format="mseed")
    except Exception:
        print("missing one component")
        continue

    st.sort(keys=["starttime"])
    sutc = st[-1].stats.starttime
    day = sutc.date.isoformat()

    # skip already done files
    if os.path.isfile(os.path.join(outdir, f"cZZ_{day}.npy")):
        print("skipping file")
        continue

    # cut windows with transient signals
    if day in cut_win_dates.keys():
        for win_start in cut_win_dates[day]:
            win_end = win_start + wdur
            st.cutout(starttime=win_start, endtime=win_end)

    # remove short segments and merge
    for tr in st:
        if (tr.stats.npts * tr.stats.delta) < psd_win_dur*2:
            st.remove(tr)

    if not st.select(component="Z"):
        print("missing component; skipping")
        continue
    elif not st.select(component="1"):
        print("missing component; skipping")
        continue
    elif not st.select(component="2"):
        print("missing component; skipping")
        continue

    st.merge()

    # make sure all components span the same time range
    st.sort(keys=["starttime"])
    sutc = st[-1].stats.starttime
    st.sort(keys=["endtime"])
    eutc = st[0].stats.endtime
    st = st.slice(starttime=sutc, endtime=eutc)

    # get components
    trZ = st.select(component="Z")[0]
    tr1 = st.select(component="1")[0]
    tr2 = st.select(component="2")[0]

    # the PSD quality-check is conducted in the specified frequency band.
    # the resulting cross/auto spectral densities are not filtered!!!
    # see QC_daily_spectra in classes.py
    try:
        daynoise = DayNoise(
            tr1=tr1,
            tr2=tr2,
            trZ=trZ,
            trP=obs.Trace(),
            window=psd_win_dur,
            overlap=0.25,
        )

        daynoise.QC_daily_spectra(
            pd=[fqmin, fqmax],
            smooth=True,
            tol=1.5,
            alpha=0.05,
            fig_QC=False,
        )

        daynoise.average_daily_spectra(
            calc_rotation=True,
            fig_average=False,
            fig_coh_ph=False,
        )
    except Exception:
        print("something went wrong")
        raise Exception

    # NOTE: the coherence at titlt angle is the mean of the coherence!
    # for some reason, the coherence is hard-coded to be compute between
    # 0.005 and 0.035 Hz; see lines 368 and 429 of utils.py
    tilt = daynoise.rotation.tilt
    coh = daynoise.rotation.coh_value
    faxis = daynoise.f

    cZZ = daynoise.power.cZZ
    cHH = daynoise.rotation.cHH
    cHZ = daynoise.rotation.cHZ

    if show_figure:
        Co = np.abs(cHZ)**2 / (cHH*cZZ)
        Ad = np.abs(cHZ) / cHH
        Ph = np.angle(cHZ/cHH)  # -pi to pi radians 

        idx = (faxis >= 0.0) & (faxis <= 0.2)
        fig,ax = plt.subplots(3)
        ax[0].plot(faxis[idx], Co[idx])
        ax[1].plot(faxis[idx], Ad[idx])
        ax[2].plot(faxis[idx], Ph[idx])
        plt.show()
        plt.close()

    if save:
        np.save(os.path.join(outdir,f"cZZ_{day}.npy"),cZZ)
        np.save(os.path.join(outdir,f"cHH_{day}.npy"),cHH)
        np.save(os.path.join(outdir,f"cHZ_{day}.npy"),cHZ)
        np.save(os.path.join(outdir,f"tilt_coh_{day}.npy"),np.array([tilt, coh]))
        np.save("freqaxis.npy", faxis)
        print(f"done {fileZ}")
