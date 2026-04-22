# This code was used to stack noise correlations for the inversion

import matplotlib.pyplot as plt
import numpy as np
import os
from mpi4py import MPI
from obspy import read, Stream, Trace
from sanpy.base.project_functions import load_project
from sanpy.base.functions import distribute_objects
from stockwell import st as swell


def tfpws_stack(data, dt, f1, f2, v):
    ntr = data.shape[0]
    npts = data.shape[1]
    t = np.arange(npts) * dt
    df = 1.0 / (t[-1]-t[0])

    f1samp = int(f1/df)
    f2samp = int(f2/df)

    for i in range(0, ntr):
        d = data[i,:]
        s = swell.st(d, f1samp, f2samp)

        if i == 0:
            f = np.linspace(f1, f2, len(s[:,0]))
            T, F = np.meshgrid(t, f)
            c = np.zeros_like(s)

        ph = s / np.abs(s) * np.exp(2j*np.pi*F*T)
        c += ph

    c /= ntr

    stk = np.mean(data, axis=0)
    ds = swell.st(stk, f1samp, f2samp)

    stk_tfpws = ds * np.abs(c) ** v
    stk_tfpws = swell.ist(stk_tfpws, f1samp, f2samp)

    return stk_tfpws


def plot_corr(st, tr_stk, lags, good_tr_idx=None, swin=None, nwin=None,
              pair=None):

    ntraces = len(st)

    st.normalize()
    tr_stk.normalize()

    sd = st[0].stats.starttime
    ed = st[-1].stats.endtime
    recdays = np.arange(sd, ed+86400.0, dtype='datetime64[D]')

    fig, ax = plt.subplots(2, figsize=(8,10), height_ratios=(5,1))

    for i in range(0, ntraces):
        if not good_tr_idx:
            c = "k"
            a = 0.5
        elif i in good_tr_idx:
            c = "k"
            a = 0.5
        else:
            c = "r"
            a = 0.1

        tr = st[i]
        j = tr.stats.starttime.date - sd.date
        x = tr.data*2.0 + j.days
        ax[0].plot(lags, x, c=c, linewidth=1.0, alpha=a)

    ax[1].plot(lags, tr_stk.data, c="k")

    st.stack()
    st.normalize()
    ax[1].plot(lags,st[0].data,c="b")

    if swin:
        for i in range(0, 2):
            ax[i].axvline(x=-swin[0],c="b")
            ax[i].axvline(x=-swin[1],c="b")
            ax[i].axvline(x=swin[0],c="b")
            ax[i].axvline(x=swin[1],c="b")

    if nwin:
        for i in range(0, 2):
            ax[i].axvline(x=-nwin[0],c="r")
            ax[i].axvline(x=-nwin[1],c="r")
            ax[i].axvline(x=nwin[0],c="r")
            ax[i].axvline(x=nwin[1],c="r")

    if pair:
        ax[0].set_title(pair)

    plt.show()
    plt.close()
    return None


comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nproc = comm.Get_size()

# general parameters
cmp = "RR"

project_file = \
"/data/valeroe/red_sea_obs_2/noise_correlations_tomo/zf_correlation.pkl"

data_path = \
f"/data/valeroe/red_sea_obs_2/noise_correlations_tomo/daily_corr/{cmp}"

output_path = \
f"/data/valeroe/red_sea_obs_2/noise_correlations_tomo/stack_tfpws/all/{cmp}"

save_stk = True


#  EDIT WITH CAUTION
# =======================================================
# sanpy project
p = load_project(project_file)
pairs = p.pairs_list
pairs.sort()

pairs_per_core = distribute_objects(pairs, nproc, myrank)

for pair in pairs_per_core:
    s1, s2 = pair.split("_")

    # skip autocorrelation
    if s1 == s2:
        continue

    # read waveforms
    ddir = os.path.join(data_path, pair, "*")

    try:
        st_all = read(ddir, format="sac")
    except:
        print(f"error reading data of {pair}")
        continue

    st_all.sort(keys=["starttime"])

    dt = st_all[0].stats.delta
    nt = st_all[0].stats.npts

    maxtime = ((nt - 1) / 2) * dt
    lags = np.linspace(-maxtime, maxtime, nt)

    # filter waveforms
    fmin = 1./50
    fmax = 1.0

    st_all.detrend("demean")
    st_all.detrend("linear")
    st_all.taper(0.05)
    st_all.filter("bandpass", freqmin=fmin, freqmax=fmax,
                  corners=2, zerophase=True)

    # tfpws stacking
    data = [tr.data for tr in st_all]
    data = np.array(data)
    stacked_corr = tfpws_stack(data, dt, fmin, fmax, 2)

    ntraces = data.shape[0]

    # write
    if save_stk:
        tr_stk = Trace(data=stacked_corr, header=st_all[0].stats)
        tr_stk.stats.sac.user7 = ntraces
        tr_stk.write(os.path.join(output_path,f"{pair}_{cmp}.sac"),format="sac")

    print("done", pair, myrank)

print(f"core {myrank} done")
