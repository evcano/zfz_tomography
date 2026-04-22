import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import yaml
from scipy.signal import correlate, correlation_lags, tukey
from obspy import read, Stream, UTCDateTime
from obspy.signal.filter import envelope
from sanpy.base.project_functions import load_project


def compute_snr(wf, lags, swin, nwin):
    swin_s = np.argmin(abs(lags - swin[0]))
    swin_e = np.argmin(abs(lags - swin[1]))

    nwin_s = np.argmin(abs(lags - nwin[0]))
    nwin_e = np.argmin(abs(lags - nwin[1]))

    sig = wf[swin_s:swin_e].copy()
    noi = wf[nwin_s:nwin_e].copy()

    a = np.max(np.abs(sig))
    b = np.sqrt(np.mean(np.square(noi)))
    snr = a / b

    if snr < 0:
        snr = -10.0 * np.log10(-snr)
    elif snr > 0:
        snr = 10.0 * np.log10(snr)

    return snr


def determine_good_correlations(st, lags, wsig, wnoi, thr):
    wsig_pos = wsig
    wnoi_pos = wnoi

    wsig_neg = [-wsig[1], -wsig[0]]
    wnoi_neg = [-wnoi[1], -wnoi[0]]

    all_snr_pos = []
    all_snr_neg = []

    for tr in st:
        snr_pos = compute_snr(tr.data, lags, wsig_pos, wnoi_pos)
        snr_neg = compute_snr(tr.data, lags, wsig_neg, wnoi_neg)

        all_snr_pos.append(snr_pos)
        all_snr_neg.append(snr_neg)

    mean_snr_pos = np.mean(np.array(all_snr_pos)) * thr
    mean_snr_neg = np.mean(np.array(all_snr_neg)) * thr

    for tr in st:
        snr_pos = compute_snr(tr.data, lags, wsig_pos, wnoi_pos)
        snr_neg = compute_snr(tr.data, lags, wsig_neg, wnoi_neg)

        if snr_pos < mean_snr_pos or snr_neg < mean_snr_neg:
            st.remove(tr)

    return st

def determine_windows(lags, pair_dist, vmin, vmax):
    # on causal branch
    t1 = pair_dist / vmax
    t2 = pair_dist / vmin
    wsig = [t1, t2]
    wnoi = [t2*1.05, min(t2+t2-t1,lags[-2])]

    return wsig, wnoi


def estimate_clock_drift(daysno, tshifts, ccoefs):
    """ estimate linear clock drift using weighted least squares """
    G = np.ones((daysno.size, 2))
    G[:,0] = daysno[:]
    W = np.diag(ccoefs)
    d = tshifts

    Gw = np.matmul(W,G)
    dw = np.matmul(W,d)
    GtWG = np.matmul(G.T, Gw)
    GtWG_inv = np.linalg.inv(GtWG)

    m = np.dot(np.matmul(GtWG_inv, G.T), dw)
    m1 = m[0]  # clock drift per day
    m2 = m[1]

    return m1, m2


def plot_corr(st, tr_ref, lags, swin=None, nwin=None,
              mark_date=None, savefig=None):
    fig, ax = plt.subplots(1, figsize=(5,10))

    st.sort(keys=["starttime"])
    st.normalize()
    tr_ref.normalize()

    sd = st[0].stats.starttime
    ed = st[-1].stats.endtime
    recdays = np.arange(sd, ed+86400.0, dtype='datetime64[D]')

    data = np.zeros((recdays.size, st[0].stats.npts))

    for tr in st:
        idx = np.argwhere(recdays == tr.stats.starttime.date)
        data[idx, :] = tr.data

    ax.imshow(data, cmap="seismic", origin="lower",
              extent=(lags[0], lags[-1], recdays[0], recdays[-1]),
              aspect="auto", interpolation=None)

    if markdates:
        for date in markdates:
            ax.axhline(y=date.date, c="k")

    if swin:
        ax.axvline(x=-swin[0],c="b")
        ax.axvline(x=-swin[1],c="b")

        ax.axvline(x=swin[0],c="b")
        ax.axvline(x=swin[1],c="b")

    if nwin:
        ax.axvline(x=-nwin[0],c="r")
        ax.axvline(x=-nwin[1],c="r")

        ax.axvline(x=nwin[0],c="r")
        ax.axvline(x=nwin[1],c="r")

    if savefig:
        plt.savefig(savefig)
    else:
        plt.show()
    plt.close()

    return


def xcor_shift(s1, s2, dt, nt, maxlag=4):
    xcor = correlate(s2, s1, mode="full")
    xcor /= (np.linalg.norm(s1, ord=2) * np.linalg.norm(s2, ord=2))

    lags = correlation_lags(s2.size, s1.size, mode="full") * dt
    idx = np.where(np.abs(lags) > maxlag)
    xcor[idx] = 0.0

    ts = (np.argmax(xcor) - nt + 1) * dt  # tr > 0 means s2 is ahead s1
    cc = np.max(xcor)

    return ts, cc


#  EDIT WITH CAUTION
# =======================================================
# read parameters
parfile = sys.argv[1]

with open(parfile, "r") as _file:
    pars = yaml.safe_load(_file)

obs_sta = pars["obs_sta"]
pair = pars["pair"]
fqmin = pars["fqmin"]
fqmax = pars["fqmax"]
wsig = pars["wsig"]
wnoi = pars["wnoi"]
apply_window = pars["apply_window"]
vmin = pars["vmin"]
vmax = pars["vmax"]
markdates = pars["markdates"]
days2stack = pars["days2stack"]
data_path = pars["data_path"]
project_path = pars["project_path"]
skew_path = pars["skew_path"]
figpath = pars["figpath"]

if markdates:
    markdates = [UTCDateTime(x) for x in markdates]

# default matplotlib epoch
reference_date = UTCDateTime("1970-01-01T00:00")

# read obs skew info
if obs_sta:
    skew_info = np.loadtxt(skew_path, dtype="str", skiprows=(1))
    staidx = np.argwhere(skew_info[:, 0] == obs_sta)[0][0]

    deployment_time = skew_info[staidx, 1]
    recover_time = skew_info[staidx, 2]
    skew = skew_info[staidx, 3]

    deployment_time = UTCDateTime(deployment_time)
    recover_time = UTCDateTime(recover_time)
    recording_time = recover_time - deployment_time  # in seconds

    skew = float(skew)
    dpd = skew / recording_time  # clock drift per second
    dpd = dpd * 86400.0  # clock drift per day

    print(f"{obs_sta} drift per day {dpd} and skew {skew}")

# read waveforms
ddir = os.path.join(data_path, pair, "*")

st_all = read(ddir, format="sac")
st_all.sort(keys=["starttime"])
st_all = st_all[1:]  # skip correlation of first day

st_all.detrend("demean")
st_all.taper(0.05)
st_all.filter("bandpass",freqmin=fqmin,freqmax=fqmax,corners=2,zerophase=True)

dt = st_all[0].stats.delta
nt = st_all[0].stats.npts

maxtime = ((nt - 1) / 2) * dt
lags = np.linspace(-maxtime, maxtime, nt)

# determine signal and noise windows
if not wsig and not wnoi and apply_window:
    p = load_project(project_path)
    pair_dist = p.pairs[pair]["dis"] * 1000.0  # in meters
    wsig, wnoi = determine_windows(lags, pair_dist, vmin, vmax)

# define taper around arrival
if apply_window:
    taper = np.zeros(nt)
    i1 = np.argmin(np.abs(lags-wsig[0]))
    i2 = np.argmin(np.abs(lags-wsig[1]))

    taper[i1:i2] = tukey(i2-i1, alpha=0.1)  # causal branch

    if apply_window == "both":
        taper += taper[::-1]
    if apply_window == "acausal":
        taper = taper[::-1]
else:
    taper = np.ones(nt)

# reference correlation
maxstacktime = st_all[0].stats.starttime + 86400.0 * days2stack
st_ref = Stream()

for tr in st_all:
    if tr.stats.starttime < maxstacktime:
        st_ref += tr
    else:
        break

print("no of waveforms ", len(st_all))
print("nstack on st ref ", len(st_ref))
st_ref.stack()
tr_ref = st_ref[0]

# correlate
ccoefs = []
tshifts = []
days = []
daysno = []

for tr_day in st_all:
    ts, cc = xcor_shift(tr_ref.data * taper,
                        tr_day.data * taper,
                        dt,
                        nt
                       )

    ccoefs.append(cc)
    tshifts.append(ts)
    days.append(tr_day.stats.starttime.date)

    dn = tr_day.stats.starttime.date - reference_date.date
    daysno.append(dn.days)

ccoefs = np.array(ccoefs)
tshifts = np.array(tshifts)
days = np.array(days)
daysno = np.array(daysno)

# eliminate outliers
idx = np.argwhere((tshifts > -2.0) & (tshifts < 2.0))
ccoefs2 = ccoefs[idx].flatten()
tshifts2 = tshifts[idx].flatten()
days2 = days[idx].flatten()
daysno2 = daysno[idx].flatten()

# FIGURES
#-----------------------------
xlim1 = UTCDateTime("2021-11-01").date
xlim2 = UTCDateTime("2023-01-01").date

fig, ax = plt.subplots(figsize=(20,10))

sc = ax.scatter(days, tshifts, c=ccoefs, s=20, cmap="viridis", vmin=0, vmax=1.0)
ax.plot(days2, tshifts2, "k", alpha=0.5, linewidth=0.5)

if markdates:
    ndates = len(markdates)

    # plot skew linear fit
    for i in range(0, ndates+1):
        if i == 0:
            sd = daysno2[0]

            ed = markdates[i].date - reference_date.date
            ed = ed.days
        elif i <= ndates-1:
            sd = markdates[i-1].date - reference_date.date
            sd = sd.days

            ed = markdates[i].date - reference_date.date
            ed = ed.days
        else:
            sd = markdates[-1].date - reference_date.date
            sd = sd.days

            ed = daysno2[-1]

        i1 = np.argmin(np.abs(daysno2 - sd))
        i2 = np.argmin(np.abs(daysno2 - ed))

        if (i2 - i1) < 2:
            continue

        try:
            cjump = abs(tshifts2[i2+1] - tshifts2[i2])
            print(f"{ed} clock jump {cjump}")
        except:
            pass

        m1, m2 = estimate_clock_drift(daysno2[i1:i2],
                                      tshifts2[i1:i2],
                                      ccoefs2[i1:i2],
                                     )

        dps = m1 / 86400.0

        if i == 0:
            label = "linear drift using weighted lsqr"
        else:
            label = None

        ax.plot(days2[i1:i2], daysno2[i1:i2]*m1 + m2, "r",
                label=label)

        print(f"estimated linear drift per day/second {m1}/{dps}")

    for date in markdates:
        ax.axvline(x=date.date, c="g", label=f"{date.date} ({date.julday})")

if obs_sta:
    ax.axvline(x=deployment_time.date, c="k")
    ax.axvline(x=recover_time.date, c="k",
                  label="deployment and recovery day")

    sta1, sta2 = pair.split("_")
    if obs_sta in sta1:
        sign = 1.0
    elif obs_sta in sta2:
        sign = -1.0

    recdays = np.arange(deployment_time, recover_time+86400.0, dtype='datetime64[D]')
    recdaysno = np.arange(0, recdays.size)

    ax.plot(recdays, recdaysno * dpd * sign,
            "b", label="linear drift using measured skew",
              )

ax.set_title(f"{obs_sta} clock drift using {pair}", fontweight="bold", fontsize=15)
ax.set_xlabel("day",fontsize=15)
ax.set_ylabel("clock difference [s]",fontsize=15)
ax.tick_params(axis="both", labelsize=15)

ax.set_xlim(xlim1, xlim2)
ax.set_ylim(-4, 4)

ax.legend(loc="upper left",fontsize=14)
ax.grid()
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label("correlation coefficient", fontsize=15)

plt.show()
plt.close()

plot_corr(st_all, tr_ref, lags)#, wsig, wnoi, markdates)
