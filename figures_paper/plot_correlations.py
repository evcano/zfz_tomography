import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from matplotlib.lines import Line2D
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from obspy import Stream, read

rc('text', usetex=True)

def plot_correlations(ax, data_path, data_format, cmp, pairs=None,
                      branch="both", maxtime=None, bpass=None,
                      global_normalization=False, yaxis=None, apparent_velocity=False,
                      shownames=False,
                      gain=1.0, lw=1.0, alpha=1.0):

    island = ["ZF.BREEM","ZF.QUMAN"]

    land = ["ZF.LAVA","ZF.KHUF","SA.WJHS","SA.EWJHS",
            "KL.DEEP","KL.EAST","KL.SOUTH","KL.WEST"]

    # read data
    st = Stream()

    if len(pairs) > 0:
        files = [f"{x}_{cmp}.{data_format}" for x in pairs]
        for f in files:
            stpath = os.path.join(data_path, f)
            if os.path.isfile(stpath):
                st += read(stpath, format=data_format)
    else:
        stpath = os.path.join(data_path, '*')
        st += read(stpath, format=data_format)

    ntr = len(st)

    # filter data
    if bpass:
        st.detrend("linear")
        st.detrend("demean")
        st.taper(0.1)
        st.filter('bandpass',
                  freqmin=bpass[0],
                  freqmax=bpass[1],
                  corners=2,
                  zerophase=True)

    # cut data
    if maxtime:
        maxlag = ((st[0].stats.npts - 1) / 2) * st[0].stats.delta
        for i in range(0, ntr):
            st[i] = st[i].slice(st[i].stats.starttime + maxlag - maxtime,
                                st[i].stats.starttime + maxlag + maxtime)

    maxtime = ((st[0].stats.npts - 1) / 2) * st[0].stats.delta
    lags = np.linspace(-maxtime, maxtime, st[0].stats.npts)

    # normalize data
    if global_normalization:
        st.normalize(global_max=True)
    else:
        st.normalize(global_max=False)

    # sort data according to interstation distance
    distances = []
    for tr in st:
        distances.append(tr.stats.sac.dist)

    idx = np.argsort(np.array(distances))
    distances = np.sort(distances)
    print(distances.min(), distances.max())

    data = np.zeros((ntr, st[0].stats.npts))
    ax2_labels = []

    for i, j in enumerate(idx):
        data[i, :] = st[j].data
        ax2_labels.append(f"{st[j].stats.sac.kevnm}_{st[j].stats.sac.kstnm}")

    # select branch
    if branch == "causal":
        data = data[:,np.where(lags>=0.0)]
        data = data.reshape(data.shape[0], data.shape[2])
        lags = lags[lags>=0.0]
    elif branch == "acausal":
        data = data[:,np.where(lags<=0.0)]
        data = data.reshape(data.shape[0], data.shape[2])
        data = np.fliplr(data)  # time reverse
        lags = lags[lags<=0.0]
        lags = -lags[::-1]  # time reverse

    # FIGURE
    # ------------------------------------------------------------------
    offset = 0

    # plot data
    for i in range(0, ntr):
        sta1, sta2 = ax2_labels[i].split("_")

        if sta1 in land or sta2 in land:
            c = "blueviolet"
        elif sta1 in island or sta2 in island:
            c = "crimson"
        else:
            c = "darkcyan"

        y = data[i,:] * gain

        if yaxis and yaxis == 'dis':
            y += distances[i]
        else:
            y += offset
            offset = np.max(y)

        ax.plot(lags, y, c=c, lw=lw, alpha=alpha)

        # station-pair name label
        if shownames:
            lab = ax2_labels[i]
            lab = lab.replace("ZF.","")
            ax.text(maxtime*0.2, np.mean(y)+0.2, lab,
                    fontsize=6.)

    if apparent_velocity:
        if yaxis and yaxis == 'dis':
            for s in apparent_velocity:
                ax.plot(s*np.array(distances), distances, 'b', lw=0.5,
                         alpha=alpha)

                #ax.text(s*distances[-30],
                #        distances[-2],
                #        f"$V_{{app}}$ = {1/s:.1f} km$\cdot$s$^{{-1}}$",
                #        fontsize=8
                #       )

                if branch == "both":
                    ax.plot(-s*np.array(distances), distances, 'b',
                             lw=0.4,alpha=alpha)

                    ax.text(-s*distances[-5],
                            distances[-5],
                             f"$V_{{app}}$ = {1/s:.1f} km$\cdot$s$^{{-1}}$",
                            fontsize=9,
                           )

    print('{} correlations plotted'.format(ntr))

    return ax

# PARAMETERS
# =========
data_path = "../noise_correlations/stacked_correlations_sorted_west_to_east/{cmp}"
figfile = f"./noise_correlations.png"

# data info
data_format = 'sac'
data_type = 'correlations'
bpass = [6., 12.]
maxlag = 180.0

global_normalization = False
yaxis = 'dis'

# figure settings
lw = 0.5
gain = 1.5
alpha = 0.6

fs = 9
fs2 = 7
figsize = (7.0, 6.5)
dpi = 400

xticks = [-maxlag,-100,-50,0,50,100,maxlag]
yticks = [15,50,75,100,125,150,175,200,230]

# -----------------------------------------------
ignore_sta = ["KL.DEEP","KL.EAST","KL.SOUTH","KL.WEST","SA.WJHS"]

# DONT EDIT BELOW THIS LINE
# =========================
fig, ax = plt.subplots(1,2)

for i, cmp in enumerate(["ZZ", "TT"]):
    if cmp == "ZZ":
        appvel = [1/2.1]
    else:
        appvel = [1/2.1]

    # list station pairs
    pairs = []
    files = os.listdir(data_path.format(cmp=cmp))

    for fname in files:
        s1, s2, _ = fname.split("_")

        if s1 in ignore_sta or s2 in ignore_sta:
            continue
        else:
            x = f"{s1}_{s2}"
            pairs.append(x)

    # plot figure
    ax[i] = plot_correlations(
        ax=ax[i],
        data_path=data_path.format(cmp=cmp),
        data_format=data_format,
        cmp=cmp,
        pairs=pairs,
        maxtime=maxlag,
        bpass=[1./bpass[1], 1./bpass[0]],
        global_normalization=global_normalization,
        yaxis=yaxis,
        apparent_velocity=appvel,
        gain=gain,
        lw=lw,
        alpha=alpha,
        shownames=False,
    )


    # figure settings
    labelsX = None
    labelsY = None
    yaxlabel = "Interstation distance [km]"

    if cmp == "TT":
        labelsY = []
        yaxlabel = None

    ax[i].set_xlim(-maxlag, maxlag)
    ax[i].set_ylim(10.0, 270.0)

    # ticks
    ax[i].set_xticks(xticks,labels=labelsX)
    ax[i].set_yticks(yticks,labels=labelsY)
    ax[i].tick_params(labelsize=fs)

    # labels
    ax[i].set_xlabel('Lag [s]', fontsize=fs)
    ax[i].set_ylabel(yaxlabel, fontsize=fs)

    # frame
    ax[i].spines["top"].set_visible(False)
    ax[i].spines["left"].set_bounds(10,yticks[-1])
    ax[i].spines["right"].set_bounds(10,yticks[-1])

    if cmp == "ZZ":
        let = r"\textbf{(a)}"
    else:
        let = r"\textbf{(b)}"

    ax[i].set_title(let + f" {cmp} noise correlations", fontsize=fs,
                    x=0.01,ha="left")

# INSET FIGURE ZZ ------------------------------------------------
ax2 = inset_axes(ax[0],
                 width="100%",
                 height="17%",
                 loc="center",
                 bbox_to_anchor=(0.0,0.42,1,1),
                 bbox_transform=ax[0].transAxes,
                )

pairs = ["ZF.OBS07_ZF.QUMAN",
         "ZF.OBS11_ZF.OBS04",
         "ZF.OBS06_ZF.KHUF",
        ]

ax2 = plot_correlations(
    ax=ax2,
    data_path=data_path.format(cmp="ZZ"),
    data_format=data_format,
    cmp="ZZ",
    pairs=pairs,
    maxtime=maxlag,
    bpass=[1./bpass[1], 1./bpass[0]],
    global_normalization=False,
    yaxis="no",
    gain=1.0,
    lw=lw,
    alpha=0.7,
    shownames=True,
)

# ticks
ax2.set_xticks(ticks=[],labels=[])
ax2.set_yticks(ticks=[],labels=[])
ax2.tick_params(labelsize=fs2, size=1.)

# lims
ax2.set_xlim(-maxlag,maxlag)
ax2.set_ylim(-1, 3.5)

# spines
for aname in ["top","bottom","left","right"]:
    ax2.spines[aname].set_visible(False)
ax2.patch.set_alpha(0.0)

# INSET FIGURE TT ------------------------------------------------
ax3 = inset_axes(ax[1],
                 width="100%",
                 height="17%",
                 loc="center",
                 bbox_to_anchor=(0.,0.42,1,1),
                 bbox_transform=ax[1].transAxes,
                )

pairs = ["ZF.OBS07_ZF.QUMAN",
         "ZF.OBS11_ZF.OBS04",
         "ZF.OBS10_ZF.KHUF",
        ]

ax3 = plot_correlations(
    ax=ax3,
    data_path=data_path.format(cmp="TT"),
    data_format=data_format,
    cmp="TT",
    pairs=pairs,
    maxtime=maxlag,
    bpass=[1./bpass[1], 1./bpass[0]],
    global_normalization=False,
    yaxis="no",
    gain=1.0,
    lw=lw,
    alpha=0.7,
    shownames=True,
)

# ticks
ax3.set_xticks(ticks=[],labels=[])
ax3.set_yticks(ticks=[],labels=[])
ax3.tick_params(labelsize=fs2, size=1.)

# lims
ax3.set_xlim(-maxlag,maxlag)
ax3.set_ylim(-1, 3.5)

# spines
for aname in ["top","bottom","left","right"]:
    ax3.spines[aname].set_visible(False)
ax3.patch.set_alpha(0.0)

# LEGEND WITH LABELS -------------------------------------------
lhandles = []
lhandles.append(Line2D([0],[0],color="blueviolet",linewidth=1.0,
                label='OBS-Land'))
lhandles.append(Line2D([0],[0],color="crimson",linewidth=1.0,
                label='OBS-Island'))
lhandles.append(Line2D([0],[0],color="darkcyan",linewidth=1.0,
                label='OBS-OBS'))

leg = ax3.legend(handles=lhandles,
                 loc="center",
                 bbox_to_anchor=(-0.6,0.,1,1),
                 fontsize=fs2,
                 framealpha=1.,
                 edgecolor="k",
                 fancybox=False,
                )

leg.get_frame().set_linewidth(0.5)

for line in leg.get_lines():
    line.set_linewidth(1.0)

# SAVE FIGURE ---------------------------------------------------
fig.subplots_adjust(wspace=2.0)
fig.set_size_inches(figsize)
fig.tight_layout()
fig.savefig(figfile, dpi=dpi)
