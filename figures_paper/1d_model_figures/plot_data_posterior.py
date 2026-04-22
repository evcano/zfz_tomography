import os
import numpy as np
import os.path as op
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from BayHunter import ModelMatrix
from BayHunter import Model
from BayHunter import PlotFromStorage
from BayHunter import Targets
from matplotlib import rc


rc('text',usetex=True)

def plot_posterior_noise(ax, noise):
    bins = 20
    formatter = "%d"

    ax.hist(noise,bins=bins,color="blue",edgecolor="white")

    median = np.median(noise)
    ax.text(0.5,0.8, f"median: {median:.4f}"+" km$\cdot$s$^{-1}$",
            fontsize=9, color="k",
            transform=ax.transAxes)
    ax.axvline(median,color="k",lw=1,ls=":")


def plot_last_datafit(ax, pobj, chains, mpath, vpvs=1.73):
    targets = Targets.JointTarget(targets=pobj.targets)

    for n, target in enumerate(targets.targets):
        ax[n,0].plot(target.obsdata.x,target.obsdata.y,
                     label="Observed", marker="x",ms=1,
                     color="blue",linewidth=0.7,zorder=10000)

    base = cm.get_cmap(name="rainbow")
    colist = base(np.linspace(0,1,len(chains)))

    for i, c in enumerate(chains):
        fname = op.join(mpath,f"c{c:03}_p2models.npy")
        models = np.load(fname)
        lastmodel = models[-1]

        vp, vs, h = Model.get_vp_vs_h(lastmodel,vpvs)
        rho = vp*0.32 + 0.77

        for n, target in enumerate(targets.targets):
            xmod, ymod = target.moddata.plugin.run_model(
                h=h,vp=vp,vs=vs,rho=rho)

            yobs = target.obsdata.y

            ax[n,0].plot(xmod,ymod,c=colist[i],linewidth=0.5,alpha=0.3)


resultsdir = "../1d_inversion/results/data"

configfile = op.join(resultsdir,"obs_config.pkl")

thining = 1  # keep at 1
nchains = 100
outliers = [38,48,60]

vpvs = 1.73
depint = np.arange(0,21,0.5)

# figure settings
figname = "1d_data_posterior.png"
dpi = 400

figsize = (7,5)
fs = 9
fs2 = 8

xticks = np.arange(2,11)
xticks2 = np.arange(0.0,0.25,0.04)
#----------------------------------------------
chains = [x for x in np.arange(0,nchains) if x not in outliers]

# setup figure
fig, ax = plt.subplots(4,2,figsize=figsize)

# load vs posterior
fname = op.join(resultsdir, f"c_models.npy")
models = np.load(fname)
models = models[::thining]

# load noise posterior
fname = op.join(resultsdir, f"c_noise.npy")
noise = np.load(fname)
noise = noise.T

# setup plotting object
pobj = PlotFromStorage(configfile)

# ORDER OF DATA IS: LOVE-GROUP,LOVE-PHASE,RAYLEIGH-GROUP,RAYLEIGH-PHASE
# plot data fit
plot_last_datafit(ax, pobj, chains, resultsdir, vpvs=1.73)


# plot noise posterior
for j, i in enumerate(range(1,8,2)):
    data = noise[i,:]
    plot_posterior_noise(ax[j,1], data)

    ax[j,1].set_yticks([])

titles = [
    "Love waves group-velocity",
    "Love waves phase-velocity",
    "Rayleigh waves group-velocity",
    "Rayleigh waves phase-velocity",
]

tlet1 = [r"\textbf{(a)}",r"\textbf{(c)}",r"\textbf{(e)}",r"\textbf{(g)}"]
tlet2 = [r"\textbf{(b)}",r"\textbf{(d)}",r"\textbf{(f)}",r"\textbf{(h)}"]


# figures style
for i in range(0,4):
    ymin = np.around(ax[i,0].get_ylim()[0],2) - 0.3
    ymax = np.around(ax[i,0].get_ylim()[1],2) + 0.3
    yticks = np.linspace(ymin,ymax,4)

    ax[i,0].set_xlim(xticks[0],xticks[-1])
    ax[i,0].set_ylim(yticks[0], yticks[-1])

    ax[i,0].set_xticks(xticks,labels=[])
    ax[i,0].set_yticks(np.around(yticks,2))
    ax[i,0].tick_params(labelsize=fs)


    ax[i,0].set_xlabel("")
    ax[i,0].set_ylabel("Velocity [km$\cdot$s$^{-1}$]",
                       fontsize=fs)

    ax[i,0].set_title(tlet1[i]+" Data fit: "+titles[i],fontsize=fs,
                      x=0.01,ha="left")

    if i == 3:
        ax[i,0].set_xticks(xticks,labels=xticks)
        ax[i,0].set_xlabel("Period [s]",fontsize=fs)

    ax[i,0].legend(loc="upper left",
                   fancybox=False,
                   edgecolor="k",
                   framealpha=1.0,
                   fontsize=fs,
                  )

# ---
for i in range(0,4):
    ax[i,1].set_xlim(xticks2[0],xticks2[-1])

    ax[i,1].set_xticks(xticks2,labels=[])

    ax[i,1].tick_params(labelsize=fs)

    ax[i,1].set_title(tlet2[i]+" Data noise: "+titles[i],fontsize=fs,
                      x=0.01,ha="left")

    if i == 3:
        ax[i,1].set_xticks(xticks2,labels=xticks2)
        ax[i,1].set_xlabel("$\sigma$ [km$\cdot$s$^{-1}$]",fontsize=fs)

fig.tight_layout()

if figname:
    fig.savefig(figname)
else:
    plt.show()

plt.show()
