import os
import numpy as np
import os.path as op
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as pltticker
from BayHunter import ModelMatrix
from BayHunter import Model
from BayHunter.Plotting import vs_round
from matplotlib import rc


rc('text',usetex=True)

def plot_models(ax,mpath,mtype,chains,wl=0.0):
    base = cm.get_cmap(name="rainbow")
    colist = base(np.linspace(0,1,len(chains)))

    for i,c in enumerate(chains):
        fname = os.path.join(mpath,f"model_{c:03}_{mtype}.dat")
        X = np.loadtxt(fname)
        X[:,0] += wl
        ax.plot(X[:,1],X[:,0],c=colist[i],linewidth=0.6,alpha=0.6)
    ax.invert_yaxis()

    if wl > 0.0:
        ax.fill_between(x=X[:,0],y1=0,y2=wl,color="lightblue")


def plot_bestmodels_hist(ax1, ax2, models, dep_int=None,wl=0.0):
    """
    2D histogram with 30 vs cells and 50 depth cells.
    As plot depth is limited to 100 km, each depth cell is a 2 km.

    pinterf is the number of interfaces to be plot (derived from gradient)
    """
    if dep_int is None:
        dep_int = np.linspace(0, 100, 201)  # interppolate depth to 0.5 km.
        depbins = np.linspace(0, 100, 101)  # 1 km bins
    else:
        maxdepth = int(np.ceil(dep_int.max()))
        interp = dep_int[1] - dep_int[0]
        dep_int = np.arange(dep_int[0], dep_int[-1] + interp / 2., interp / 2.)
        depbins = np.arange(0, maxdepth + 2*interp, interp)  # interp km bins

    # get interfaces, #first
    models2 = ModelMatrix._replace_zvnoi_h(models)
    models2 = np.array([model[~np.isnan(model)] for model in models2])

    yinterf = np.array([np.cumsum(model[int(model.size/2):-1])
                        for model in models2])
    yinterf = np.concatenate(yinterf)
    yinterf += wl

    vss_int, deps_int = ModelMatrix.get_interpmodels(models, dep_int)
    deps_int += wl

    singlemodels = ModelMatrix.get_singlemodels(models, dep_int=depbins)

    vss_flatten = vss_int.flatten()
    vsinterval = 0.025  # km/s, 0.025 is assumption for vs_round

    # vs bins
    vs_histmin = vs_round(vss_flatten.min())-2*vsinterval - 1
    vs_histmax = vs_round(vss_flatten.max())+3*vsinterval
    vsbins = np.arange(vs_histmin, vs_histmax, vsinterval) # some buffer

    # initiate plot
    data2d, xedges, yedges = np.histogram2d(vss_flatten,
                                            deps_int.flatten(),
                                            bins=(vsbins, depbins),
                                            density=False,
                                           )

    ax1.imshow(data2d.T,
               extent=(xedges[0], xedges[-1],yedges[0], yedges[-1]),
               origin='lower',
               vmax=len(models),
               aspect='auto',
              )

    # plot mean
    vs, dep = singlemodels["mean"]
    dep += wl

    ax1.plot(vs, dep, color="green",
             lw=1., alpha=0.6, label="Average")

    # plot standard deviations
    stdmin = singlemodels["stdminmax"][0][0,:]
    stdmax = singlemodels["stdminmax"][0][1,:]
    dep = singlemodels["stdminmax"][1]

    ax1.plot(stdmin,dep,color="red",
             linewidth=0.8,alpha=0.6,
             label="Standard\ndeviation")

    ax1.plot(stdmax,dep,color="red",
             linewidth=1.,alpha=0.6)

    # histogram for interfaces
    data = ax2.hist(yinterf,
                    bins=depbins,
                    orientation='horizontal',
                    color='lightgray',
                    alpha=0.7,
                    edgecolor='k',
                   )
    bins, lay_bin, _ = np.array(data).T
    center_lay = (lay_bin[:-1] + lay_bin[1:]) / 2.

    ax1.invert_yaxis()
    ax2.invert_yaxis()

    plt.show()

def plot_posterior_layers(ax, models):
    models = np.array([model[~np.isnan(model)] for model in models])
    layers = np.array([(model.size/2 - 1) for model in models])
    bins = np.arange(np.min(layers), np.max(layers)+2)-0.5
    formatter = "%d"
    ax.hist(layers,bins=bins,color="blue",edgecolor="white")
    median = np.median(layers)
    ax.text(0.5,0.8,f"median: {median}",
            fontsize=9, color="k",
            transform=ax.transAxes)
    ax.axvline(median,color="k",lw=1,ls=":")

def plot_iter(ax,resultsdir,chains,iburn,imax,ptype):
    base = cm.get_cmap(name="rainbow")
    colist = base(np.linspace(0,1,len(chains)))

    for i,c in enumerate(chains):
        for p in ["p1","p2"]:
            fname = op.join(resultsdir, f"c{c:03}_{p}{ptype}.npy")
            data = np.load(fname)
            if ptype == "misfits":
                data = data.T[-1]
            if p == "p1":
                x = np.linspace(0,iburn,data.size)
                a = 0.1
            else:
                x = np.linspace(iburn,imax,data.size)
                a = 0.5
            ax.plot(x,data,color=colist[i],alpha=a)

    ax.axvline(iburn,color="k",ls=":",alpha=0.7)


resultsdir = "../1d_inversion/results/data"
modelsdir = "./models_per_chain"
wl = 1.2

thining = 1  # keep thining at 1, since i am reading thined chains
nchains = 100
outliers = [38,48,60]  # 38 and 48 are outliers, 60 had an errror
iburn = 2000000
imax = 3000000

vpvs = 1.73
depint = np.arange(0,24,0.5)

# figure settings
figname = "1d_model_posterior.png"
dpi = 300

figsize = (7,5)
fs = 9
fs2 = 8
xticks = [2, 3, 4, 5, 6, 7]
yticks = np.arange(0,23,2)

#----------------------------------------------
# setup figure
fig, ax = plt.subplot_mosaic(
    [['a','b','c','d1'],
     ['a','b','c','d2'],
     ['a','b','c','d3']],
    figsize=figsize,
    layout='constrained',
    gridspec_kw={
        "width_ratios":[1,1,0.2,1],
    },
)

chains = [x for x in np.arange(0,nchains) if x not in outliers]

# load vs posterior
fname = op.join(resultsdir, f"c_models.npy")
models = np.load(fname)
models = models[::thining]

hur = 11. + wl
# -------------------------------------------------
# best models
plot_models(ax['a'], modelsdir, "final", chains,wl=wl)

ax['a'].grid(linewidth=0.5,alpha=0.5)

ax['a'].fill_between([0,100],[hur,hur],[26,26],
                     color="gray",
                     edgecolor=None,
                     alpha=0.4,)

ax['a'].set_xlim(xticks[0]-0.5,xticks[-1])
ax['a'].set_ylim(yticks[-1],yticks[0])

ax['a'].set_xticks(xticks)
ax['a'].set_yticks(yticks)

ax['a'].xaxis.set_minor_locator(pltticker.MultipleLocator(base=0.5))
ax['a'].yaxis.set_minor_locator(pltticker.MultipleLocator(base=0.5))

ax['a'].set_xlabel("$V_s$ [km $\cdot$ s$^{-1}$]")
ax['a'].set_ylabel('Depth [km]',
                   fontsize=fs)

ax['a'].set_title(r'\textbf{(a)} $V_s$ last models',fontsize=fs,
                  x=0.01,ha="left")

# vs and interface posterior
plot_bestmodels_hist(ax['b'],ax['c'],models, depint, wl=wl)

ax['b'].grid(linewidth=0.5,alpha=0.3)

ax['b'].fill_between([0,100],[hur,hur],[26,26],
                     color="gray",
                     edgecolor=None,
                     alpha=0.4,)

ax['b'].legend(loc=3,
               fancybox=False,
               framealpha=1.0,
               edgecolor="k",
               fontsize=fs,
              )

ax['b'].set_xlabel("$V_s$ [km $\cdot$ s$^{-1}$]")

ax['b'].set_xlim(xticks[0]-0.5,xticks[-1])
ax['b'].set_ylim(yticks[-1],yticks[0])

ax['b'].set_xticks(xticks)
ax['b'].set_yticks(yticks, labels=[])

ax['b'].xaxis.set_minor_locator(pltticker.MultipleLocator(base=0.5))
ax['b'].yaxis.set_minor_locator(pltticker.MultipleLocator(base=0.5))

mtext = f"{len(models)} models\nfrom {len(chains)} chains"
ax['b'].text(3.5,3.,mtext,fontsize=fs,color="w")

mtext = "High uncertainty\nregion"
ax['b'].text(2, 15.,mtext,fontsize=fs,color="w")

ax['b'].set_title(r'\textbf{(b)} $V_s$ posterior',fontsize=fs,
                  x=0.01,ha="left")

ax['c'].set_ylim(yticks[-1],yticks[0])
ax['c'].set_xticks([])
ax['c'].set_yticks(yticks, labels=[])
ax['c'].yaxis.set_minor_locator(pltticker.MultipleLocator(base=0.5))

ax['c'].set_title(r'\textbf{(c)} Interface posterior',fontsize=fs)

if wl > 0.0:
    ax['b'].fill_between(x=ax['b'].get_xlim(),y1=0,y2=wl,color="lightblue")
    ax['c'].fill_between(x=ax['c'].get_xlim(),y1=0,y2=wl,color="lightblue")

# posterior layers
plot_posterior_layers(ax['d1'],models)
xticks2 = np.arange(2,14,2)

ax['d1'].set_xticks(xticks2)
ax['d1'].set_yticks([])

ax['d1'].set_xlim(xticks2[0],xticks2[-1])

ax['d1'].set_xlabel("No. of layers",fontsize=fs)
ax['d1'].set_title(r'\textbf{(d)} No. of layers posterior',fontsize=fs,
                   x=0.01,ha="left")

# likelihood
plot_iter(ax['d2'],resultsdir,chains,iburn,imax,"likes")

ax['d2'].text(iburn/2, 50,"Burn-in",
              color="k", ha="center",va="top",
              fontsize=fs2)

ax['d2'].text(1.01*(iburn+imax)/2, 50,"Exploration",
              color="k", ha="center",va="top",
              fontsize=fs2)

ax['d2'].set_xlabel("Step",fontsize=fs)
ax['d2'].set_ylabel("Likelihood",fontsize=fs)

ax['d2'].set_ylim(20, 120)
ax['d2'].set_title(r'\textbf{(e)} Likelihood',fontsize=fs,
                   x=0.01,ha="left")

# misfits
plot_iter(ax['d3'],resultsdir,chains,iburn,imax,"misfits")

ax['d3'].text(iburn/2, 0.15, "Burn-in",
              color="k", ha="center",va="top",
              fontsize=fs2)

ax['d3'].text(1.01*(iburn+imax)/2, 0.15,"Exploration",
              color="k", ha="center",va="top",
              fontsize=fs2)

ax['d3'].set_xlabel("Step",fontsize=fs)
ax['d3'].set_ylabel("Misfit",fontsize=fs)

ax['d3'].set_ylim(0.1,0.35)
ax['d3'].set_title(r'\textbf{(f)} Joint misfit',fontsize=fs,
                   x=0.01,ha="left")

for l in ax.keys():
    ax[l].tick_params(labelsize=fs)

if figname:
    fig.savefig(figname,dpi=dpi)
else:
    plt.show()
