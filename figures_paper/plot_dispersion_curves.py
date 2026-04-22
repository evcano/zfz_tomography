import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from matplotlib import rc


rc("text", usetex=True)

def detect_jumps(per, vel, dtp):
    # based on velocity
    pdiff = np.diff(vel)
    idx_jump = np.argwhere(pdiff > 0.3)
    nseg = idx_jump.size + 1

    if nseg > 1:
        seg_idx = []
        i2 = -1
        for i in range(0, nseg):
            i1 = i2 + 1
            if i == nseg - 1:
                i2 = per.size - 1
            else:
                i2 = idx_jump[i][0]
            seg_idx.append(np.arange(i1, i2 + 1,dtype=int))
    else:
        seg_idx = [np.arange(0, per.size, dtype=int)]

    return seg_idx


def reverse_pair_name(pair):
    s1, s2 = pair.split("_")
    pair = f"{s2}_{s1}"
    return pair


def label_correlation(fname):
    obs = [f"ZF.OBS{x:02}" for x in range(1,13)]
    obs.extend(["ZF.NORTH", "ZF.SOUTH"])

    island = ["ZF.BREEM", "ZF.QUMAN"]

    coast = ["ZF.LAVA","ZF.KHUF","SA.EWJHS"]


    fname = os.path.basename(fname)
    fname = fname.replace(f"{prefix}","")
    fname = fname.replace(".dat","")
    s1, s2, _ = fname.split("_")

    if s1 not in obs  and s2 not in obs:
        raise Exception("no obs involved in correlation")
    elif s1 in obs  and s2 in obs:
        label = "obs"
    elif s1 in island or s2 in island:
        label = "island"
    elif s1 in coast or s2 in coast:
        label = "land"

    pair = f"{s1}_{s2}"

    return label, pair



curvepath0 = "../dispersion_curves/disp_cur_{cmp}/final_curves/all"

# CURVES SETTINGS
pmin = 3.
pmax = 12.
dtp = 0.5

# FIGURE SETTINGS
figname = "dispersion_curves.pdf"
figsize = (7.0, 4.)
dpi = 400
fs = 9
fs2 = 7

xticks = np.arange(2,15,2)
yticks = np.arange(0,5,dtype=float)
yticks2 = [0,50,100]

# ---------------------------------------------------
colors = {"obs":"darkcyan","island":"crimson","land":"blueviolet"}
linestyles = {"GDisp.":"solid", "CDisp.T.":"dashed"}#(0, (5,3))}
vtype = {"GDisp.": "Group velocity", "CDisp.T.": "Phase velocity"}

per_axis = np.arange(pmin, pmax+dtp, dtp)

# initialize figures
fig, ax = plt.subplots(2,2,gridspec_kw={"height_ratios":[2,1]})

for k, cmp in enumerate(["ZZ", "TT"]):
    curvepath = curvepath0.format(cmp=cmp)

    for prefix in ["GDisp.","CDisp.T."]:
        files = glob(os.path.join(curvepath,f"{prefix}*"))
        files.sort()

        vel_arr = np.zeros((len(files), per_axis.size))
        vel_arr[:] = np.nan

        # plot curves
        for i, fname in enumerate(files):
            label, pair = label_correlation(fname)

            label2 = label.capitalize()
            label2 = label2.replace("Obs","OBS")
            vtype2 = vtype[prefix].replace('velocity','')
            llabel = f"OBS-{label2} {vtype2}"

            X = np.loadtxt(fname)
            per = X[:,0]
            vel = X[:,1]

            ii = (per >= pmin) & (per <= pmax)
            per = per[ii]
            vel = vel[ii]

            for j in range(per.size):
                idx = np.where(per_axis == per[j])
                vel_arr[i,idx] = vel[j]

            # segment dispersion curves to avoid plotting gaps
            seg_idx = detect_jumps(per, vel, dtp)
            for sidx in seg_idx:
                ax[0,k].plot(per[sidx], vel[sidx],
                             linestyle=linestyles[prefix], c=colors[label],
                             linewidth=0.3, alpha=0.5,
                             label=llabel,
                            )

        # compute mean
        vel_mean = np.nanmean(vel_arr, axis=0)
        counter = np.sum(~np.isnan(vel_arr), axis=0)

        ## plot mean curve
        #ax[0,k].plot(per_axis, vel_mean, linewidth=0.5, c="k",
        #             linestyle=linestyles[prefix])

        # histogram of rays
        if vtype[prefix] == "Group velocity":
            marker = "."
        else:
            marker = "x"

        ax[1,k].plot(per_axis[::2], counter[::2], marker, lw=1.0, markersize=3.0,
                     color="k",label=vtype[prefix])

        # top subplot settings -------------------------------------------------
        ax[0,k].set_xticks(xticks)
        ax[0,k].set_yticks(yticks)
        ax[0,k].tick_params(labelsize=fs)

        ax[0,k].set_xlabel("Period [s]", fontsize=fs)
        ax[0,k].set_ylabel("Velocity [km/s]", fontsize=fs)

        ax[0,k].set_xlim(xticks[0],xticks[-1])
        ax[0,k].set_ylim(yticks[0], yticks[-1])

        ax[0,k].spines["left"].set_bounds(yticks[0],yticks[-1])

        if k == 1:
            ax[0,k].set_yticks(yticks,labels=[])
            ax[0,k].set_ylabel(None)

            # legend
            handles, llabels = ax[0,k].get_legend_handles_labels()
            by_label = dict(zip(llabels,handles))
            handles = by_label.values()
            llabels = by_label.keys()
            llabels, handles = zip(*sorted(zip(llabels, handles), key=lambda t: t[0]))

            leg = ax[0,k].legend(handles,
                                 llabels,
                                 loc="upper right",
                                 bbox_to_anchor=(0.2,0.2,1,1),
                                 fontsize=fs2,
                                 framealpha=0.8,
                                 edgecolor="k",
                                 fancybox=False,
                                )

            for line in leg.get_lines():
                line.set_linewidth(1.0)

        # bottom subplot settings -------------------------------------------------
        ax[1,k].set_xticks(xticks)
        ax[1,k].set_yticks(yticks2)
        ax[1,k].tick_params(labelsize=fs)

        ax[1,k].set_xlabel("Period [s]", fontsize=fs)
        ax[1,k].set_ylabel("Count", fontsize=fs)

        ax[1,k].set_xlim(xticks[0],xticks[-1])
        ax[1,k].set_ylim(yticks2[0], yticks2[-1])

        ax[1,k].spines["left"].set_bounds(yticks2[0],yticks2[-1])

        if k == 1:
            ax[1,k].set_yticks(yticks2,labels=[])
            ax[1,k].set_ylabel(None)

            # legend
            leg = ax[1,k].legend(loc="upper right",
                                 bbox_to_anchor=(0.15,0.1,1,1),
                                 fontsize=fs2,
                                 framealpha=1,
                                 edgecolor="k",
                                 fancybox=False,
                                 )

ax[0,0].set_title(r"\textbf{(a)} Rayleigh wave dispersion curves",fontsize=fs,
                  loc="left")
ax[0,1].set_title(r"\textbf{(b)} Love wave dispersion curves",fontsize=fs,
                  loc="left")
ax[1,0].set_title(r"\textbf{(c)} Number of Rayleigh wave measurements",fontsize=fs,
                  loc="left")
ax[1,1].set_title(r"\textbf{(d)} Number of Love wave measurements",fontsize=fs,
                  loc="left")

# savefigures
fig.set_size_inches(figsize)
fig.tight_layout()
fig.savefig(figname, dpi=dpi)
