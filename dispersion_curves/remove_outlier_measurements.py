import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from sanpy.base.project_functions import load_project
from sanpy.util.raymap import plot_ray_coverage


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
        label = "coast"

    pair = f"{s1}_{s2}"

    return label, pair


project_path = "../../../noise_correlations_tomo/zf_correlation.pkl"

curvepath = "./disp_cur_TT/manually_checked/coast"
outpath = "./disp_cur_TT/final_curves/coast"

thr = 3.0  # no of std
pmin = 1.
pmax = 15.
dtp = 0.5

# ---------------------------------------------------
P = load_project(project_path)
per_axis = np.arange(pmin, pmax+dtp, dtp)

# loop over group and phase curves
for prefix in ["GDisp.", "CDisp.T."]:
    files = glob(os.path.join(curvepath,f"{prefix}*"))
    files.sort()

    pairs = []
    labels = []
    vel_arr = np.zeros((len(files), per_axis.size))
    vel_arr[:] = np.nan

    # read curves and store in a single array
    for i, fname in enumerate(files):
#        label, pair = label_correlation(fname)
#        pairs.append(pair)
#        labels.append(label)
#
        X = np.loadtxt(fname)
        per = X[:,0]
        vel = X[:,1]
        snr = X[:,2]

        for j in range(per.size):
            idx = np.where(per_axis == per[j])
            vel_arr[i,idx] = vel[j]

    # compute mean and std
    vel_mean = np.nanmean(vel_arr, axis=0)
    vel_std = np.nanstd(vel_arr, axis=0)

    # remove outliers
    for i in range(len(files)):
        for j in range(per_axis.size):
            if vel_arr[i,j] > (vel_mean[j] + thr*vel_std[j]):
                vel_arr[i,j] = np.nan
            elif vel_arr[i,j] < (vel_mean[j] - thr*vel_std[j]):
                vel_arr[i,j] = np.nan

    # save curves without outliers
    for i in range(len(files)):
        outfile = os.path.basename(files[i])
        outfile = os.path.join(outpath,outfile)

        out_per = per_axis[~np.isnan(vel_arr[i,:])]
        out_vel = vel_arr[i, ~np.isnan(vel_arr[i,:])]
        out_arr = np.array([out_per, out_vel]).T

        np.savetxt(outfile, out_arr, fmt="%.6f")


    # count measurements at each period
    counter = np.sum(~np.isnan(vel_arr), axis=0)

    # plot curves
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    for i in range(len(files)):
        ax1.plot(per_axis, vel_arr[i,:], c="k", alpha=0.2)
        ax1.plot(per_axis, vel_arr[i,:], '.', c="k", alpha=0.2)

    ax1.plot(per_axis, vel_mean, c="b")
    ax1.plot(per_axis, vel_mean + thr*vel_std, c="r")
    ax1.plot(per_axis, vel_mean - thr*vel_std, c="r")
    ax2.bar(per_axis, counter, width=dtp, align="center", color="b", edgecolor="k")

    ax1.set_xlim(1, 20)
    ax1.set_ylim(0.5, 5)
    ax1.set_title(f"{prefix} dispersion curves")

    ax2.set_xlim(1, 20)
    ax2.set_ylim(0, 100)
    ax2.set_title(f"No. of rays per period ({prefix})")

    plt.show()
    plt.close()
