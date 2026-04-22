import matplotlib.pyplot as plt
import numpy as np
import os


figpath = "./figures"

burnin = 500000
chains_no = np.arange(1,21)
steps_list = [8,9,10,11]

# -----------------------------------------------------------
freqs = np.arange(3,13)

n0_rgrp = []
n1_rgrp = []
pidx_rgrp = []

n0_rpha= []
n1_rpha = []
pidx_rpha = []

n0_lgrp = []
n1_lgrp = []
pidx_lgrp = []

n0_lpha= []
n1_lpha = []
pidx_lpha = []

for i in chains_no:
    fname = f"./final_chains_2million/chain_{i:02}/Results/samples_1.txt"

    step = np.loadtxt(fname,usecols=0,dtype=int)
    accepted = np.loadtxt(fname,usecols=1,dtype=str)
    vindex = np.loadtxt(fname,usecols=2,dtype=int)
    noise0 = np.loadtxt(fname,usecols=13, dtype=float)
    noise1 = np.loadtxt(fname,usecols=14, dtype=float)

    totalsamp = step.size

    # remove burn-in
    step = step[burnin:]
    accepted = accepted[burnin:]
    vindex = vindex[burnin:]
    noise0 = noise0[burnin:]
    noise1 = noise1[burnin:]

    # keep only accepted samples
    idx = np.where(accepted=="T")
    step = step[idx]
    accepted = accepted[idx]
    vindex = vindex[idx]
    noise0 = noise0[idx]
    noise1 = noise1[idx]

    accesamp = step.size

    print(f"chain {i} accepted samples: {accesamp}/{totalsamp}")

    # loop over rgrp, rpha, lgrp, lpha noise
    for i in steps_list:
        idx = np.where(step==i)

        n0= noise0[idx]
        n1= noise1[idx]
        pidx = vindex[idx]

        if i == 8:
            n0_rgrp.extend(n0)
            n1_rgrp.extend(n1)
            pidx_rgrp.extend(pidx)
        elif i == 9:
            n0_rpha.extend(n0)
            n1_rpha.extend(n1)
            pidx_rpha.extend(pidx)
        elif i == 10:
            n0_lgrp.extend(n0)
            n1_lgrp.extend(n1)
            pidx_lgrp.extend(pidx)
        elif i == 11:
            n0_lpha.extend(n0)
            n1_lpha.extend(n1)
            pidx_lpha.extend(pidx)

n0_rgrp = np.array(n0_rgrp)
n1_rgrp = np.array(n1_rgrp)
pidx_rgrp = np.array(pidx_rgrp)

n0_rpha = np.array(n0_rpha)
n1_rpha = np.array(n1_rpha)
pidx_rpha = np.array(pidx_rpha)

n0_lgrp = np.array(n0_lgrp)
n1_lgrp = np.array(n1_lgrp)
pidx_lgrp = np.array(pidx_lgrp)

n0_lpha = np.array(n0_lpha)
n1_lpha = np.array(n1_lpha)
pidx_lpha = np.array(pidx_lpha)

noise_pars = [
    ["Rayleigh group", n0_rgrp, n1_rgrp, pidx_rgrp],
    ["Rayleigh phase", n0_rpha, n1_rpha, pidx_rpha],
    ["Love group", n0_lgrp, n1_lgrp, pidx_lgrp],
    ["Love phase", n0_lpha, n1_lpha, pidx_lpha],
]

# figures settings
fs = 24

n0min = 0.00001
n0max = 0.2
n0step = 0.005

n1min = 0.0
n1max = 2.0
n1step = 0.05

bins_n0 = np.arange(n0min, n0max, n0step)
bins_n1 = np.arange(n1min, n1max, n1step)

## loop over frequency
for pars in noise_pars:
    label = pars[0]
    n0 = pars[1]
    n1 = pars[2]
    pidx = pars[3]

    fig, ax = plt.subplots(2,10)

    for i, fq in enumerate(freqs):
        idx = np.argwhere(pidx == i+1)

        ax[0,i].hist(n0[idx],
                   bins_n0,
                   histtype="bar",
                   color="lightskyblue",
                   edgecolor="blue",
                   linewidth=1.0,
                  )

        ax[0,i].axvline(np.mean(n0[idx]),c="k")
        ax[0,i].axvline(np.median(n0[idx]),c="b")

        ax[1,i].hist(n1[idx],
                   bins_n1,
                   histtype="bar",
                   color="lightskyblue",
                   edgecolor="blue",
                   linewidth=1.0,
                  )

        ax[1,i].axvline(np.mean(n1[idx]),c="k")
        ax[1,i].axvline(np.median(n1[idx]),c="b")

        # labels
        ax[0,i].set_xlabel(f"n0 {fq} s", fontsize=fs)
        ax[0,i].set_ylabel("count", fontsize=fs)
        ax[0,i].set_xlim(n0min, n0max)
        ax[0,i].tick_params(axis="x", labelsize=fs)
        ax[0,i].tick_params(axis="y", labelsize=fs)

        ax[1,i].set_xlabel(f"n1 {fq} s", fontsize=fs)
        ax[1,i].set_ylabel("count", fontsize=fs)
        ax[1,i].set_xlim(n1min, n1max)
        ax[1,i].tick_params(axis="x", labelsize=fs)
        ax[1,i].tick_params(axis="y", labelsize=fs)

    fig.set_size_inches(60,30)
    fig.suptitle(f"Data noise {label}", fontsize=fs)
    fname = f"{label}.png"
    fig.savefig(os.path.join(figpath,fname), bbox_inches="tight")

    print(label)
