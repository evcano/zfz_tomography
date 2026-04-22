import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc


rc("text", usetex=True)

figsize = (6,3.5)
fs = 9

# NOTE: vindex indicates the period of noise parameter
fig, ax = plt.subplots(3)
fig.set_size_inches(figsize)

for i in range(1,21):
    fname = f"./final_chains_2million/chain_{i:02}/Results/samples_1.txt"

    step = np.loadtxt(fname,usecols=0,dtype=int)
#    accepted = np.loadtxt(fname,usecols=1,dtype=str)
#    vindex = np.loadtxt(fname,usecols=2,dtype=int)
    ncells = np.loadtxt(fname,usecols=3,dtype=int)
    misfit = np.loadtxt(fname,usecols=4,dtype=float)
    unweighted_misfit = np.loadtxt(fname,usecols=5,dtype=float)
    like = np.loadtxt(fname,usecols=6,dtype=float)
#    coor1 = np.loadtxt(fname,usecols=7,dtype=float)
#    coor2 = np.loadtxt(fname,usecols=8,dtype=float)
#    coor3 = np.loadtxt(fname,usecols=9, dtype=float)
#    val1 = np.loadtxt(fname,usecols=10, dtype=float)
#    val2 = np.loadtxt(fname,usecols=11, dtype=float)
#    val3 = np.loadtxt(fname,usecols=12, dtype=float)
#    noise0 = np.loadtxt(fname,usecols=13, dtype=float)
#    noise1 = np.loadtxt(fname,usecols=14, dtype=float)


    print("number of samples: ",len(step))

    ax[0].plot(like, label=str(i),lw=0.6,alpha=0.5)
    ax[1].plot(misfit, label=str(i),lw=0.6,alpha=0.5)
    ax[2].plot(ncells, label=str(i),lw=0.6,alpha=0.5)

for i in range(0,3):
    ax[i].axvline(500000, color="k", linestyle=":", linewidth=1.0)
    ax[i].set_xlim(-50000, 2e6)
    ax[i].tick_params(labelsize=10)

ax[0].text(170000, 4600, "Burn-in", fontsize=fs-1)
ax[1].text(170000, 2600, "Burn-in", fontsize=fs-1)
ax[2].text(170000, 220, "Burn-in", fontsize=fs-1)

ax[0].text(800000, 4600, "Exploration", fontsize=fs-1)
ax[1].text(800000, 2600, "Exploration", fontsize=fs-1)
ax[2].text(800000, 220, "Exploration", fontsize=fs-1)

ax[0].set_ylim(3900,5000)
ax[1].set_ylim(600, 4000)
ax[2].set_ylim(80,300)

ax[0].set_title(r"\textbf{(a)} Negative likelihood",fontsize=fs,loc="left")
ax[1].set_title(r"\textbf{(b)} Joint misfit",fontsize=fs,loc="left")
ax[2].set_title(r"\textbf{(c)} Number of cells",fontsize=fs,loc="left")

ax[0].set_ylabel("Negative likelihood",fontsize=fs)
ax[1].set_ylabel("Joint misfit",fontsize=fs)
ax[2].set_ylabel("Number of cells",fontsize=fs)

ax[0].set_xticklabels([])
ax[1].set_xticklabels([])

ax[2].set_xlabel("Step",fontsize=fs)

fig.tight_layout()
fig.savefig("./3dlike.png",dpi=300)
plt.close()
raise Exception

# keep only accepted samples
idx = np.where(accepted=="T")

step = step[idx]
accepted = accepted[idx]
vindex = vindex[idx]
ncells = ncells[idx]
misfit = misfit[idx]
unweighted_misfit = unweighted_misfit[idx]
like = like[idx]
coor1 = coor1[idx]
coor2 = coor2[idx]
coor3 = coor3[idx]
val1 = val1[idx]
val2 = val2[idx]
val3 = val3[idx]
noise0 = noise0[idx]
noise1 = noise1[idx]


# loop over rgrp, rpha, lgrp, lpha noise
for i in range(8,12):
    idx = np.where(step==i)

    n0= noise0[idx]
    n1= noise1[idx]
    vindtmp = vindex[idx]

    # loop over frequency
    for j in range(1,11):
        idx = np.argwhere(vindtmp==j)
        print(i, j, len(idx))

        plt.hist(n0[idx])
        plt.show()

        plt.hist(n1[idx])
        plt.show()
