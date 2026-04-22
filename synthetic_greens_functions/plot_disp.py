import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc

rc("text", usetex=True)

refpath = "/home/eduardo/projects/red_sea_obs_gji/dispersion_curves/EGFAnalysisTimeFreq/ZAFRAN/reference_curves"
files = ["./SREGN.ASC", "./SLEGN.ASC"]


fig, axes = plt.subplots(1,2,figsize=(9.5, 5),constrained_layout=True)

for ii, fname in enumerate(files):
    X = np.loadtxt(fname, skiprows=(1))
    
    cols = ["r","g","b"]
    ax = axes[ii]


    if fname == "./SLEGN.ASC":
        title = "b) Love wave dispersion curves"
        refc1 = np.loadtxt(os.path.join(refpath, "tt-group-obs-twindow.txt"))
        refc2 = np.loadtxt(os.path.join(refpath, "tt-phase-obs-refcurve.txt"))
    elif fname == "./SREGN.ASC":
        title = "a) Rayleigh wave dispersion curves"
        refc1 = np.loadtxt(os.path.join(refpath, "zz-group-obs-twindow.txt"))
        refc2 = np.loadtxt(os.path.join(refpath, "zz-phase-obs-refcurve.txt"))

    for mn in [0, 1, 2]:
        idx = np.where(X[:,0]==mn)
    
        t = X[idx,2]
        c = X[idx,4]
        u = X[idx,5]
    
        t = t[0]
        c = c[0]
        u = u[0]
    
        ax.plot(t,u,"-", c=cols[mn], label=f"Mode {mn} group velocity", linewidth=0.7)
        ax.plot(t,c,"--",c=cols[mn], label=f"Mode {mn} phase velocity", linewidth=0.7)

    ax.plot(refc1[:,0],refc1[:,2],c="k", linewidth=0.7, label="Group velocity bound") 
    ax.plot(refc1[:,0],refc1[:,3],c="k", linewidth=0.7)

    ax.plot(refc2[:,0],refc2[:,1],"--", c="k", linewidth=0.7, label="Phase velocity reference") 
    ax.plot(refc2[:,0],refc2[:,1],"--", c="k", linewidth=0.7)

    ax.set_xlabel("Period [s]")
 
    ax.set_ylim(0, 5) 
    ax.set_xlim(0, 17)

    if ii == 1:
        ax.legend(loc="lower right",bbox_to_anchor=(0.2,0.05,1,1),fancybox=False)
        ax.set_yticklabels([])
    else:
        ax.set_ylabel("Velocity [km/s]")

    ax.set_title(title)

plt.savefig(fname+".png")
plt.show()
plt.close()
