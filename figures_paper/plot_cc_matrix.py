import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.colors import ListedColormap
from obspy import read


datapath = "../noise_correlations_tomo/daily_corr"

stations = [
 "ZF.OBS01",
 "ZF.OBS02",
 "ZF.OBS03",
 "ZF.OBS04",
 "ZF.OBS05",
 "ZF.OBS06",
 "ZF.OBS07",
 "ZF.OBS08",
 "ZF.OBS09",
 "ZF.OBS10",
 "ZF.OBS11",
 "ZF.OBS12",
 "ZF.NORTH",
 "ZF.SOUTH",
 "ZF.QUMAN",
 "ZF.BREEM",
 "ZF.LAVA",
 "ZF.KHUF",
 "SA.EWJHS",
]

nsta = len(stations)

fig, axes = plt.subplots(1,2, figsize=(9.,5.5), constrained_layout=True)

for cmp in ["ZZ","TT"]:
    pairs = []
    counts = []

    for sta1 in stations:
        for sta2 in stations:
            pair1 = f"{sta1}_{sta2}"
            pair2 = f"{sta2}_{sta1}"

            dir1 = os.path.join(datapath, cmp, pair1)
            dir2 = os.path.join(datapath, cmp, pair2)
    
            if os.path.exists(dir1):
                count = len([f for f in os.listdir(dir1) if os.path.isfile(os.path.join(dir1, f))])
            elif os.path.exists(dir2):
                count = len([f for f in os.listdir(dir2) if os.path.isfile(os.path.join(dir2, f))])
            else:
                count = np.nan

            if sta1 == sta2:
                count = 0

            if count == 0:
                count = np.nan

            pairs.extend([pair1,pair2])
            counts.extend([count, count])

    X = np.zeros((nsta,nsta))
    for i, sta1 in enumerate(stations):
        for j, sta2 in enumerate(stations):
            pair = f"{sta1}_{sta2}"
            idx = pairs.index(pair)
            X[i,j] = counts[idx]

    counts = np.array(counts)
    print(np.nanmin(counts[:]))
    print(np.nanmax(counts[:]))

    vmin = 50
    vmax = 400
    ncolors = int((vmax-vmin) / 10)
    cmap = ListedColormap(plt.cm.Spectral(np.linspace(0, 1, ncolors)))

    ticks = np.arange(0.5, nsta+0.5, 1)

    if cmp == "ZZ":
        ax = axes[0]
        title = f"a) Number of daily {cmp} noise correlations"
    else:
        ax = axes[1]
        title = f"b) Number of daily TT (RR) noise correlations"

    im = ax.pcolormesh(X, edgecolors="k", linewidth=1, cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_xticks(ticks, labels=stations)
    ax.set_yticks(ticks, labels=stations)
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_title(title)

    if cmp == "RR" or cmp == "ZZ":
        pass
    else:
        ax.set_yticklabels([])
        cbar = plt.colorbar(mappable=im, label="Count", shrink=0.7)

    plt.setp(ax.get_xticklabels(), rotation=90, ha='right')

plt.savefig(f"{cmp}_matrix.png",dpi=200)
plt.show()
plt.close()
