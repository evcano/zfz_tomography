import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
from glob import glob
from mpl_point_clicker import clicker
from scipy.io import loadmat
from scipy.interpolate import interp1d, CubicSpline
from scipy.signal import find_peaks



def save_curves(path, prefix, curves):
    # reference curve
    for curve in curves:
        k = curve[2]
        if k == "low":
            low = curve[1]
        elif k == "mid":
            ref = curve[1]
        elif k == "upp":
            upp = curve[1]
        per = curve[0]

    var = (upp - ref) + (ref - low)
    var *= 0.5

    out_twindow = np.array([per, ref, low, upp]).T
    oname = os.path.join(path,f"{prefix}-twindow.txt")
    np.savetxt(oname, out_twindow, fmt="%.5f")

    out_refcurve = np.array([per, ref, var]).T
    oname = os.path.join(path,f"{prefix}-refcurve.txt")
    np.savetxt(oname, out_refcurve, fmt="%.5f")

    return

def pick_ref_curves(klicker, cper, vel_axis, A):
    positions = klicker.get_positions()
    curves = []

    pmin = 0
    pmax = 100
    for key in positions.keys():
        curve = positions[key]
        x = curve[:,0]
        pmin = max(pmin, x.min())
        pmax = min(pmax, x.max())

    imin = np.argmin(np.abs(cper - pmin))
    imax = np.argmin(np.abs(cper - pmax))
    x = cper[imin:imax+1]

    for key in positions.keys():
        curve = positions[key]
        if not curve.any():
            continue

        # interpolate the picked curve into the FTAN grid
        ip1d = interp1d(x=curve[:,0], y=curve[:,1],
                        kind="linear", bounds_error=None,
                        fill_value="extrapolate")

        y = ip1d(x)

        curves.append([x, y, key])

    return curves


def find_ridges(cper, vel_axis, A):
    nfilt = A.shape[1]
    ridge_cper = []
    ridge_gvel = []
    for i in range(nfilt):
        peaks, properties = find_peaks(A[:,i])
        for p in peaks:
            ridge_cper.append(cper[i])
            ridge_gvel.append(vel_axis[p])

    return ridge_cper, ridge_gvel


def plot_ftan(cper, vel_axis, A, curves=None):
    fig, ax = plt.subplots(1, 1, figsize=(10,5))

    # FTAN diagram
    ax.contourf(cper, vel_axis, A, cmap="magma", levels=100)

    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

    ax.set_xlabel("period [s]")
    ax.set_ylabel("velocity [km/s]")

    # local maxima
    ridges_per, ridges_vel = find_ridges(cper, vel_axis, A)
    ax.scatter(ridges_per, ridges_vel, c="g", s=1.0)

    if curves:
        for curve in curves:
            x = curve[0]
            y = curve[1]
            k = curve[2]
            ax.plot(x, y, "-", c="k", label=k)
        ax.legend()

    ax.set_xlim(1.0, 21.0)

    return fig, ax

suffix = "zz-group-coast"

dpath = "./disp_img_TT/obs"
corrfile = "./TT/good_corr_obs.txt"

veltype = "group"
pick_curve = False

curveprefix = ""
curvepath = "./"
outpath = "./average_disp_img"

# ---------------------------------------------------
vaxis0 = np.arange(0.5, 5.002, 0.002)

if veltype == "group":
    prefix1 = "GImg"
    prefix2 = "GroupVImg"
elif veltype == "phase":
    prefix1 = "PImg"
    prefix2 = "PhaseVImg"

files = np.loadtxt(corrfile, dtype=str)
files.sort()

S = np.zeros((2251,75))

for f in files:
    fname = os.path.join(dpath, f"{prefix1}.{f}.mat")

    # disp image
    X = loadmat(fname)
    X = X[prefix2]
    for i in range(X.shape[1]):
        X[:,i] /= np.max(abs(X[:,i]))

    # vel axis
    vaxis = loadmat(fname.replace(prefix1,"VPoint"))
    vaxis = vaxis["VPoint"].flatten()

    # time axis
    taxis = loadmat(fname.replace(prefix1,"TPoint"))
    taxis = taxis["TPoint"].flatten()

    # stack
    if X.shape[0] != vaxis0.size:
        idx1 = np.argmin(abs(vaxis0 - vaxis.min()))
        idx2 = np.argmin(abs(vaxis0 - vaxis.max()))
        S[idx1:idx2+1,:] += X
    else:
        S += X

# normalize
for i in range(S.shape[1]):
    S[:,i] /= np.max(abs(S[:,i]))

## save average velocity images
#np.save(os.path.join(outpath,f"S_{suffix}"),S)
#np.save(os.path.join(outpath,f"taxis_{suffix}"),taxis)
#np.save(os.path.join(outpath,f"vaxis_{suffix}"),vaxis)

fig, ax = plot_ftan(taxis, vaxis, S)

if pick_curve:
    klicker = clicker(ax, ["upp","mid","low"], legend_loc="upper right")
    plt.show()
    plt.close()

    curves = pick_ref_curves(klicker, taxis, vaxis, S)

    if curves:
        fig, ax = plot_ftan(taxis, vaxis, S, curves=curves)
        plt.show()
        plt.close()

        oname = os.path.join(dpath,f"mean-vel-{veltype}.txt")
        save_curves(curvepath, curveprefix, curves)
else:
    plt.show()
    plt.close()
