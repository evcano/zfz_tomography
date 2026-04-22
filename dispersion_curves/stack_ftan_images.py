import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
from glob import glob
from mpl_point_clicker import clicker
from scipy.io import loadmat
from scipy.interpolate import interp1d, CubicSpline
from scipy.signal import find_peaks



def pick_dispersion_curve(klicker, cper, vel_axis, A):
    positions = klicker.get_positions()

    for key in positions.keys():
        curve = positions[key]
        if not curve.any():
            continue

        # interpolate the picked curve into the FTAN grid
        imin = np.argmin(np.abs(cper - curve[:,0].min()))
        if cper[imin] < curve[:,0].min():
            imin += 1

        imax = np.argmin(np.abs(cper - curve[:,0].max()))
        if cper[imax] > curve[:,0].max():
            imax -= 1

        tmp_curve_cper= cper[imin:imax+1]
        ip1d = interp1d(x=curve[:,0], y=curve[:,1],
                        kind="linear", bounds_error=None)

        tmp_curve_gvel = ip1d(tmp_curve_cper)

        # find the closest ridge to the curve
        ridge_cper, ridge_gvel = find_ridges(cper, vel_axis, A)
        npoints = tmp_curve_gvel.size

        curve_cper = np.array([])
        curve_gvel = np.array([])

        for i in range(npoints):
            perdif = np.abs(ridge_cper - tmp_curve_cper[i])
            veldif = np.abs(ridge_gvel - tmp_curve_gvel[i])
            idx = np.argmin(perdif**2.0 + veldif**2.0)

            curve_cper = np.append(curve_cper, ridge_cper[idx])
            curve_gvel = np.append(curve_gvel, ridge_gvel[idx])

        _, idx = np.unique(curve_cper, return_index=True)
        curve_cper = curve_cper[idx]
        curve_gvel = curve_gvel[idx]
        curve = [curve_cper, curve_gvel]

    return curve


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

def plot_ftan(cper, vel_axis, A, curve=None):
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

    if curve:
        ax.plot(curve[0], curve[1], "-", c="k")

    ax.set_xlim(1.0, 21.0)
    return fig, ax


dpath = "./ftan_TT_checked"
veltype = "phase"
pick_curve = True

vaxis0 = np.arange(0.5, 5.002, 0.002)

# ---------------------------------------------------
obs_list = [f"ZF.OBS{x:02}" for x in range(1,13)]
obs_list.extend(["ZF.NORTH", "ZF.SOUTH"])

island_list = ["ZF.BREEM", "ZF.QUMAN"]

coast_list = ["ZF.LAVA","ZF.KHUF","SA.EWJHS"]

inland_list = ["KL.WEST","KL.EAST","KL.SOUTH","KL.DEEP"]
# ---------------------------------------------------
if veltype == "group":
    prefix1 = "GImg"
    prefix2 = "GroupVImg"
elif veltype == "phase":
    prefix1 = "PImg"
    prefix2 = "PhaseVImg"

for stkonly in ["obs","island","coast","inland"]:
    print(stkonly)

    files = glob(os.path.join(dpath,f"{prefix1}*"))
    files.sort()

    S = np.zeros((2251,37))

    for fname in files:
        fname2 = os.path.basename(fname)
        fname2 = fname2.replace(f"{prefix1}.","")
        fname2 = fname2.replace(".dat.mat","")
        s1, s2, _ = fname2.split("_")

        if stkonly == "obs":
            if s1 not in obs_list or s2 not in obs_list:
                continue
        elif stkonly == "island":
            if s1 not in island_list and s2 not in island_list:
                continue
        elif stkonly == "coast":
            if s1 not in coast_list and s2 not in coast_list:
                continue
        elif stkonly == "inland":
            if s1 not in inland_list and s2 not in inland_list:
                continue

        X = loadmat(fname)
        X = X[prefix2]
        for i in range(X.shape[1]):
            X[:,i] /= np.max(abs(X[:,i]))

        vaxis = loadmat(fname.replace(prefix1,"VPoint"))
        vaxis = vaxis["VPoint"].flatten()

        taxis = loadmat(fname.replace(prefix1,"TPoint"))
        taxis = taxis["TPoint"].flatten()

        if X.shape[0] != vaxis0.size:
            idx1 = np.argmin(abs(vaxis0 - vaxis.min()))
            idx2 = np.argmin(abs(vaxis0 - vaxis.max()))
            S[idx1:idx2+1,:] += X
        else:
            S += X


    for i in range(S.shape[1]):
        S[:,i] /= np.max(abs(S[:,i]))

    fig, ax = plot_ftan(taxis, vaxis, S)

#    if pick_curve:
#        klicker = clicker(ax, ["curve"], legend_loc="upper right")
#        plt.show()
#        plt.close()
#
#        curve = pick_dispersion_curve(klicker, taxis, vaxis, S)
#
#        if curve:
#            x = curve[0]
#            y = curve[1]
#
#            spl = interp1d(x, y)
#            x2 = np.arange(x.min(), x.max(), 0.1)
#            y2 = spl(x2)
#
#            fig, ax = plot_ftan(taxis, vaxis, S, curve=[x2,y2])
#            plt.show()
#            plt.close()
#
#            x = curve[0]
#
#            oname = os.path.join(dpath,f"mean-vel-{veltype}-{stkonly}.txt")
#            oarr = np.array([x2, y2]).T
#            np.savetxt(oname, oarr, fmt="%.4f")
#    else:
#        plt.show()
#        plt.close()

    if pick_curve:
        klicker = clicker(ax, ["up", "mid", "low"], legend_loc="upper right")
        plt.show()
        plt.close()

        positions = klicker.get_positions()
        for key in positions.keys():
            curve = positions[key]
            if not curve.any():
                continue

            x = curve[0]
            y = curve[1]

            iobj = interp1d(x, y, kind="linear")
            x2 = np.arange(x.min(), x.max(), 0.1)
            y2 = iobj(x2)

            plt.plot(x2,y2)
        plt.show()
