import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
from scipy.interpolate import interp1d
from scipy.io import loadmat
from obspy.geodetics import gps2dist_azimuth
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector


class SelectFromCollection:
    """
    Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        Axes to interact with.
    collection : `matplotlib.collections.Collection` subclass
        Collection you want to select from.
    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to *alpha_other*.
    """

    def __init__(self, ax, collection, alpha_other=0.1):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


def segment_array(a):
    x = [s for s in np.ma.clump_unmasked(np.ma.masked_invalid(a))]
    return x


def quality_control(per, vel, snr, s1coor, s2coor, snrthr=0.0, nlambda=2):
    qlt = np.ones(per.size)

    # check progapaged wavelength
    stadis, _, _ = gps2dist_azimuth(lon1=s1coor[0],
                                    lat1=s1coor[1],
                                    lon2=s2coor[0],
                                    lat2=s2coor[1]
                                   )
    stadis /= 1000.0  # km

    wavelength = per * vel
    idx = np.where(wavelength * nlambda > stadis)
    qlt[idx] = 0.0

    # check snr
    idx = np.where(snr < snrthr)
    qlt[idx] = 0.0

    return qlt


imgpath = "./disp_img_ZZ/coast"
curvepath = "./disp_cur_ZZ/automatically_picked/coast"
outpath = "./disp_cur_ZZ/manually_checked/coast"

veltype = "group"

snrthr = 6
nlambda = 2
minpts = 4

# ---------------------------------------------------
if veltype == "group":
    prefix1 = "GDisp"
    prefix2 = "GImg"
    prefix3 = "GroupVImg"
elif veltype == "phase":
    prefix1 = "CDisp.T"
    prefix2 = "PImg"
    prefix3 = "PhaseVImg"

files = glob(os.path.join(curvepath, f"{prefix1}*"))
files.sort()

# loop over curves
for fim in files:
    fname  = os.path.basename(fim)
    fname = fname.replace(f"{prefix1}.", "")

    # read image
    fname = f"{prefix2}.{fname}.mat"
    A = loadmat(os.path.join(imgpath, fname))
    A = A[prefix3]

    vaxis = loadmat(os.path.join(imgpath, fname.replace(prefix2,"VPoint")))
    vaxis = vaxis["VPoint"].flatten()

    taxis = loadmat(os.path.join(imgpath, fname.replace(prefix2,"TPoint")))
    taxis = taxis["TPoint"].flatten()

    # prepare plot
    fig, ax = plt.subplots(figsize=(15,15))
    ax.set_title(fname)
    ax.contourf(taxis, vaxis, A, cmap="jet", levels=100)

    # read station coordinates
    try:
        # lon, lat
        sta1coor = np.loadtxt(fim, max_rows=1)
        sta2coor = np.loadtxt(fim, skiprows=1, max_rows=1)
    except:
        continue

    # read curve
    try:
        X = np.loadtxt(fim, skiprows=2)
    except:
        continue

    cper = X[:,0]  # filter central period
    oper = X[:,1]  # observed period
    vel = X[:,2]   # velocity curve
    snr = X[:,3]   # spectral snr
    qlt1 = X[:,4]  # quality determined by automatic code, 0 bad, 1 good

    # check SNR and maximum solved period again
    # this time, we check the maximum solved period using the observed period
    qlt2 = quality_control(oper, vel, snr, sta1coor, sta2coor,
                           snrthr=snrthr, nlambda=nlambda)

    # keep good periods
    qlt2[qlt2==0] = np.nan

    # get good curve segments
    slices_list = segment_array(qlt2)

    # loop over segments
    vel_itp = np.zeros(oper.size)
    snr_itp = np.zeros(oper.size)

    for slice_obj in slices_list:
        x = oper[slice_obj]
        y = vel[slice_obj]
        z = snr[slice_obj]

        # discard short segments
        if np.sum(~np.isnan(x)) < minpts:
            continue

        # using the observed period, interpolate velocity curve and
        # spectral snr to a uniform period grid (filter central periods)
        iobj = interp1d(x, y, kind="cubic", bounds_error=False, fill_value=0.0)
        y = iobj(cper)

        iobj = interp1d(x, z, kind="cubic", bounds_error=False, fill_value=0.0)
        z = iobj(cper)

        vel_itp += y
        snr_itp += z

    pts = ax.scatter(cper, vel_itp, c="k")
    selector = SelectFromCollection(ax, pts)

    def accept(event):
        if event.key == "enter":
            print("Selected points:")
            idx = selector.ind

            per_out = cper[idx]
            vel_out = vel_itp[idx]
            snr_out = snr_itp[idx]
            outarr = np.array([per_out, vel_out, snr_out]).T

            ofile = fname
            ofile = ofile.replace(prefix2, prefix1)
            ofile = ofile.replace(".mat","")
            ofile = os.path.join(outpath, ofile)

            np.savetxt(ofile, outarr, fmt="%.4f")
            print(f"{ofile} saved")

            selector.disconnect()
            fig.canvas.draw()
            plt.close()

    fig.canvas.mpl_connect("key_press_event", accept)
    plt.show()
