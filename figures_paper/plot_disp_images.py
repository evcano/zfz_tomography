import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
from glob import glob
from obspy.geodetics import gps2dist_azimuth
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.io import loadmat
from cmcrameri import cm


rc('text', usetex=True)


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


def segment_array(a):
    x = [s for s in np.ma.clump_unmasked(np.ma.masked_invalid(a))]
    return x


def read_average_velimage(impath,suffix):
    S_file = os.path.join(impath, f"S_{suffix}.npy")
    taxis_file = os.path.join(impath, f"taxis_{suffix}.npy")
    vaxis_file = os.path.join(impath, f"vaxis_{suffix}.npy")

    S = np.load(S_file)
    taxis = np.load(taxis_file)
    vaxis = np.load(vaxis_file)

    return S, taxis, vaxis


def read_velimage(impath, fname, vtype):
    if vtype == "group":
        imtype = "GImg"
        imtype2 = "GroupVImg"
    else:
        imtype = "PImg"
        imtype2 = "PhaseVImg"

    S_file = os.path.join(impath,f"{imtype}.{fname}.mat")
    taxis_file = os.path.join(impath,f"TPoint.{fname}.mat")
    vaxis_file = os.path.join(impath,f"VPoint.{fname}.mat")

    S = loadmat(S_file)
    S = S[imtype2]

    taxis = loadmat(taxis_file)
    taxis = taxis["TPoint"].flatten()

    vaxis = loadmat(vaxis_file)
    vaxis = vaxis["VPoint"].flatten()

    return S, taxis, vaxis


vtype = "group"

wdir = "../dispersion_curves"
avgimpath = f"{wdir}/average_disp_img"
refcurvepath = f"{wdir}/reference_curves"

ctypes = ["obs","island","coast"]
pairs = ["ZF.OBS06_ZF.OBS04","ZF.OBS07_ZF.QUMAN","ZF.OBS10_ZF.KHUF"]

# figure settings
figname = f"./dispersion_images_{vtype}.png"
figsize = (7., 7.2)
fs = 9
dpi = 400
xticks = np.arange(2,15,2)
yticks = np.arange(0,6,dtype=float)

# ---------------------------------------------------
fig = plt.figure(layout="constrained")
subfigs = fig.subfigures(2,1)

subfigs[0].suptitle(r"\textbf{Rayleigh waves}", x=0.1, fontsize=fs+2)
subfigs[1].suptitle(r"\textbf{Love waves}", x=0.1, fontsize=fs+2, y=1.1)

ax0 = subfigs[0].subplots(2,3)
ax1 = subfigs[1].subplots(2,3)

letters1 = [
    [r"\textbf{(a)}",r"\textbf{(b)}",r"\textbf{(c)}"],
    [r"\textbf{(d)}",r"\textbf{(e)}",r"\textbf{(f)}"],
           ]

letters2 = [
    [r"\textbf{(g)}",r"\textbf{(h)}",r"\textbf{(i)}"],
    [r"\textbf{(j)}",r"\textbf{(k)}",r"\textbf{(l)}"],
           ]

# loop over subfigures (Rayleigh and Love waves)
for k, ax in enumerate([ax0, ax1]):
    if k == 0:
        cmp = "ZZ"
        letters = letters1
        wave = "Rayleigh"
    elif k == 1:
        cmp = "TT"
        letters = letters2
        wave = "Love"
    # loop over subfigure row
    for j in range(0,2):
        # loop over subfigure columns
        for i in range(0,3):
            ctype = ctypes[i]
            pair = pairs[i]

            suffix = f"{cmp.lower()}-{vtype}-{ctype}"
            rcfile = f"{suffix}-twindow.txt"
            RC  = np.loadtxt(os.path.join(refcurvepath,rcfile))

            # plot dispersion curve boundaries
            if j == 0:
                S, taxis, vaxis = read_average_velimage(avgimpath,suffix)
                ax[j,i].contourf(taxis,vaxis,S,cmap=cm.lapaz,levels=100)

                ctype=ctype.replace("coast","land")
                ctype=ctype.capitalize()
                if ctype == "Obs":
                    ctype=ctype.upper()

                title = (letters[j][i]+" Stack of "
                         f"OBS-{ctype} station pairs")

            # plot examples of dispersion curves
            elif j == 1:
                impath = f"{wdir}/disp_img_{cmp}/{ctype}"
                fname = f"{pair}_{cmp}.dat"
                S, taxis, vaxis = read_velimage(impath,fname,vtype)
                CS = ax[j,i].contourf(taxis,vaxis,S,cmap=cm.lapaz,levels=100)

                curvepath = f"{wdir}/disp_cur_{cmp}/automatically_picked/{ctype}"
                if vtype == "group":
                    cfile = f"GDisp.{pair}_{cmp}.dat"
                else:
                    cfile = f"CDisp.T.{pair}_{cmp}.dat"

                # dispersion curve
                sta1coor = np.loadtxt(os.path.join(curvepath,cfile),
                                      max_rows=1)

                sta2coor = np.loadtxt(os.path.join(curvepath,cfile),
                                      max_rows=1, skiprows=1)

                C = np.loadtxt(os.path.join(curvepath,cfile),
                               skiprows=2)

                cper = C[:,0]
                oper = C[:,1]
                vel = C[:,2]
                snr = C[:,3]
                qlt0 = C[:,4]

                qlt = quality_control(oper,vel,snr,sta1coor,sta2coor,
                                      snrthr=6, nlambda=2)
                qlt[qlt==0] = np.nan
                qlt[qlt0==0] = np.nan
                seg_idx = segment_array(qlt)

                for ii, sidx in enumerate(seg_idx):
                    if ii == 0:
                        label = "Dispersion curve"
                    else:
                        label = None
                    ax[j,i].plot(cper[sidx],vel[sidx],c="k",
                                 label=label)

                titpair = pair
                titpair = titpair.replace("_","-")
                titpair = titpair.replace("ZF.","")
                title = (f"{letters[j][i]} "
                         f"{titpair}")

            # always plot boundaries
            if vtype == "group":
                ax[j,i].plot(RC[:,0],RC[:,2],"green",label="Bound")
                ax[j,i].plot(RC[:,0],RC[:,3],"green")
            else:
                ax[j,i].plot(RC[:,0],RC[:,1],"green",label="Ref. curve")

            # FIGURE ------------------------------------------------
            xlabel = "Period [s]"
            ylabel = f'Velocity [km/s]'
            xticklabel = None
            yticklabel = None

            if i > 0:
                ylabel = None
                yticklabel = []
            if j == 0:
                xlabel = None
                xticklabel = []

            # ticks
            ax[j,i].set_xticks(xticks,labels=xticklabel)
            ax[j,i].set_yticks(yticks,labels=yticklabel)
            ax[j,i].tick_params(labelsize=fs)

            # labels
            ax[j,i].set_xlabel(xlabel, fontsize=fs)
            ax[j,i].set_ylabel(ylabel, fontsize=fs)

            ax[j,i].set_title(title,fontsize=fs,c="k",x=0.01,ha="left")

            # lims
            ax[j,i].set_xlim(xticks[0],xticks[-1])
            ax[j,i].set_ylim(0.5,yticks[-1])

            # legend
            if k == 0 and j == 1 and i == 1:
                leg = ax[j,i].legend(ncols=2,
                                     fontsize=fs,
                                     loc="upper left",
                                     bbox_to_anchor=(-0.1,-1.4,1.3,1),
                                     edgecolor="k",
                                     fancybox=False,
                                     mode="expand",
                                    )

            # colorbar
            if k == 0 and j == 1 and i == 2:
                cbaxes = inset_axes(ax[j,i],
                                    width="50%",
                                    height="8%",
                                    loc="lower left",
                                    bbox_to_anchor=(0.2,-0.6, 1,1),
                                    bbox_transform=ax[j,i].transAxes,
                                   )
                cbar = fig.colorbar(CS,
                                    cax=cbaxes,
                                    orientation="horizontal",
                                    extend="both",
                                   )

                cbar.ax.set_xticks([S.min(),S.max()])
                cbar.ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
                cbar.ax.tick_params(labelsize=fs)
                cbar.ax.set_title("Amplitude",fontsize=fs,y=-2.5)

fig.set_size_inches(figsize)
fig.savefig(figname,dpi=dpi)
plt.close()
