import cartopy.crs as ccrs
import cartopy.mpl.geoaxes
import geopandas
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as pltticker
import matplotlib.tri as tr
import numpy as np
import utm
import os
import rasterio
import shapely
from cartopy import geodesic
from cmcrameri import cm
from glob import glob
from matplotlib.colors import LinearSegmentedColormap
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.patheffects import Stroke, Normal
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from plot_study_area import get_elevation_map, zebra_frame


def vel_cmap():
    colors = [
        (0.5, 0.0, 0.0),
        (1.0,1.0,0.3),
        (0.0,0.6,0.0),
        (0.1,0.9,1.0),
        (0.,0.,0.5),
    ]

    pts = [0.0,0.32,0.5,0.68,1.0]

    vcmap = LinearSegmentedColormap.from_list(
        'avel',
        list(zip(pts,colors)),
        N=30,
    )

    return vcmap


rc("text",usetex=True)

# files path
elvfile = "../bathymetry/gebco_2024_n32.0_s12.0_w30.0_e46.0.tif"
stafile = "./all_stations_coordinates.txt"

modfile = "./3d_model_figures/average_model.npy"
stdfile = "./3d_model_figures/std_model.npy"
latfile = "./3d_model_figures/lat_mesh.npy"
lonfile = "./3d_model_figures/lon_mesh.npy"
depfile = "./3d_model_figures/dep_mesh.npy"

eqfile = None #"../other_geophysical_data/2_relocated_earthquakes.csv"
gravfile = None #"../other_geophysical_data/grav_32.1.nc"  # free-air anomaly
smapfile = None #"../other_geophysical_data/MGDS_Download/test.dat" # s reflection map

# figure settings
figname = "gravity_sreflection.png"
figsize = (6.97,3.68)
fs = 9
fs2 = 6
dpi = 300

latmin = 24.
latmax = 26.
lonmin = 36.
lonmax = 37.2
mapnx = 200
mapny = 200

xticks = [36,37,37.2]
yticks = [24,25,26]
elvlevels = np.arange(-2,2.4,0.4)

std_thr = 0.5

cbar_per_title = "Vs [km/s]"
cbar_per_ticks = [1.5,2.5,3.5,4.5]
per_cont = np.arange(1,5,0.25)

cbar_gra_title = "Anomaly [mGal]"
cbar_gra_ticks = [-100,0,100]
gra_cont = np.arange(-100,120,20)

sref_cont = np.arange(500,2200,100)
cbar_sref_title = "Depth [m]"
cbar_sref_ticks = [500, 1250, 2000]

# ------------------------------------------------
# load stations
stanet = np.loadtxt(stafile, skiprows=1, usecols=0, dtype="U20")
staname = np.loadtxt(stafile, skiprows=1, usecols=1, dtype="U20")
stalat = np.loadtxt(stafile, skiprows=1, usecols=2)
stalon = np.loadtxt(stafile, skiprows=1, usecols=3)

# read gravity model
if gravfile:
    h5obj = h5py.File(gravfile, "r")
    gravlat = np.array(h5obj["lat"])
    gravlon = np.array(h5obj["lon"])
    gravz = np.array(h5obj["z"])
    gravz /= 10.0  # miligals

    i0 = np.where((gravlat >= latmin) & (gravlat <= latmax))[0]
    i1 = np.where((gravlon >= lonmin) & (gravlon <= lonmax))[0]

    gravlat = gravlat[i0[0]:i0[-1]]
    gravlon = gravlon[i1[0]:i1[-1]]
    gravz = gravz[i0[0]:i0[-1], i1[0]:i1[-1]]

# read s-reflection map
if smapfile:
    smap = np.loadtxt(smapfile)
    idx = np.where(smap[:,1]>=24)
    smap_lon = smap[idx,0].flatten()
    smap_lat = smap[idx,1].flatten()
    smap = smap[idx,2].flatten()
    smap *= -1         

    # s-reflection triangulation
    triang = tr.Triangulation(smap_lon,smap_lat)
    isbad = np.isnan(smap)                                                                                                    
    mask = np.all(np.where(isbad[triang.triangles],True,False),axis=1)
    triang.set_mask(mask)
    smap[np.isnan(smap)] = 0.0

# read model
VS = np.load(modfile)
STD = np.load(modfile.replace("average","std"))

LATMESH = np.load(latfile)
LONMESH = np.load(lonfile)
DEPMESH = np.load(depfile)

# mute water
STD[np.where(VS == 0.0)] = np.nan
VS[np.where(VS == 0.0)] = np.nan

# mute high std regions
VS[np.where(STD > std_thr)] = np.nan

# read elevation
dataset = rasterio.open(elvfile)
gebco_elevation =dataset.read(1)

extent_map = [lonmin, lonmax, latmin, latmax]
maplatax = np.linspace(latmin, latmax, mapny)
maplonax = np.linspace(lonmin, lonmax, mapnx)

elevation ,shaded, _, _ = get_elevation_map(
    maplatax, maplonax, extent_map, dataset, gebco_elevation)

if eqfile:
    eqlon = np.loadtxt(eqfile, delimiter=",", skiprows=1, usecols=7)
    eqlat = np.loadtxt(eqfile, delimiter=",", skiprows=1, usecols=8)
    eqdep = np.loadtxt(eqfile, delimiter=",", skiprows=1, usecols=9)  #km
    eqdep = -eqdep

# FIGURE
# ------------------------------------------------
setattr(cartopy.mpl.geoaxes.GeoAxes, 'zebra_frame', zebra_frame)
mercator = ccrs.PlateCarree()
fig = plt.gcf()

for fc in range(3):
    ax = plt.subplot(int(f"13{fc+1}"),projection=mercator)
    ax.set_extent(extent_map)

    if fc == 0:
        if not gravfile:
            continue

        # gravity map
        cbar_title = cbar_gra_title
        cbar_ticks = cbar_gra_ticks
        vmin = cbar_ticks[0]
        vmax = cbar_ticks[-1]

        ax.set_title(r"\textbf{(a)} Free-air gravity anomaly",
                     loc="left",
                     fontsize=fs,
        )

        CS = ax.contourf(gravlon, gravlat, gravz,
                         levels=np.arange(-100,110,10),
                         vmin=vmin,
                         vmax=vmax,
                         cmap=cm.roma,
                         transform=mercator,
                         extend="both",
                        )

        ECS = ax.contour(gravlon,gravlat,gravz,
                         levels=gra_cont,
                         colors="w",
                         linewidths=0.3,
                         linestyles="solid",
                         alpha=0.9
                        )

        ax.clabel(ECS,ECS.levels,inline=True,fontsize=fs2)

    elif fc == 1:
        # vs map
        cbar_title = cbar_per_title
        cbar_ticks = cbar_per_ticks
        vmin = cbar_ticks[0]
        vmax = cbar_ticks[-1]

        idx = np.argmin(abs(DEPMESH[0,0,:]-3.))
        modslice = VS[:,:,idx]

        # velocity slice
        cmap = vel_cmap()

        ax.set_title(r"\textbf{(b)} Shear-wave velocity",
                     loc="left",
                     fontsize=fs,
        )

        CS = ax.contourf(LONMESH[:,:,0],LATMESH[:,:,0],modslice,
                         levels=np.linspace(vmin,vmax,100),
                         cmap=cmap,
                         transform=mercator,
                         extend="both",
                        )

        ECS = ax.contour(LONMESH[:,:,0],LATMESH[:,:,0],modslice,
                         levels=per_cont,
                         colors="w",
                         linewidths=0.3,
                         linestyles="solid",
                         alpha=0.9
                        )

        ax.clabel(ECS,ECS.levels,inline=True,fontsize=fs2)
        ax.text(36.75,25.8,"3 km depth",fontsize=fs)

    elif fc == 2:
        if not smapfile:
            continue

        # s reflection
        ax.set_title(r"\textbf{(c)} S reflector",
                     loc="left",
                     fontsize=fs,
        )

        cbar_title = cbar_sref_title
        cbar_ticks = cbar_sref_ticks
        vmin = cbar_ticks[0]
        vmax = cbar_ticks[-1]

        pts = [0.,0.25,0.5,0.75,1.0]
        col = ["white","red","yellow","cyan","blue"]
        cmap = LinearSegmentedColormap.from_list("my",list(zip(pts,col)),N=200)


        CS = ax.tricontourf(triang, smap, levels=sref_cont,
                       cmap=cmap, extend="both")
                                   
        ECS = ax.tricontour(triang, smap,
                            levels=np.arange(900,2100,200),
                            colors="w",linewidths=0.3,
                            linestyles="solid",alpha=0.9)

        ax.clabel(ECS,ECS.levels,inline=True,fontsize=fs2)

    # coast line
    ax.contour(maplonax,maplatax,elevation,
               levels=[0],
               colors="k",
               linewidths=0.2,
               transform=mercator,
              )

    # Mabahis Mons
    ax.scatter(36.03, 25.47, s=11.0, marker="^", c="k")
    ax.text(36.15, 25.46, "M.M.",
            ha="center", fontsize=fs2,zorder=20000)

    # Mabahis Deep
    ax.scatter(36.11, 25.39, s=11.0, marker="o", c="k")
    ax.text(36.2, 25.27, "M.D.",
            ha="center", fontsize=fs2,zorder=20000)

    # Kebrit Deep
    ax.scatter(36.26, 24.715, s=11.0, marker="o", c="k")
    ax.text(36.3, 24.6, "K.D.",
            ha="center", fontsize=fs2,zorder=20000)

    # eq
    if eqfile:
        zz = -3
        idx = np.where((eqdep>zz-0.5)&(eqdep<zz+0.5))
        ax.scatter(eqlon[idx],eqlat[idx],marker="*",c="m",s=1)

    #  colorbar
    cbaxes = inset_axes(ax,
                        width="40%",
                        height="2%",
                        loc="lower left",
                        bbox_to_anchor=(0.05,0.88,1,1),
                        bbox_transform=ax.transAxes,
                       )

    cbar = fig.colorbar(CS,
                        cax=cbaxes,
                        orientation="horizontal",
                        extend="both",
                       )

    cbar.ax.set_title(cbar_title,loc="center",fontsize=fs)
    cbar.ax.set_xticks(cbar_ticks)
    cbar.ax.tick_params(labelsize=fs)

    # north arrow
    arrowx, arrowy, arrowlen = 0.05, 0.2, 0.08
    ax.annotate("N",
                xy=(arrowx,arrowy),
                xytext=(arrowx, arrowy-arrowlen),
                arrowprops=dict(facecolor='black', width=0.5,
                                headwidth=6, headlength=5),
                ha="center",
                va="center",
                fontsize=fs,
                xycoords=ax.transAxes,
               )

    # scalebar
    path2 = shapely.geometry.LineString([(latmin ,lonmin),(latmin,lonmin+1.0)])
    distmeter = geodesic.Geodesic().geometry_length(path2)
    scalebar = ScaleBar(distmeter,
                        units="m",
                        dimension="si-length",
                        fixed_value=50,
                        fixed_units="km",
                        box_alpha=0.0,
                        location="lower left",
                        scale_loc="bottom",
                        font_properties={"size":fs},
                       )

    ax.add_artist(scalebar)

    # ticks and labels
    ax.set_xticks(xticks, crs=mercator)
    ax.xaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°",sep=""))
    ax.xaxis.set_minor_locator(pltticker.MultipleLocator(base=0.1))

    ax.set_yticks(yticks, crs=mercator)
    ax.yaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°",sep=""))
    ax.yaxis.set_minor_locator(pltticker.MultipleLocator(base=0.1))

    xticklabels = ax.get_xticklabels()
    xticklabels[-1] = ""
    ax.set_xticklabels(xticklabels)

    lleft = True
    if fc > 0:
        lleft = False
    ax.tick_params(labelsize=fs,
                   labelbottom=True,labeltop=False,labelleft=lleft,labelright=False,
                   bottom=True,top=False,left=True,right=False)

    # zebra frame
    ax.zebra_frame(lw=3,crs=mercator,iFlag_outer_frame_in=None)

# savefig
fig.set_size_inches(figsize)
fig.tight_layout()
fig.savefig(figname, dpi=dpi)
plt.close()
