import cartopy.crs as ccrs
import cartopy.mpl.geoaxes
import geopandas
import matplotlib.pyplot as plt
import matplotlib.ticker as pltticker
import numpy as np
import os
import rasterio
import pyvista as pv
import shapely

from cartopy import geodesic
from cmcrameri import cm
from glob import glob
from matplotlib.colors import LightSource,LinearSegmentedColormap
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.patheffects import Stroke, Normal
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from obspy.geodetics.base import gps2dist_azimuth

from plot_study_area import get_elevation_map, zebra_frame


rc("text",usetex=True)


def vel_cmap():
    colors = [
            (0.5, 0.0, 0.0),    # Dark red (fast)
            (1.0, 1.0, 0.3),   # Yellow
            (0., 0.6, 0.),   # Green
            (0.1, 0.9, 1.0),   # Cyan
            (0., 0., 0.5),   # Dark blue (slow)
    ]

    pts = [0., 0.32, 0.5, 0.68, 1.0]

    vcmap = LinearSegmentedColormap.from_list('avel',
                                              list(zip(pts,colors)),
                                              N=30,
                                             )
    return vcmap


def slice_to_array(slc, normal, origin, name, xlim, ylim, nx=50, ny=50):
    """Converts a PolyData slice to a 2D NumPy array.

    It is crucial to have the true normal and origin of
    the slicing plane

    Parameters
    ----------
    slc : PolyData
        The slice to convert.
    normal : tuple(float)
        the normal of the original slice
    origin : tuple(float)
        the origin of the original slice
    name : str
        The scalar array to fetch from the slice
    ni : int
        The resolution of the array in the i-direction
    nj : int
        The resolution of the array in the j-direction

    """
    # Make structured grid
    x = np.linspace(xlim[0], xlim[1], nx)
    y = np.linspace(ylim[0], ylim[1], ny)
    z = np.array([0])
    plane = pv.StructuredGrid(*np.meshgrid(x,y,z))

    # rotate and translate grid to be ontop of the slice
    direction = normal / np.linalg.norm(normal)
    vx = np.array([0., 0., 1.])
    vx -= vx.dot(direction) * direction
    vx /= np.linalg.norm(vx)
    vy = np.cross(direction, vx)
    rmtx = np.array([vx, vy, direction])

    plane.points = plane.points.dot(rmtx)
    plane.points -= plane.center
    plane.points += origin

    # resample the data
    sampled = plane.sample(slc,)

    # fill bad data
    idx = ~sampled["vtkValidPointMask"].view(bool)
    sampled[name][idx] = np.nan

    # get arrays
    array = sampled[name].reshape(sampled.dimensions[0:2])

    # get points coordinates
    pcoor = sampled.points
    arrayX = pcoor[:,0].reshape(array.shape)
    arrayY = pcoor[:,1].reshape(array.shape)
    arrayZ = pcoor[:,2].reshape(array.shape)

    return array, arrayX, arrayY, arrayZ


# paths
stafile = "../stations_coordinates.txt"
elvfile = "../../bathymetry/gebco_2024_n32.0_s12.0_w30.0_e46.0.tif"
eqfile = None #"../../other_geophysical_data/2_relocated_earthquakes.csv"

modfile = "./average_model.npy"
stdfile = "./std_model.npy"
latfile = "./lat_mesh.npy"
lonfile = "./lon_mesh.npy"
depfile = "./dep_mesh.npy"

# vertical slices coordinates
mtype = "vel"

vslices = [
    ["v1","A", (1, 0.45, 0.0), (36.35,25,0.0)],
    ["v2","B", (1, 0.09, 0.0), (36.44,25.,0.0)],
    ["v3","C", (-0.5, 0.9, 0.0), (36.45,25.55,0.0)],
    ["v4","D", (-0.2, 1., 0.0), (36.5,24.765,0.0)],
    ["v5","E", (0.15, 1., 0.0), (36.6,24.25,0.0)],
]

zslices = [["h1",3],["h2",5],["h3",9]]

# map settings
utmno = 37
utmle = "R"

latmin = 24.
latmax = 26.
lonmin = 36.
lonmax = 37.2

mapnx = 200
mapny = 200

lonticks = [36,37,37.2]
latticks = [24,25,26]

# figure settings
figname = f"./slices_paper_{mtype}.png"

figsize = (7.7,9.6)
fs = 9
fs2 = 8
dpi = 600

vscale = 0.05
zmax = -12.2
std_thr = 0.5

cbar_vel = [1.5,2.5,3.5,4.5]
cbar_per = [-25,-15,0,15,25]
cbar_std = [0,0.5,1.0]
xlim = (0,250)

subtit = [
    r"\textbf{(d)}",
    r"\textbf{(e)}",
    r"\textbf{(f)}",
    r"\textbf{(g)}",
    r"\textbf{(h)}",
    r"\textbf{(a)}",
    r"\textbf{(b)}",
    r"\textbf{(c)}",
]

# ------------------------------------------------
if mtype == "std":
    eqfile = None

# read hypocenters
if eqfile:
    eqlon = np.loadtxt(eqfile, delimiter=",", skiprows=1, usecols=7)
    eqlat = np.loadtxt(eqfile, delimiter=",", skiprows=1, usecols=8)
    eqdep = np.loadtxt(eqfile, delimiter=",", skiprows=1, usecols=9)  #km
    eqdep = -eqdep

# read stations 
stanet = np.loadtxt(stafile, skiprows=1, usecols=0, dtype="U20")
staname = np.loadtxt(stafile, skiprows=1, usecols=1, dtype="U20")
stalat = np.loadtxt(stafile, skiprows=1, usecols=2)
stalon = np.loadtxt(stafile, skiprows=1, usecols=3)

# ------------------------------------------------
# read model
VS = np.load(modfile)
STD = np.load(stdfile)
LATMESH = np.load(latfile)
LONMESH = np.load(lonfile)
DEPMESH = np.load(depfile)
DEPMESH = -DEPMESH  # z decreases with depth

# cut to X km depth
zax = DEPMESH[0,0,:].copy()
idx = np.argwhere(zax >= zmax)
idx = idx.max()

VS = VS[:,:,:idx]
STD = STD[:,:,:idx]
LATMESH = LATMESH[:,:,:idx]
LONMESH = LONMESH[:,:,:idx]
DEPMESH = DEPMESH[:,:,:idx]

# scale depth
DEPMESH *= vscale

# mute water
STD[np.where(VS==0.0)] = np.nan
VS[np.where(VS==0.0)] = np.nan

# mute high std
VS[np.where(STD > std_thr)] = np.nan

# ------------------------------------------------
# read elevation
dataset = rasterio.open(elvfile)
gebco_elevation =dataset.read(1)

extent_map = [lonmin, lonmax, latmin, latmax]
latax = np.linspace(latmin, latmax, mapny)
lonax = np.linspace(lonmin, lonmax, mapnx)

elevation ,shaded, cmaptopo, cmaptopo_norm = get_elevation_map(
    latax, lonax, extent_map, dataset, gebco_elevation)

LONELV, LATELV = np.meshgrid(lonax, latax)

# scale elevation
elevation *= vscale

# ------------------------------------------------
# model colorbar
if mtype == "vel":
    MOD = VS
    cbar_ticks = cbar_vel
    cbar_title = "Vs [km/s]"
    cmap = vel_cmap()

elif mtype == "std":
    MOD = STD
    cmap = cm.oslo
    cbar_ticks = cbar_std
    cbar_title = "$\sigma$ [km/s]"

vmin = cbar_ticks[0]
vmax = cbar_ticks[-1]

# setup figure
setattr(cartopy.mpl.geoaxes.GeoAxes, 'zebra_frame', zebra_frame)

mercator = ccrs.PlateCarree()

fig, ax_all = plt.subplot_mosaic(
    [
     ["h1","h1","h2","h2","h3","h3"],
     [".","v1","v1","v1","v1","."],
     [".","v2","v2","v2","v2","."],
     ["v3","v3","v4","v4","v5","v5"],
    ],
    figsize=figsize,
    gridspec_kw={
        "height_ratios":[3.0, 1, 1, 1],
    },
    per_subplot_kw={
        "h1":{"projection": mercator},
        "h2":{"projection": mercator},
        "h3":{"projection": mercator},
    },
    constrained_layout=False,
)

# VERTICAL SLICES
nslices = len(vslices)

subcount = -1
for slc in vslices:
    subcount += 1

    ax = ax_all[slc[0]]
    scode = slc[1]
    normal = slc[2]
    origin = slc[3]

    # pyvista datasets
    model_pv = pv.StructuredGrid(LONMESH,LATMESH,DEPMESH)
    model_pv.point_data["values"] = MOD.flatten(order="F")
    if mtype != "std":
        model_pv.hide_points(np.argwhere(STD.flatten(order="F")>std_thr))

    topo_pv = pv.StructuredGrid(LONELV,LATELV,elevation)
    topo_pv.point_data["values"] = elevation.flatten(order="F")

    # pyvista slice
    model_slice = model_pv.slice(normal=normal,origin=origin)
    topo_slice = topo_pv.slice(normal=normal,origin=origin)

    # slice to numpy array
    arr, arrX, arrY, arrZ = slice_to_array(
        model_slice,
        normal=normal,
        origin=origin,
        name="values",
        xlim=(lonmin,lonmax),
        ylim=(latmin,latmax),
        nx=800,
        ny=800,
    )

    arrZ /= vscale # rescale z axis

    # delete nan values
    delrow = []
    for i in range(0, arr.shape[0]):
        if np.all(np.isnan(arr[i,:])):
            delrow.append(i)

    delcol = []
    for i in range(0, arr.shape[1]):
        if np.all(np.isnan(arr[:,i])):
            delcol.append(i)

    arr = np.delete(arr, delrow, axis=0)
    arrX = np.delete(arrX, delrow, axis=0)
    arrY = np.delete(arrY, delrow, axis=0)
    arrZ = np.delete(arrZ, delrow, axis=0)

    arr = np.delete(arr, delcol, axis=1)
    arrX = np.delete(arrX, delcol, axis=1)
    arrY = np.delete(arrY, delcol, axis=1)
    arrZ = np.delete(arrZ, delcol, axis=1)

    # make slices shorter
    tmpdic = {
    "A": [20,800],
    "B": [120,800],
    "C": [5,800],
    "D": [30,450],
    "E": [70,400],
    }

    if mtype == "std":
        tmpdic["C"] = [5,370]

    p1 = tmpdic[scode][0]
    p2 = tmpdic[scode][1]
    arr  = arr[:,p1:p2]
    arrX = arrX[:,p1:p2]
    arrY = arrY[:,p1:p2]
    arrZ = arrZ[:,p1:p2]
    
    # get slice end points
    lon1 = arrX[0,0]
    lon2 = arrX[0,-1]
    lat1 = arrY[0,0]
    lat2 = arrY[0,-1]
    slc.append([lon1,lat1,lon2,lat2])

    # calculate distance
    dist, _, _ = gps2dist_azimuth(lat1,lon1,lat2,lon2)  # in m
    dist /= 1000  # in km

    # slice axis
    x = np.linspace(0, dist, arr.shape[1])
    y = np.linspace(arrZ.min(), arrZ.max(), arr.shape[0])

    # topo slice
    # ----------------
    topo_x = topo_slice.points[:,0]
    topo_y = topo_slice.points[:,1]
    topo_z = topo_slice.points[:,2] / vscale

    topo_lon1 = topo_x[0]
    topo_lon2 = topo_x[-1]
    topo_lat1 = topo_y[0]
    topo_lat2 = topo_y[-1]

    topo_dist, _, _ = gps2dist_azimuth(topo_lat1,topo_lon1,topo_lat2,topo_lon2)  # in m
    topo_dist /= 1000  # in km
    topo_xax = np.linspace(0,topo_dist,topo_x.shape[0])

    off,_,_ = gps2dist_azimuth(topo_lat1,topo_lon1,lat1,lon1)
    off /= 1000  # km
    topo_xax -= off

    # project earthquakes
    if eqfile:
        eqdis = np.zeros(eqlon.shape)
        for e in range(0, eqlon.size):
            x1 = eqlon[e] - origin[0]
            y1 = eqlat[e] - origin[1]
            dd = x1*normal[0] + y1*normal[1]
            # point projected to the plane
            x2 = eqlon[e] - dd*normal[0]
            y2 = eqlat[e] - dd*normal[1]
            # distance of projection along slice
            eqdis[e],_,_ = gps2dist_azimuth(y2,x2,lat1,lon1)
            eqdis[e] /= 1000

            # distance between hypocenter and projection to plane
            pdis,_,_ = gps2dist_azimuth(y2,x2,eqlat[e],eqlon[e])
            pdis /= 1000.0
            if pdis > 2.0:
                eqdis[e] = np.nan

    # FIGURE
    # ---------------
    ax.plot(topo_xax,topo_z,"k",linewidth=0.7)
    im = ax.pcolormesh(x, y, arr,
                       cmap=cmap, vmin=vmin, vmax=vmax)

    # model contours
    contours = [1.2,2,2.5,3,3.5,4,4.5]
    ECS = ax.contour(x,y,arr,
                     levels=contours,
                     colors="white",
                     linewidths=0.2,
                     alpha=0.8,
                    )

    # model contours
    ax.clabel(ECS,
              ECS.levels,
              inline=True,
              fontsize=fs2-1,
             )

    # hypocenters
    if eqfile:
        ax.scatter(eqdis,eqdep,marker="*",c="m",s=0.15,alpha=0.8)

    #  colorbar
    if scode == "A":
        cbaxes = inset_axes(ax,
                            width="40%",
                            height="95%",
                            loc="lower right",
                            bbox_to_anchor=(1.0, -1.2, 0.1, 1.8),
                            bbox_transform=ax.transAxes,
                           )

        
        cbar = fig.colorbar(im,
                            cax=cbaxes,
                            orientation="vertical",
                            extend="both",
                           )

        cbar.set_label(cbar_title, fontsize=fs-1)
        cbar.ax.tick_params(labelsize=fs-1)

    # limits
    ax.set_xlim(-5,dist+10)
    ax.set_ylim(-12, 1.5)
    ax.tick_params(labelsize=fs)

    ax.text(-5.0, 2.5, f"{subtit[subcount]}",
            fontsize=fs)

    ax.text(0.0, 0.1, f"{scode}",
            fontsize=fs)

    ax.text(dist, 0.1, f"{scode}'",
            fontsize=fs)

    if scode in ["A","B","C"]:
        ax.set_ylabel("Depth [km]",fontsize=fs)
    else:
        ax.set_yticklabels([])

    if scode in ["C","D","E"]:
        ax.set_xlabel("Distance [km]",fontsize=fs)
    ax.tick_params(labelbottom=False)
    ax.tick_params(labelbottom=True)

    ax.xaxis.set_major_locator(pltticker.MultipleLocator(base=50))
    ax.xaxis.set_minor_locator(pltticker.MultipleLocator(base=10))
    ax.yaxis.set_minor_locator(pltticker.MultipleLocator(base=1))

# HORIZONTAL SLICES
elevation /= vscale
DEPMESH /= vscale
DEPMESH = -DEPMESH
zaxis = DEPMESH[0,0,:]
if eqfile:
    eqdep = -eqdep

for slc in zslices:
    subcount += 1

    ax = ax_all[slc[0]]
    ax.set_extent(extent_map)

    z = slc[1]
    idx = np.argmin(abs(zaxis-z))
    mod_slice = MOD[:,:,idx]

    # model
    CS = ax.contourf(LONMESH[:,:,0], LATMESH[:,:,0], mod_slice,
                     levels=np.linspace(vmin,vmax,100),
                     cmap=cmap,
                     transform=mercator,
                     extend="both",
                    )

    # model contours
    ECS = ax.contour(LONMESH[:,:,0],LATMESH[:,:,0],mod_slice,
                     levels=np.arange(1,5,0.2),
                     colors="white",
                     linewidths=0.2,
                     alpha=0.9,
                    )

    ax.clabel(ECS,
              ECS.levels,
              inline=True,
              fontsize=fs2-1,
             )

    # elevation curves
    ECS = ax.contour(lonax,latax,elevation,
                     levels=np.arange(-2,2.4,0.4),
                     colors="k",
                     linewidths=0.2,
                     linestyles="dashed",
                     alpha=0.3,
                     transform=mercator,
                    )

    # coast line
    ax.contour(lonax,latax,elevation,
               levels=[0],
               colors="k",
               linewidths=0.2,
               transform=mercator,
              )

    # cities etc
    ax.text(36.76, 25.34, "\it{Al Wajh}\n\it{platform}",
            ha="left", fontsize=fs2*0.95, zorder=20000,color="grey")

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

    # hypocenters
    if eqfile:
        idx = np.where((eqdep>z-0.5) & (eqdep<z+0.5))
        ax.scatter(eqlon[idx],eqlat[idx],marker="*",c="m",s=0.15,alpha=0.8)


    # plot vertical slices profiles
    for zslc in vslices:
        code = zslc[1]
        lon1 = zslc[-1][0]
        lat1 = zslc[-1][1]
        lon2 = zslc[-1][2]
        lat2 = zslc[-1][3]

        if code == "E":
            xoff = 0.05
            yoff = -0.06
        elif code == "C":
            xoff = 0.03
            yoff = 0.0
        else:
            xoff = 0.0
            yoff = 0.0

        ax.plot([lon1,lon2],[lat1,lat2],transform=mercator,
                c="k",linewidth=0.4,alpha=1.0,linestyle="--")

        ax.text(lon1+xoff,lat1+yoff,code,fontsize=fs2)
        ax.text(lon2-xoff,lat2+yoff,code+"'",fontsize=fs2)

    ax.text(36.1,25.8,f"{z} km depth",fontsize=fs)

    # north arrow
    arrowx, arrowy, arrowlen = 0.08, 0.25, 0.08
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
    ax.set_xticks(lonticks, crs=mercator)
    ax.xaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°",sep=""))
    ax.xaxis.set_minor_locator(pltticker.MultipleLocator(base=0.1))

    ax.set_yticks(latticks, crs=mercator)
    ax.yaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°",sep=""))
    ax.yaxis.set_minor_locator(pltticker.MultipleLocator(base=0.1))

    ax.set_title(subtit[subcount],fontsize=fs,c="k",loc="left")

    if z != 3:
        ax.set_yticklabels([])

    # zebra frame
    ax.zebra_frame(lw=3,crs=mercator,iFlag_outer_frame_in=None)

    ax.tick_params(labelsize=fs,
                   labelbottom=True, labeltop=False,labelleft=True,labelright=False,
                   bottom=True,top=False,left=True,right=False)

    xticklabels = ax.get_xticklabels()
    xticklabels[-1] = ""
    ax.set_xticklabels(xticklabels)

fig.subplots_adjust(hspace=0.3, wspace=0.1)
fig.set_size_inches(figsize)
fig.savefig(figname, dpi=dpi, bbox_inches="tight")
