import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.geoaxes
import itertools
import geopandas
import matplotlib.patches as mpatches
import matplotlib.ticker as pltticker
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import shapely
import os
from cartopy import geodesic
from glob import glob
from matplotlib.colors import LightSource,LinearSegmentedColormap,Normalize
from matplotlib.lines import Line2D
from matplotlib.patheffects import Stroke, Normal
from matplotlib import rc
from matplotlib_scalebar.scalebar import ScaleBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


rc('text', usetex=True)

def zebra_frame(self, lw=3, crs=None, zorder=None, iFlag_outer_frame_in = None):
    # Alternate black and white line segments
    bws = itertools.cycle(["k", "w"])
    self.spines["geo"].set_visible(False)

    if iFlag_outer_frame_in is not None:
        #get the map spatial reference        
        left, right, bottom, top = self.get_extent()
        crs_map = self.projection
        xticks = np.arange(left, right+(right-left)/9, (right-left)/8)
        yticks = np.arange(bottom, top+(top-bottom)/9, (top-bottom)/8)
        #check spatial reference are the same           
        pass
    else:
        crs_map =  crs
        xticks = sorted([*self.get_xticks()])
        xticks = np.unique(np.array(xticks))
        yticks = sorted([*self.get_yticks()])
        yticks = np.unique(np.array(yticks))

    for ticks, which in zip([xticks, yticks], ["lon", "lat"]):
        for idx, (start, end) in enumerate(zip(ticks, ticks[1:])):
            bw = next(bws)
            if which == "lon":
                xs = [[start, end], [start, end]]
                ys = [[yticks[0], yticks[0]], [yticks[-1], yticks[-1]]]
            else:
                xs = [[xticks[0], xticks[0]], [xticks[-1], xticks[-1]]]
                ys = [[start, end], [start, end]]

            # For first and last lines, used the "projecting" effect
            capstyle = "butt" if idx not in (0, len(ticks) - 2) else "projecting"
            for (xx, yy) in zip(xs, ys):
                self.plot(xx,
                          yy,
                          color=bw,
                          linewidth=max(0, lw - self.spines["geo"].get_linewidth()*2),
                          clip_on=False,
                          transform=crs_map,
                          zorder=zorder,
                          solid_capstyle=capstyle,
                          # Add a black border to accentuate white segments
                          path_effects=[
                              Stroke(linewidth=lw, foreground="black"),
                              Normal(),
                          ],
                         )


def plot_borders(ax):
    borders = cfeature.NaturalEarthFeature(
        category="cultural",
        name="admin_0_boundary_lines_land",
        scale="10m",
        facecolor="none")
    ax.add_feature(borders, edgecolor="black", lw=0.5, zorder=4)


def plot_study_area_limits(ax, extent_map, mercator):
    lonmin = extent_map[0]
    lonmax = extent_map[1]
    latmin = extent_map[2]
    latmax = extent_map[3]
    ax.plot([lonmin,lonmin],[latmin,latmax],linewidth=1.2,color="g",
            transform=mercator)
    ax.plot([lonmin,lonmax],[latmax,latmax],linewidth=1.2,color="g",
            transform=mercator)
    ax.plot([lonmax,lonmax],[latmax,latmin],linewidth=1.2,color="g",
            transform=mercator)
    ax.plot([lonmin,lonmax],[latmin,latmin],linewidth=1.2,color="g",
            transform=mercator)


def get_elevation_map(latax, lonax, extent_map, dataset, elevation):
    def _topo_colormap():
        N = 200
        pts = [0.0, 0.25, 0.50, 0.501, 0.75, 1.0]
        col = ["darkblue", "lightblue", "white",
               "beige", "peru", "sienna"]
        cmap = LinearSegmentedColormap.from_list('topocmap',list(zip(pts,col)),N=N)
        return cmap

    nx = lonax.size
    ny = latax.size
    elevation2 = np.zeros((ny,nx))

    for i in range(nx):
        for j in range(ny):
            lat = latax[j]
            lon = lonax[i]
            row, col= rasterio.transform.rowcol(dataset.transform, lon, lat)
            elevation2[j,i] = elevation[row, col] / 1000.0  # to km

    cmap = _topo_colormap()
    vmax = np.abs(elevation2).max()
    cmap_norm = Normalize(vmin=-vmax,vmax=vmax)

    ls = LightSource(azdeg=315.0, altdeg=45)
    shaded = ls.shade(elevation2, cmap=cmap, norm=cmap_norm,
                      blend_mode="overlay", vert_exag=0.01)

    return elevation2, shaded, cmap, cmap_norm


if __name__ == "__main__":
    stafile = "./all_stations_coordinates.txt"
    elvfile = "../bathymetry/gebco_2024_n32.0_s12.0_w30.0_e46.0.tif"
    rsrift_file = None #"../other_geophysical_data/RedSea_Dealunay/mmc1/axis-line.shp"
    figfile = "./study_region.png"

    # number of points to sample elevation
    nx = 384
    ny = 480

    # map limits
    latmin = 23.5
    latmax = 26.5
    lonmin = 34.5
    lonmax = 38.0

    # figure settings
    fs = 9
    fs2 = 6
    figsize = (6.0, 4.5)
    dpi = 400

    xticks = [34.5,35,36,37,38]
    yticks = [23.5,24,25,26,26.5]
    cbar_ticks = [-2, 0, 2]
    majloc = 1.0
    minloc = 0.1

    # inset figure settings
    latmin2 = 12.0
    latmax2 = 31.99
    lonmin2 = 30.0
    lonmax2 = 45.99

    yticks2 = [12, 21, 31]
    xticks2 = [30, 37, 45]

    nx2 = 200
    ny2 = 400

    # ---------------------------------
    # load geological features
    if rsrift_file:
        rsrift = geopandas.read_file(rsrift_file)

        # for Dealunay rift, remove rift inside ZFZ
        idx = []
        for i,tmp in enumerate(rsrift['geometry']):
            if (tmp.bounds[1] > 24.0) & (tmp.bounds[1] < 25.5):
                idx.append(i)
        rsrift['geometry'] = rsrift['geometry'].drop(index=idx)

    # ---------------------------------
    # load stations
    stanet = np.loadtxt(stafile, skiprows=1, usecols=0, dtype="U20")
    staname = np.loadtxt(stafile, skiprows=1, usecols=1, dtype="U20")
    stalat = np.loadtxt(stafile, skiprows=1, usecols=2)
    stalon = np.loadtxt(stafile, skiprows=1, usecols=3)

    # load elevation
    dataset = rasterio.open(elvfile)
    gebco_elevation = dataset.read(1)

    # MAP OF ZFZ ------------------------------------------------
    # set map limits
    extent_map = [lonmin, lonmax, latmin, latmax]
    latax = np.linspace(latmin,latmax,ny)
    lonax = np.linspace(lonmin,lonmax,nx)

    elevation, shaded, cmap, cmap_norm = get_elevation_map(
        latax, lonax, extent_map, dataset, gebco_elevation)

    # MAP OF RED SEA -------------------------------------------
    extent_map2 = [lonmin2, lonmax2, latmin2, latmax2]
    latax2 = np.linspace(latmin2, latmax2, ny2)
    lonax2 = np.linspace(lonmin2, lonmax2, nx2)

    elevation2, shaded2, cmap2, cmap_norm2 = get_elevation_map(
        latax2, lonax2, extent_map2, dataset, gebco_elevation)

    # FIGURE ---------------------------------------------------
    setattr(cartopy.mpl.geoaxes.GeoAxes, 'zebra_frame', zebra_frame)

    # set figure
    mercator = ccrs.PlateCarree()
    ax = plt.axes(projection=mercator)
    ax.set_extent(extent_map)

    # plot elevation
    im = ax.imshow(shaded, cmap, cmap_norm,
                   extent=extent_map,
                   origin="lower",
                   transform=mercator)

    # plot coastline
    ax.contour(elevation, levels=[0], colors="k",
               linewidths=0.4, extent=extent_map)

    # plot geological features
    if rsrift_file:
        ax.add_geometries(rsrift["geometry"],crs=mercator,
                           facecolor=(1,1,1,0),edgecolor="red",
                           linewidth=0.8)

    # plot stations
    sta_elevations = []
    for i in range(0,len(stanet)):
        if stanet[i] == "ZF":
            c = "white"
        elif stanet[i] == "SA":
            c = "green"
        elif stanet[i] == "KL":
            c = "blue"
            continue

        ax.scatter(stalon[i], stalat[i], s=14.0,
                   marker="s", c=c, edgecolor="k", linewidth=0.4, zorder=10000)

        ax.text(stalon[i]+0.05, stalat[i]-0.02,f"{staname[i]}",
                ha="left", fontsize=fs2,zorder=10000)

        jj = np.argmin(np.abs(latax2-stalat[i]))
        ii = np.argmin(np.abs(lonax2-stalon[i]))
        sta_elevations.append(elevation2[jj,ii])

    # cities etc
    ax.scatter(36.4689, 26.2366, s=14.0, marker="o", c="y",
               edgecolor="k", linewidth=0.4)
    ax.text(36.5, 26.27, "Al Wajh",
            ha="left", fontsize=fs2,zorder=20000)

    ax.scatter(37.2651, 25.05, s=14.0, marker="o", c="y",
               edgecolor="k", linewidth=0.4)
    ax.text(37.269, 25.08, "Umluj",
            ha="left", fontsize=fs2,zorder=20000)

    ax.scatter(36.1958, 23.6097, s=14.0, marker="o", c="y",
               edgecolor="k", linewidth=0.4)
    ax.text(36.25, 23.6097, "Zabargad Island",
            ha="left", fontsize=fs2,zorder=20000)

    ax.text(36.76, 25.34, "\it{Al Wajh}\n\it{platform}",
            ha="left", fontsize=fs2*0.95, zorder=20000,color="grey")

    # Mabahis Mons
    ax.scatter(36.03, 25.47, s=12.0, marker="^", c="k")
    ax.text(35.92, 25.46, "Mabahis\nMons",
            ha="center", fontsize=fs2,zorder=20000)

    # Mabahis Deep
    ax.scatter(36.11, 25.39, s=12.0, marker="o", c="k")
    ax.text(35.95, 25.27, "Mabahis\nDeep",
            ha="center", fontsize=fs2,zorder=20000)

    # Kebrit Deep
    ax.scatter(36.26, 24.715, s=12.0, marker="o", c="k")
    ax.text(36.1, 24.7, "Kebrit\nDeep",
            ha="center", fontsize=fs2,zorder=20000)

    # north arrow
    arrowx, arrowy, arrowlen = 0.05, 0.2, 0.07
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

    # colorbar
    cbaxes = inset_axes(ax,
                        width="20%",
                        height="2%",
                        loc="lower left",
                        bbox_to_anchor=(0.76, 0.05, 1, 1),
                        bbox_transform=ax.transAxes,
                       )

    fig = plt.gcf()
    cbar = fig.colorbar(im, cax=cbaxes, orientation="horizontal")
    cbar.ax.set_title("Elevation [km]", loc="center", fontsize=fs)
    cbar.ax.set_xticks(cbar_ticks)
    cbar.ax.tick_params(labelsize=fs)

    # legend
    lhandles = []

    lhandles.append(Line2D([0],[0],color="r",linewidth=0.8,
                           label="Rift axis"))

    lhandles.append(Line2D([0],[0],marker="o",markersize=4,
                           markerfacecolor="y",markeredgecolor="k",
                           linestyle="",markeredgewidth=0.4,
                           label=f"City",)
                   )

    lhandles.append(Line2D([0],[0],marker="o",markersize=4,
                           markerfacecolor="k",markeredgecolor="k",
                           linestyle="",
                           label=f"Submarine deep",)
                   )

    lhandles.append(Line2D([0],[0],marker="^",markersize=4,
                           markerfacecolor="k",markeredgecolor="k",
                           linestyle="",
                           label=f"Submarine volcano",)
                   )

    mkcol = {"KAUST":"b","SGS":"g","ZF":"w"}
    for llabel in ["SGS","ZF"]:
        lhandles.append(Line2D([0],[0],marker="s",markersize=4,
                               markerfacecolor=mkcol[llabel],markeredgecolor="k",
                               linestyle="",markeredgewidth=0.4,
                               label=f"{llabel} network station",)
                       )

    leg = ax.legend(handles=lhandles,
                    loc="upper right",
                    fontsize=fs2,
                    fancybox=False,
                    edgecolor="k",
                    framealpha=1.0,
                   )

    # ticks and labels
    ax.set_xticks(xticks,crs=mercator)
    ax.xaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°",sep=""))
    ax.xaxis.set_minor_locator(pltticker.MultipleLocator(base=0.1))
    ax.tick_params(axis="x", labelsize=fs,
                   labelbottom=True, labeltop=False,
                   bottom=True, top=False)

    ax.set_yticks(yticks, crs=mercator)
    ax.yaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°",sep=""))
    ax.yaxis.set_minor_locator(pltticker.MultipleLocator(base=0.1))
    ax.tick_params(axis="y", labelsize=fs,
                   labelleft=True, labelright=False,
                   left=True, right=False)

    # zebra frame
    ax.zebra_frame(lw=3, crs=mercator, iFlag_outer_frame_in=None)

    # INSET MAP OF READ SEA -----------------------------------------
    ax2 = inset_axes(ax,
                     width="45%",
                     height="45%",
                     loc="lower left",
                     bbox_to_anchor=(-0.06, 0.52, 1, 1),
                     bbox_transform=ax.transAxes,
                     axes_class=cartopy.mpl.geoaxes.GeoAxes,
                     axes_kwargs=dict(projection=mercator),
                    )

    ax2.set_extent(extent_map2)

    # plot elevation
    im2 = ax2.imshow(shaded2, cmap2, cmap_norm2,
                     extent=extent_map2,
                     origin="lower",
                     transform=mercator)

    # plot coastline
    ax2.contour(elevation2, levels=[0], colors="k",
               linewidths=0.2, extent=extent_map2)

    # ticks and labels
    ax2.set_xticks(xticks2, crs=mercator)
    ax2.xaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°",sep=""))
    ax2.tick_params(axis="x", labelsize=fs2,
                    labelbottom=True, labeltop=False,
                    bottom=True, top=False,
                    length=1.5, pad=0.5)

    ax2.set_yticks(yticks2, crs=mercator)
    ax2.yaxis.set_major_formatter(pltticker.EngFormatter(unit=u"°",sep=""))
    ax2.tick_params(axis="y", labelsize=fs2,
                    labelleft=False, labelright=True,
                    left=False, right=True,
                    length=1.5, pad=0.5)

    # plot borders
    plot_borders(ax2)

    # plot countries
    ax2.text(38.0,27.0,"Saudi Arabia",fontsize=fs2)
    ax2.text(30.5,26.0,"Egypt",fontsize=fs2)
    ax2.text(31.0,18.0,"Sudan",fontsize=fs2)

    # plot geological features
    if rsrift_file:
        ax2.add_geometries(rsrift["geometry"],crs=mercator,
                           facecolor=(1,1,1,0),edgecolor="red",
                           linewidth=0.4)

    # plot study region box
    plot_study_area_limits(ax2, extent_map, mercator)

    # INSET MAP OF READ SEA -----------------------------------------

    # save figure
    fig.set_size_inches(figsize)
    fig.tight_layout()
    if figfile:
        fig.savefig(figfile, dpi=dpi)
    else:
        plt.show()

    plt.close()
