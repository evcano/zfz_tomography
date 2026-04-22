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
from obspy.geodetics.base import gps2dist_azimuth


rc('text', usetex=True)


def plot_rays(ax,stanet,staname,stalon,stalat,staelv,pairs,mercator):
    stations = [f"{stanet[i]}.{staname[i]}" for i in range(0,stanet.size)]
    stations = np.array(stations)

    for p in pairs:
        s1, s2 = p.split("_")
        if s1 == s2:
            continue

        i1 = np.argwhere(stations == s1)[0][0]
        i2 = np.argwhere(stations == s2)[0][0]

        s1lon = stalon[i1]
        s1lat = stalat[i1]
        s2lon = stalon[i2]
        s2lat = stalat[i2]
        s1elv = staelv[i1]
        s2elv = staelv[i2]

        ax.plot([s1lon,s2lon],[s1lat,s2lat],
                color="k",alpha=0.8,linewidth=0.4,
                transform=mercator,
               )
        dist, _, _ = gps2dist_azimuth(s1lat,s1lon,s2lat,s2lon)
        elvdif = np.abs(s1elv - s2elv)
        plt.plot(dist, elvdif*1000, "x")
        ds = np.sqrt(dist**2 + (elvdif*1000)**2)
        print(dist/ds, p, dist, elvdif*1000)


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
    ax.plot([lonmin,lonmin],[latmin,latmax],linewidth=1.2,color="w",
            transform=mercator)
    ax.plot([lonmin,lonmax],[latmax,latmax],linewidth=1.2,color="w",
            transform=mercator)
    ax.plot([lonmax,lonmax],[latmax,latmin],linewidth=1.2,color="w",
            transform=mercator)
    ax.plot([lonmin,lonmax],[latmin,latmin],linewidth=1.2,color="w",
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
    stafile = "./stations_coordinates.txt"
    elvfile = "../bathymetry/gebco_2024_n32.0_s12.0_w30.0_e46.0.tif"
    obsdir = "../3d_inversion/observed_times"
    figfile = "./raycoverage.png"

    # number of points to sample elevation
    nx = 384 
    ny = 480

    # map limits
    latmin = 24.
    latmax = 26.2
    lonmin = 36.
    lonmax = 37.8

    # figure settings
    fs = 9
    fs2 = 6
    figsize = (7.0, 4.5)
    dpi = 400

    yticks = [24,25,26,26.2]
    xticks = [36,37,37.8]

    cbar_ticks = [-2, 0, 2]
    majloc = 1.0
    minloc = 0.1

    # ---------------------------------

    # load stations
    stanet = np.loadtxt(stafile, skiprows=1, usecols=0, dtype="U20")
    staname = np.loadtxt(stafile, skiprows=1, usecols=1, dtype="U20")
    stalat = np.loadtxt(stafile, skiprows=1, usecols=2, dtype=float)
    stalon = np.loadtxt(stafile, skiprows=1, usecols=3, dtype=float)

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

    # FIGURE ---------------------------------------------------
    setattr(cartopy.mpl.geoaxes.GeoAxes, 'zebra_frame', zebra_frame)

    for fi, cmp_ in enumerate(["ZZ","TT"]):
        # set figure
        mercator = ccrs.PlateCarree()
        ax = plt.subplot(int(f"12{fi+1}"), projection=mercator)
        ax.set_extent(extent_map)

        # plot elevation
        im = ax.imshow(shaded, cmap, cmap_norm,
                       extent=extent_map,
                       origin="lower",
                       transform=mercator)

        # plot coastline
        ax.contour(elevation, levels=[0], colors="k",
                   linewidths=0.4, extent=extent_map)

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

            if staname[i] == "LAVA":
                looff = -0.1
            else:
                looff = 0.05

            jj = np.argmin(np.abs(latax-stalat[i]))
            ii = np.argmin(np.abs(lonax-stalon[i]))
            sta_elevations.append(elevation[jj,ii])

        sta_elevations = np.array(sta_elevations)

        # plot rays
        pairs = os.listdir(os.path.join(obsdir,cmp_))
        pairs = [x.replace(".dat","") for x in pairs]
        pairs = [x.replace("_TT","") for x in pairs]
        pairs = [x.replace("_ZZ","") for x in pairs]
        pairs = [x.replace("CDisp.T.","") for x in pairs]
        pairs = [x.replace("GDisp.","") for x in pairs]
        pairs = list(set(pairs))
        pairs.sort()

        plot_rays(ax,stanet,staname,stalon,stalat,sta_elevations,pairs,mercator)
        print("rays plotted")

        # cities etc

        ax.scatter(37.2651, 25.05, s=14.0, marker="o", c="y",
                   edgecolor="k", linewidth=0.4)
        ax.text(37.269, 25.08, "Umluj",
                ha="left", fontsize=fs2,zorder=20000)

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


        # north arrow
        arrowx, arrowy, arrowlen = 0.08, 0.23, 0.08
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
                            bbox_to_anchor=(0.72, 0.05, 1, 1),
                            bbox_transform=ax.transAxes,
                           )

        fig = plt.gcf()
        cbar = fig.colorbar(im, cax=cbaxes, orientation="horizontal")
        cbar.ax.set_title("Elevation [km]", loc="center", fontsize=fs)
        cbar.ax.set_xticks(cbar_ticks)
        cbar.ax.tick_params(labelsize=fs)

        # ticks and labels
        if fi == 0:
            llab = True
        else:
            llab = False

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
                       labelleft=llab, labelright=False,
                       left=llab, right=False)

        # zebra frame
        ax.zebra_frame(lw=3, crs=mercator, iFlag_outer_frame_in=None)

        if fi == 0:
            title_ = r"\textbf{(a)} Ray coverage of Rayleigh waves"
        else:
            title_ = r"\textbf{(b)} Ray coverage of Love waves"

        ax.set_title(title_,loc="left",fontsize=fs)

    # save figure
    fig.set_size_inches(figsize)
    fig.tight_layout()

    if figfile:
        fig.savefig(figfile, dpi=dpi)
    else:
        plt.show()

    plt.close()
