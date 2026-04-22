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


