import matplotlib.pyplot as plt
import numpy as np
import rasterio
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


def array_polygon(sta_list, sta_lat, sta_lon, stabound):
    points = []
    for sta in stabound:
        idx = np.argwhere(sta_list == sta).flatten()
        points.append((sta_lon[idx], sta_lat[idx]))
    polygon = Polygon(points)
    return polygon


bfile = "./gebco_2023_n27.0_s23.0_w35.0_e38.0.tif"

#all stations
#stafile = "../3d_inversion/stations_coordinates.txt"
#
#stabound = ["EWJHS", "QUMAN", "KHUF", "OBS04",
#            "LAVA", "OBS01", "OBS02", "OBS11", "NORTH", "OBS12", "EWJHS"]

# obs only
stafile = "../3d_inversion/obs_coordinates.txt"
stabound = ["OBS12", "OBS09", "OBS03", "OBS01", "OBS02", "OBS11", "NORTH", "OBS12"]

nlon = 100
nlat = 100

# ---------------------------------
# read stations coordinates
net_list = np.loadtxt(stafile, usecols=(0), dtype="U20")
sta_list = np.loadtxt(stafile, usecols=(1), dtype="U20")
sta_lat= np.loadtxt(stafile, usecols=(2))
sta_lon= np.loadtxt(stafile, usecols=(3))

lonmin = sta_lon.min()
lonmax = sta_lon.max()
latmin = sta_lat.min()
latmax = sta_lat.max()

# get polygon surrounding the stations
polygon = array_polygon(sta_list, sta_lat, sta_lon, stabound)

# read bathymetry
dataset = rasterio.open(bfile)
topography = dataset.read(1)

# get bathymetry of points of interest
lon_ax = np.linspace(lonmin, lonmax, nlon)
lat_ax = np.linspace(latmin, latmax, nlat)
topography2 = np.zeros((nlat, nlon))

transformer = rasterio.transform.AffineTransformer(dataset.transform)

for i, lon in enumerate(lon_ax):
    for j, lat in enumerate(lat_ax):
        if polygon.contains(Point(lon, lat)):
            row, col= rasterio.transform.rowcol(dataset.transform, lon, lat)
            topography2[j,i] = topography[row, col]
        else:
            topography2[j,i] = np.nan

# info
print("min elevation: ", np.nanmin(topography2[:]))
print("max elevation: ", np.nanmax(topography2[:]))
print("average elevation: ", np.nanmean(topography2[:]))

# figure
fig, ax = plt.subplots()

cf = ax.contourf(lon_ax, lat_ax, topography2, levels=100,
                 cmap="gist_earth", vmin=-1800, vmax=0)

ax.plot(sta_lon, sta_lat, ".", c="r")
for i in range(sta_list.size):
    ax.text(sta_lon[i], sta_lat[i], sta_list[i])

x, y = polygon.exterior.xy
ax.plot(x, y)

plt.colorbar(cf)
plt.show()
plt.close()
