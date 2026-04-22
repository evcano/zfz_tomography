import matplotlib.pyplot as plt
import numpy as np
import os
import utm


# NOTE: elevation of sources/recivers is not used when modelling surface waves
stafile = "./stations_coordinates.txt"
outsources_file = "./mctomo_inputfiles/sreceivers.dat"

# zabargard UTM zone
utmno = 37
utmle = "R"

extdist_x = 8.0  # extend domain by this amount (km)
extdist_y = 9.0

nx = 180
ny = 245
nz = 81

# ------------------------------------
# read stations coordinates
net_list = np.loadtxt(stafile, usecols=(0), dtype="U20")
sta_list = np.loadtxt(stafile, usecols=(1), dtype="U20")
sta_lat = np.loadtxt(stafile, usecols=(2))
sta_lon = np.loadtxt(stafile, usecols=(3))

print(np.mean(sta_lon), np.mean(sta_lat))

lonmin = sta_lon.min()
lonmax = sta_lon.max()
latmin = sta_lat.min()
latmax = sta_lat.max()

# convert station coordinates to cartesian system
sta_x = []
sta_y = []

for i in range(sta_lon.size):
    tmp = utm.from_latlon(sta_lat[i], sta_lon[i],utmno, utmle)
    sta_x.append(tmp[0])
    sta_y.append(tmp[1])

sta_x = np.array(sta_x) / 1000.0  # to km
sta_y = np.array(sta_y) / 1000.0  # to km

# save sources file
nsources = sta_list.size
ndim = 3
with open(outsources_file, "w") as _file:
    _file.write(f"{nsources} {ndim}\n")
    for i in range(nsources):
        _file.write(f"{sta_x[i]:.5f} {sta_y[i]:.5f} 0.0\n")

# info 
xmin = np.floor(sta_x.min()) - extdist_x
xmax = np.ceil(sta_x.max()) + extdist_x
xdis = xmax - xmin

ymin = np.floor(sta_y.min()) - extdist_y
ymax = np.ceil(sta_y.max()) + extdist_y
ydis = ymax - ymin

print("xmin/xmax: ", xmin, xmax)
print("ymin/ymax: ", ymin, ymax)
print("xdis/ydis: ", xdis, ydis)

xpts = np.linspace(xmin, xmax, nx)
ypts = np.linspace(ymin, ymax, ny)

print("nx: ", nx)
print("ny: ", ny)
print("dx: ", xpts[1]-xpts[0])
print("dy: ", ypts[1]-ypts[0])

# figure
fig, ax = plt.subplots()

ax.plot(sta_x, sta_y, "o")
for i in range(len(sta_list)):
    ax.text(sta_x[i], sta_y[i], sta_list[i])
for x in xpts:
    ax.axvline(x, color="gray", alpha=0.5)
for y in ypts:
    ax.axhline(y, color="gray", alpha=0.5)

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
plt.show()
