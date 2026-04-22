import matplotlib.pyplot as plt
import numpy as np
import rasterio
import utm
from scipy.ndimage import gaussian_filter


bfile = "../bathymetry/gebco_2023_n27.0_s23.0_w35.0_e38.0.tif"
stafile = "./mctomo_inputfiles/sreceivers.dat"
outwater_file = "./mctomo_inputfiles/waterFile.dat"

# zabargard UTM zone
utmno = 37
utmle = "R"

# mesh info in utm coordinates (km)
# dx and dy are around 1km
xmin = 198.0
xmax = 378.0
ymin = 2660.0
ymax = 2905.0
nx = 180
ny = 245

# evcano: 
# 3 s * 2 km/s  = 6km min wav
# 12s * 4.5 km/s = 54 max wav

# according to green, the shortest topography wavelenght
# that affects surface waves is 2.5 * wav

# for wav=6km, that would be 15km
# gaussian filter scalelength should be at least 15 km

# the scalelength of a gaussian filter is the
# full-width of the filter and is roughly: sigma * sqrt(8)

# give the discretization of my mesh, 1 pixel is 1 km
# sigma = 7 km, gives scalelength of ~20 km

sigma = 7  # in pixels, controls the smoothing
ndecimals = 3  # decimals to round values 

# ---------------------------------
# get mesh points
xaxis = np.linspace(xmin,xmax,nx)
yaxis = np.linspace(ymin,ymax,ny)

# get elevation of points of interest
dataset = rasterio.open(bfile)
elevation = dataset.read(1)

elevation2 = np.zeros((ny,nx))
transformer = rasterio.transform.AffineTransformer(dataset.transform)

# same as MCTOMO, inner axis is the fastest
# k is depth, j is ycoor and i is xcoor
for i in range(nx):
    for j in range(ny):
        tmp = utm.to_latlon(xaxis[i]*1000.0,yaxis[j]*1000.0,utmno,utmle)
        lat = tmp[0]
        lon = tmp[1]

        row, col= rasterio.transform.rowcol(dataset.transform, lon, lat)
        elevation2[j,i] = elevation[row, col] / 1000.0  # to km

# smooth the elevation
elevation2 = gaussian_filter(elevation2,sigma=sigma,order=0,mode='reflect')

# set topography to 0
elevation2[elevation2>=0.0] = 0.0

# set bathymetry to positve values
elevation2[elevation2<0.0] = np.abs(elevation2[elevation2<0.0])

# round elevation
elevation2 = np.round(elevation2, decimals=ndecimals)

# transform to vector
x = []
y = []
z = []

k = 0

for i in range(nx):
    for j in range(ny):
        k += 1
        x.append(xaxis[i])
        y.append(yaxis[j])
        z.append(elevation2[j,i])

x = np.array(x)
y = np.array(y)
z = np.array(z)

# save file
with open(outwater_file, "w") as _file:
    for zv in z:
        _file.write(f"{zv:.3f}\n")

# info
print("min elevation: ", np.min(z))
print("min water depth: ", np.min(z[z!=0.0]))
print("max water depth: ", np.max(z))
print("average elevation: ", np.mean(z))

# figure
fig, ax = plt.subplots()
cf = ax.scatter(x, y, c=z, s=15.0)
ax.scatter(x[z==0.0],y[z==0.0],c="brown",s=15.0)

try:
    stacoor = np.loadtxt(stafile,skiprows=1)
    ax.plot(stacoor[:,0], stacoor[:,1], ".", c="r")
except:
    pass

plt.colorbar(cf)
plt.show()
plt.close()
