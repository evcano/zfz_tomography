import matplotlib.pyplot as plt
import numpy as np
import os
import plotly.graph_objects as go
from scipy.interpolate import interp1d
from cmcrameri import cm


def add_water_layer(xaxis,yaxis,zaxis,vs,waterDepth):
    dz = zaxis[-1] - zaxis[-2]
    new_zaxis = np.arange(0.0, zaxis.max()+waterDepth.max() + dz, dz)

    vs2 = np.zeros((yaxis.size, xaxis.size, new_zaxis.size))

    for i in range(0,xaxis.size):
        for j in range(0,yaxis.size):
            new_z = np.append([0.0, waterDepth[j,i]], zaxis + waterDepth[j,i])
            new_vs = np.append([0.0, 0.0], vs[j,i,:])

            iobj = interp1d(new_z, new_vs, bounds_error=False, fill_value=new_vs[-1])
            new_vs = iobj(new_zaxis)

            vs2[j,i,:] = new_vs

    return vs2, new_zaxis


def get_the_slice(x, y, z, surface, alpha=1.0):
    z[np.where(surface==0.0)] = np.nan

    surf = go.Surface(x=x,
                      y=y,
                      z=z,
                      surfacecolor=surface,
                      connectgaps=False,
                      coloraxis="coloraxis",
                      opacity=alpha,
                     )

    return surf


def get_lims_colors(surface):
    return np.min(surface), np.max(surface)


def mpl_to_plotly(cmap, pl_entries):
    scale = np.linspace(0, 1, pl_entries)
    colors = (cmap(scale)[:, :3]*255).astype(np.uint8)
    pl_colorscale = [[s, f'rgb{tuple(color)}']
                     for s, color in zip(scale, colors)]
    return pl_colorscale

# Unweighted_average.dat contains both VP and VS models, first VP then VS
# The lines of the file correspond to a (x,y) point, columns to a (z) point
# The lines of the file loop over x then over y
# E.G.,
# first line corresponds to (x=1,y=1) and contains (z) values from z=1 to z=zmax
# second line corresponds to (x=1,y=2) and contains (z) values from z=1 to z=max
# and so on ...

waterfile = "./mctomo_inputfiles/waterFile.dat"
mpath = "/data/valeroe/red_sea_obs/3d_inversion/final_chains/all_chains/Results"
figpath = "/data/valeroe/red_sea_obs/3d_inversion/figures"
mtype = "average"

xmin = 198.0
xmax = 378.0
ymin = 2660.0
ymax = 2905.0
zmin = 0.0
zmax = 20.0

# ------------------------------------------------
# read stations
stafile = "./mctomo_inputfiles/sreceivers.dat"
stacoor = np.loadtxt(stafile,skiprows=1)
# ------------------------------------------------

# read model
mpath = os.path.join(mpath, f"Unweighted_{mtype}.dat")

nz, ny, nx2 = np.loadtxt(mpath, max_rows=1, dtype=int)
MODEL = np.loadtxt(mpath, skiprows=1)  # 2*nx*ny rows and nz columns

# get vs model
nx = int(nx2 / 2)
nxny = int(nx*ny)
vp = MODEL[0:nxny, :]
vs = MODEL[nxny:, :]

vsmin = vs.min()
vsmax = vs.max()

# flatten model to vectors
# the vectors contain model parameters in the following order
# loop over x; loop over y; loop over z
vs = vs.flatten(order="C")  # row-wise flattening
vs = np.reshape(vs, (nx,ny,nz), order="C")  # convert model to 3D array
vs = np.transpose(vs, (1, 0, 2))  # change axes to y,x,z
# ------------------------------------------------

# create coordinates axes
xaxis = np.linspace(xmin, xmax, nx)
yaxis = np.linspace(ymin, ymax, ny)
zaxis = np.linspace(zmin, zmax, nz)

# ------------------------------------------------

# read water layer depth
waterDepth = np.loadtxt(waterfile, dtype=float)
waterDepth = np.reshape(waterDepth, (nx, ny), order="C")
waterDepth = waterDepth.T  # change axes to y,x

# add water layer
vs, zaxis = add_water_layer(xaxis,yaxis,zaxis,vs,waterDepth)

# ------------------------------------------------

# FIGURES
if mtype == "std":
    cmap = cm.vik.resampled(25)
elif mtype == "average":
    cmap = cm.roma.resampled(25)

alpha = 1.0
vsmax = 5.0

pl_colorscale = mpl_to_plotly(cmap, 200)

coloraxis = {'autocolorscale': False,
             'colorscale':pl_colorscale,
             'cmin': vsmin,
             'cmax': vsmax,
            }

all_slices = []

# z slices
for z in [5]:
    idx = np.argmin(abs(zaxis-z))

    xgrd, ygrd = np.meshgrid(xaxis, yaxis)
    zgrd = z * np.ones(xgrd.shape)

    slice_z = get_the_slice(xgrd, ygrd, zgrd, vs[:,:,idx] ,alpha)

    all_slices.append(slice_z)

# x slices
for x in [250]:
    idx = np.argmin(abs(xaxis-x))

    ygrd, zgrd = np.meshgrid(yaxis, zaxis)
    xgrd = x * np.ones(ygrd.shape)

    slice_x = get_the_slice(xgrd, ygrd, zgrd, vs[:,idx,:].T, alpha)
    all_slices.append(slice_x)

# y slices
for y in [2750]:
    idx = np.argmin(abs(yaxis-y))

    xgrd, zgrd = np.meshgrid(xaxis, zaxis)
    ygrd = y * np.ones(xgrd.shape)

    slice_y = get_the_slice(xgrd, ygrd, zgrd, vs[idx,:,:].T, alpha)

    all_slices.append(slice_y)

# stations
sta3d = go.Scatter3d(x=stacoor[:,0],y=stacoor[:,1],z=stacoor[:,0]*0.0,
                     mode='markers',
                     marker=dict(size=2.0, color="magenta"),
                    )

all_slices.append(sta3d)

# 3d figure
fig1 = go.Figure(data=all_slices)

fig1.update_layout(
    title_text="slices",
    scene_xaxis_range=[200,300],
    scene_yaxis_range=[2650,2850],
    scene_zaxis_range=[10, 0.0],
    scene=dict(aspectratio=dict(x=1,y=2,z=0.5)),
    coloraxis=coloraxis,
)

fig1.write_html('/home/valeroe/myfig.html')
