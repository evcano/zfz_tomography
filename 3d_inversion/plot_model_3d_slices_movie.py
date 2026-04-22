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


def frame_args(duration):
    return {
        "frame": {"duration": duration},
        "mode": "immediate",
        "fromcurrent": True,
        "transition": {"duration": duration, "easing": "linear"},
    }


def get_updatemenus():
    updatemenus = [
        {
            "buttons": [
                {
                    "args": [None, frame_args(50)],
                    "label": "&#9654;",
                    "method": "animate",
                },
                {
                    "args": [[None], frame_args(0)],
                    "label": "&#9724;",
                    "method": "animate",
                },
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 70},
            "type": "buttons",
            "x": 0.1,
            "y": 0,
        }
    ]

    return updatemenus


def get_sliders(fig):
    sliders = [
        {
            "pad": {"b": 10, "t": 60},
            "len": 0.9,
            "x": 0.1,
            "y": 0,
            "steps": [{"args": [[f.name], frame_args(0)],
                       "label": f"{f.name} km",
                        "method": "animate",
                      }
                        for f in fig.frames
                     ],
        }
    ]

    return sliders


def get_the_slice(x, y, z, surface, alpha=1.0):
    z[np.where(surface==0.0)] = np.nan  # mask water
    customdata = surface.T  # plotly has a bug that transposes surfacecolor
    htemplate = 'x=%{x:.2f},<br>y=%{y:.2f},<br>z=%{z:.2f},<br>c=%{customdata:.2f}'

    surf = go.Surface(x=x,
                      y=y,
                      z=z,
                      surfacecolor=surface,
                      customdata=customdata,  # to show correct info in hover
                      connectgaps=False,
                      coloraxis="coloraxis",
                      opacity=alpha,
                      hovertemplate=htemplate,
                      hoverlabel=dict(bgcolor='rgba(255,255,255,0.1)'),
                      name="",
                     )

    return surf


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
mpath = "/data/valeroe/red_sea_obs/3d_inversion/final_chains_3million/all_chains/Results"
figpath = "/data/valeroe/red_sea_obs/3d_inversion/figures"

mtype = "average"
maxis = "x"

# figures limits
xmin2 = 190
xmax2 = 320
ymin2 = 2650
ymax2 = 2870
zmin2 = 0
zmax2 = 10

# mesh limits
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
stanames = np.loadtxt("./stations_coordinates.txt", usecols=(1),
                      dtype="U20")
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

# colormap
alpha = 1.0

if mtype == "std":
    cmap = cm.devon.reversed()
    cmap = cmap.resampled(21)
    vmin = 0.0
    vmax = 1.0
elif mtype == "average":
    #cmap = cm.roma.resampled(51)
    cmap = cm.roma
    vmin = 1.0
    vmax = 4.5

pl_colorscale = mpl_to_plotly(cmap, 200)

coloraxis = {'autocolorscale': False,
             'colorscale':pl_colorscale,
             'cmin': vmin,
             'cmax': vmax,
            }

# stations
sta3d = go.Scatter3d(x=stacoor[:,0],
                     y=stacoor[:,1],
                     z=stacoor[:,0]*0.0,
                     text=stanames,
                     marker=dict(size=2.0, color="magenta"),
                     mode='markers+text',
                    )

# frames
frames = []

if maxis == "z":
    for z in np.arange(zmin2, zmax2+0.5, 0.5):
        idx = np.argmin(abs(zaxis-z))

        xgrd, ygrd = np.meshgrid(xaxis, yaxis)
        zgrd = z * np.ones(xgrd.shape)

        slice_ = get_the_slice(xgrd, ygrd, zgrd, vs[:,:,idx] ,alpha)
        frames.append(go.Frame(data=slice_, name=f"{z:.1f}"))

elif maxis == "x":
    for x in np.arange(xmin2, xmax2+2, 5):
        idx = np.argmin(abs(xaxis-x))

        ygrd, zgrd = np.meshgrid(yaxis, zaxis)
        xgrd = x * np.ones(ygrd.shape)

        slice_ = get_the_slice(xgrd, ygrd, zgrd, vs[:,idx,:].T, alpha)
        frames.append(go.Frame(data=slice_, name=f"{x:.1f}"))

elif maxis == "y":
    for y in np.arange(ymin2, ymax2, 5):
        idx = np.argmin(abs(yaxis-y))

        xgrd, zgrd = np.meshgrid(xaxis, zaxis)
        ygrd = y * np.ones(xgrd.shape)

        slice_ = get_the_slice(xgrd, ygrd, zgrd, vs[idx,:,:].T, alpha)
        frames.append(go.Frame(data=slice_, name=f"{y:.1f}"))

# MOVIE
fig = go.Figure(data=[sta3d, sta3d], frames=frames)

updatemenus = get_updatemenus()

sliders = get_sliders(fig)

whr = (ymax2-ymin2) / (xmax2 - xmin2)

scene = dict(
    xaxis=dict(range=[xmin2-1, xmax2+1], autorange=False),
    yaxis=dict(range=[ymin2-1, ymax2+1], autorange=False),
    zaxis=dict(range=[zmax2+1, zmin2-1], autorange=False),
    aspectratio=dict(x=1, y=whr, z=0.5),
    camera=dict(eye=dict(x=1.25,y=-1.25,z=1.25)),
    xaxis_title="easting [km]",
    yaxis_title="northing [km]",
    zaxis_title="depth [km]",
)


fig.update_layout(
    scene=scene,
    coloraxis=coloraxis,
    sliders=sliders,
    updatemenus=updatemenus,
    title=dict(text=f"{mtype}"),
    coloraxis_colorbar=dict(title="[km/s]"),
)

fname = f"movie_{mtype}_{maxis}.html"
fig.write_html(os.path.join(figpath,fname))
