import matplotlib.pyplot as plt
import numpy as np
import utm
import os
from scipy.interpolate import interp1d


def add_water_layer(xaxis,yaxis,zaxis,vs,waterDepth):
    dz = zaxis[-1] - zaxis[-2]

    new_zmax = zaxis.max() + waterDepth.max()
    new_zaxis = np.arange(0.0, new_zmax + dz, dz)

    vs2 = np.zeros((yaxis.size, xaxis.size, new_zaxis.size))

    for i in range(0,xaxis.size):
        for j in range(0,yaxis.size):
            new_z = np.append(0.0, zaxis + waterDepth[j,i])
            new_vs = np.append(vs[j,i,0], vs[j,i,:])

            iobj = interp1d(new_z, new_vs, bounds_error=False, fill_value=new_vs[-1])
            new_vs = iobj(new_zaxis)

            new_vs[np.where(new_zaxis < waterDepth[j,i])] = 0.0
            vs2[j,i,:] = new_vs

    return vs2, new_zaxis


def get_mean_initial_model(mpath,nchains):
    for i in range(1,nchains+1):
        fname = os.path.join(mpath,f"model_{i:03}.dat")
        m = np.loadtxt(fname,usecols=(0,2),skiprows=1)
        if i == 1:
            mean = m.copy()
        else:
            mean[:,1] += m[:,1]

    mean[:,1] /= nchains
    return mean


def get_perturbations(model, zaxis, inimodel):
    iobj = interp1d(inimodel[:,0],inimodel[:,1])
    inimodel2 = iobj(zaxis)

    model2 = model.copy()
    for i in range(0,model2.shape[0]):
        for j in range(0,model2.shape[1]):
            model2[i,j,:] = np.log(np.divide(model[i,j,:],inimodel2)) * 100

    return model2


# Unweighted_average.dat contains both VP and VS models, first VP then VS
# The lines of the file correspond to a (x,y) point, columns to a (z) point
# The lines of the file loop over x then over y
# E.G.,
# first line corresponds to (x=1,y=1) and contains (z) values from z=1 to z=zmax
# second line corresponds to (x=1,y=2) and contains (z) values from z=1 to z=max
# and so on ...

mpath0 = "../3d_inversion/final_chains_2million/all_chains/Results"
waterfile = "../3d_inversion/mctomo_inputfiles/waterFile.dat"
inimodels = "../3d_inversion/initial_models"

# mesh limits
xmin = 198.0
xmax = 378.0
ymin = 2660.0
ymax = 2905.0
zmin = 0.0
zmax = 20.0

utmno = 37
utmle = "R"

# ------------------------------------------------
ini_mod = get_mean_initial_model(inimodels,20)

for mtype in ["std"]:
    # read model
    mpath = os.path.join(mpath0, f"Unweighted_{mtype}.dat")

    nz, ny, nx2 = np.loadtxt(mpath, max_rows=1, dtype=int)
    MODEL = np.loadtxt(mpath, skiprows=1)  # 2*nx*ny rows and nz columns

    # get vs model
    nx = int(nx2 / 2)
    nxny = int(nx*ny)
    vp = MODEL[0:nxny, :]
    vs = MODEL[nxny:, :]

    # reshape model array
    vs = vs.flatten(order="C")  # row-wise flattening
    vs = np.reshape(vs, (nx,ny,nz), order="C")  # convert model to 3D array
    vs = np.transpose(vs, (1, 0, 2))  # change axes to y,x,z

    # ------------------------------------------------
    # create coordinates axes
    xaxis = np.linspace(xmin, xmax, nx)
    yaxis = np.linspace(ymin, ymax, ny)
    zaxis = np.linspace(zmin, zmax, nz)

    print(xaxis[-1]-xaxis[-2])
    print(yaxis[-1]-yaxis[-2])
    print(zaxis[-1]-zaxis[-2])
    print(zaxis[0],zaxis[-1])
    print(nx,ny,nz)

    # get perturbations
    # vs_per = get_perturbations(vs,zaxis,ini_mod)

    # ------------------------------------------------
    # read water layer depth
    waterDepth = np.loadtxt(waterfile, dtype=float)
    waterDepth = np.reshape(waterDepth, (nx, ny), order="C")
    waterDepth = waterDepth.T  # change axes to y,x
    print("water min max", waterDepth.min(), waterDepth.max())

    # add water layer to model
    zaxis_bak = zaxis.copy()
    vs, zaxis = add_water_layer(xaxis,yaxis,zaxis,vs,waterDepth)
    # vs_per, _ = add_water_layer(xaxis,yaxis,zaxis_bak,vs_per,waterDepth)

    # ------------------------------------------------
    # geto geographical coordinates
    print(nx,ny,nz)
    nz = zaxis.size

    LATMESH = np.zeros((ny,nx,nz))
    LONMESH = np.zeros((ny,nx,nz))
    DEPMESH = np.zeros((ny,nx,nz))

    for i in range(nx):
        for j in range(ny):
            tmp = utm.to_latlon(xaxis[i]*1000.0,yaxis[j]*1000.0,utmno,utmle)
            LATMESH[j,i,:] = tmp[0]
            LONMESH[j,i,:] = tmp[1]
            DEPMESH[j,i,:] = zaxis
    print(nx,ny,nz)
    # ------------------------------------------------
    with open ("std_model_ascii.txt","w") as _file:
        for k in range(nz):
            for i in range(nx):
                for j in range(ny):
                    x = LONMESH[j,i,k]
                    y = LATMESH[j,i,k]
                    z = DEPMESH[j,i,k]
                    v = vs[j,i,k]

                    if z > 12.0:
                        continue

                    txt = f"{x:.6f} {y:.6f} {z:.6f} {v:.2f}\n"
                    _file.write(txt)
