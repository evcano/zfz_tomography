import os
import numpy as np
import os.path as op
import matplotlib.pyplot as plt
from BayHunter import ModelMatrix
from BayHunter import Model


resultsdir = "./results/data"
outdir = "./"

thinning = 100
vcap = 4.2

h = np.arange(0,30,0.5)
#---------------------------------------------- Model resaving and plotting
# get models from chains with highest likelihoods
fname = op.join(resultsdir, f"c_models.npy")
models = np.load(fname)
models = models[::thinning]

fig, ax = plt.subplots()

mm = ModelMatrix.get_singlemodels(models, dep_int=h)

for mtype in ["mean",]:
    vs = mm[mtype][0]
    h  = mm[mtype][1]

    if vcap:
        vs[np.where(vs>vcap)] = vcap

    vp = vs * 1.73

    # get layered version
    vp_step, vs_step, dep_step = Model.get_stepmodel_from_h(
        h, vs, vpvs=1.73, dep=h)

    if (dep_step[1] - dep_step[0]) == 0:
        vp_step = vp_step[2:]
        vs_step = vs_step[2:]
        dep_step = dep_step[2:]

    # get thickness model
    thick = np.zeros(int(dep_step.size/2))
    vs2 = thick.copy()

    a = 0
    for i in range(0, thick.size):
        b = a + 1

        thick[i] = dep_step[b] - dep_step[a]
        vs2[i] = vs_step[b]

        a = b + 1

    vp2 = vs2 * 1.73
    rho2 = vp2 * 0.32 + 0.77

    thick = np.append(1.2, thick)
    vp2 = np.append(1.5, vp2)
    vs2 = np.append(0.0, vs2)
    rho2 = np.append(1.0, rho2)

    # plot
    ax.plot(vs, h, label=mtype)
    ax.plot(vs_step, dep_step, c="k")

    # save model
    outarr = np.array([thick,vp2,vs2,rho2]).T
    outfile = op.join(outdir,f"model_{mtype}_layered.dat")
    np.savetxt(outfile, outarr, fmt="%.5f")

    outarr = np.array([h, vp, vs]).T
    outfile = op.join(outdir,f"model_{mtype}.dat")
    np.savetxt(outfile, outarr, fmt="%.5f")

ax.invert_yaxis()
ax.legend()

plt.grid()
plt.show()
plt.close()
