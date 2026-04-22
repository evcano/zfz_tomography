import os
import numpy as np
import os.path as op
import matplotlib.pyplot as plt
from BayHunter import ModelMatrix
from BayHunter import Model


resultsdir = \
"/data/valeroe/red_sea_obs_2/1d_inversion/obs_inversion/results/data"

outdir = "./models_per_chain"

nchains = 100
thinning = 100  # keep in mind I already thinned every 2 models

vpvs = 1.73
depint = np.arange(0,21,0.2)

#----------------------------------------------
for i in range(0, nchains):
    if i == 60:
        continue

    # load models
    fname = op.join(resultsdir, f"c{i:03}_p2models.npy")
    models = np.load(fname)

    # initial model
    initial_model = models[0].copy()
    vpi, vsi, hi = Model.get_vp_vs_h(initial_model, vpvs)
    vpi, vsi, hi = Model.get_stepmodel_from_h(vp=vpi, vs=vsi, h=hi)

    # last model
    last_model = models[-1].copy()
    vpf, vsf, hf = Model.get_vp_vs_h(last_model, vpvs)
    vpf, vsf, hf = Model.get_stepmodel_from_h(vp=vpf, vs=vsf, h=hf)

    # mtype model
    models = models[::thinning]
    for mtype in ["mean","median","mode"]:
        mm = ModelMatrix.get_singlemodels(models, dep_int=depint)
        vs = mm[mtype][0]
        dep = mm[mtype][1]

        # save model
        outarr = np.array([dep, vs]).T
        outfile = op.join(outdir, f"model_{i:03}_{mtype}.dat")
        np.savetxt(outfile, outarr, fmt="%.5f")

    # save initial model
    outarr = np.array([hi, vsi]).T
    outfile = op.join(outdir, f"model_{i:03}_initial.dat")
    np.savetxt(outfile, outarr, fmt="%.5f")

    # save final model
    outarr = np.array([hf, vsf]).T
    outfile = op.join(outdir, f"model_{i:03}_final.dat")
    np.savetxt(outfile, outarr, fmt="%.5f")

    print(i," done")
