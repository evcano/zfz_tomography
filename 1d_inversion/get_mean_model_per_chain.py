import os
import numpy as np
import os.path as op
import matplotlib.pyplot as plt
from BayHunter import ModelMatrix
from BayHunter import Model
from scipy.signal import convolve, windows


resultsdir = "./results/data"
outdir = "./mean_models"

nchains = 100
nmodels = 20  # number of requested models
outlier_chains = [38,48,60]

mtype = "mean"
vcap = 4.2

thinning = 100  # i already thined individual chains by 2
depint = np.arange(0,21,0.2)

#---------------------------------------------- Model resaving and plotting
# determine chains with highest likelihood
medianlikes = []

for i in range(nchains):
    fname = op.join(resultsdir, f"c{i:03}_p2likes.npy")

    if not os.path.isfile(fname):
        medianlikes.append(0.0)
    else:
        likes = np.load(fname)
        medianlikes.append(np.median(likes))

idx = np.argsort(medianlikes)  # sort from lowest to highets
idx = idx[::-1]  # from highest to lowest
chains_idx = idx[0:nmodels+1]

fig, ax = plt.subplots()

# get models from chains with highest likelihoods
for j, i in enumerate(chains_idx):
    if i in outlier_chains:
        continue

    # read models
    fname = op.join(resultsdir, f"c{i:03}_p2models.npy")
    models = np.load(fname)
    models = models[::thinning]

    # get mean model
    mm = ModelMatrix.get_singlemodels(models, dep_int=depint)
    vs = mm[mtype][0]
    dep = mm[mtype][1]

    # cap velocity
    vs[np.where(vs>vcap)] = vcap

    # smooth
    n = 13
    n2 = int((n-1)/2)
    win = windows.hann(n)
    vs = np.pad(vs, (n2,n2), mode="edge")
    vs = convolve(vs, win, "valid") / np.sum(win)

    # get vp
    vp = vs * 1.73

    # make sure depth starts at 0
    dep[0] = 0.0

    # plot
    ax.plot(vs, dep)

    # append headers required by MCTOMO
    dep = np.append(dep.size, dep)
    vp = np.append(3, vp)
    vs = np.append(0, vs)

    # save model
    outarr = np.array([dep, vp, vs]).T
    outfile = op.join(outdir,f"model_{j:03}.dat")
    np.savetxt(outfile, outarr, fmt="%.5f")


ax.invert_yaxis()
plt.grid()
plt.show()
plt.close()
