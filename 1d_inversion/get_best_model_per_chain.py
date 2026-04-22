import os
import numpy as np
import os.path as op
import matplotlib.pyplot as plt
from BayHunter import ModelMatrix
from BayHunter import Model


resultsdir = "./results/data"
outdir = "./high_like_models"

nchains = 100
nmodels = 20 # number of requested models
outlier_chains = [38,48,60]

thinning = 100
depint = np.arange(0,21,0.1)
vcap = 4.2

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

    fname = op.join(resultsdir, f"c{i:03}_p2misfits.npy")
    misfits = np.load(fname).T[-1]

    fname = op.join(resultsdir, f"c{i:03}_p2models.npy")
    models = np.load(fname)

    mm = models[np.argmin(misfits)]
    _, vs, dep = Model.get_stepmodel(mm, vpvs=1.73)

    # get vp
    vp = vs * 1.73

    # make sure depth starts at 0
    dep[0] = 0.0

    ax.plot(vs, dep)

    # append headers required by MCTOMO
    dep = np.append(dep.size, dep)
    vp = np.append(3, vp)
    vs = np.append(0, vs)

    # save model
    outarr = np.array([dep, vp, vs]).T
    outfile = op.join(outdir,f"model_{j:03}.dat")
    np.savetxt(outfile, outarr, fmt="%.5f")

plt.ylim(0,20)
ax.invert_yaxis()
plt.grid()


plt.show()
plt.close()
