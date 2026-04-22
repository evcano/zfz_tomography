# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

import os
import numpy as np
import os.path as op
import matplotlib
import matplotlib.pyplot as plt
from BayHunter import PlotFromStorage
from BayHunter import utils


# Load priors and initparams from config.ini or simply create dictionaries.
initfile = './config.ini'
nchains = 100
maxmodels = 500000

#  ---------------------------------------------- Model resaving and plotting
priors, initparams = utils.load_params(initfile)
path = initparams['savepath']

cfile = f"{initparams['station']}_config.pkl"
configfile = op.join(path, 'data', cfile)

obj = PlotFromStorage(configfile)
obj.save_final_distribution(maxmodels=maxmodels, dev=0.05)

fig = obj.plot_bestdatafits()
fig.savefig(op.join(path, "bestdata.png"))

fig = obj.plot_bestmodels()
fig.savefig(op.join(path, "bestmodels.png"))

fig = obj.plot_posterior_models1d(depint=0.2)
fig.savefig(op.join(path, "models1d.png"))

fig = obj.plot_posterior_models2d(depint=0.2)
fig.savefig(op.join(path, "models2d.png"))


fig = obj.plot_iiterlikes(nchains=nchains)
fig.savefig(op.join(path, f"likes.png"))

obj.save_plots(nchains=4)
