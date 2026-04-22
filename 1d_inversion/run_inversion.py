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
from BayHunter import PlotFromStorage
from BayHunter import Targets
from BayHunter import utils
from BayHunter import MCMC_Optimizer
from BayHunter import ModelMatrix
from BayHunter import SynthObs
import logging


# set os.environment variables to ensure that numerical computations
# do not do multiprocessing !! Essential !! Do not change !
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

# console printout formatting
formatter = ' %(processName)-12s: %(levelname)-8s |  %(message)s'
logging.basicConfig(format=formatter, level=logging.INFO)
logger = logging.getLogger()

# ------------------------------------------------------------  LOAD CONFIG
initfile = './config.ini'
priors, initparams = utils.load_params(initfile)

#  -----------------------------------------------------------  DEFINE TARGETS

# love group
xsw, ysw = np.loadtxt('./observed_obs/tt-group-obs.txt').T
target1 = Targets.LoveDispersionGroup(xsw, ysw)

# love phase
xsw, ysw = np.loadtxt('./observed_obs/tt-phase-obs.txt').T
target2 = Targets.LoveDispersionPhase(xsw, ysw)

# rayleigh group
xsw, ysw = np.loadtxt('./observed_obs/zz-group-obs.txt').T
target3 = Targets.RayleighDispersionGroup(xsw, ysw)

# rayleigh phase
xsw, ysw = np.loadtxt('./observed_obs/zz-phase-obs.txt').T
target4 = Targets.RayleighDispersionPhase(xsw, ysw)

# Join the targets. targets must be a list instance with all targets
# you want to use for MCMC Bayesian inversion.
targets = Targets.JointTarget(targets=[target1, target2, target3, target4])

#  -------------------------------------------------------  MCMC BAY INVERSION
optimizer = MCMC_Optimizer(targets, initparams=initparams, priors=priors, random_seed=None)
optimizer.mp_inversion(baywatch=False, nthreads=25)
