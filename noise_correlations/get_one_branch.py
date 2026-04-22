import matplotlib.pyplot as plt
import numpy as np
import os
from obspy import read, Stream
from glob import glob
from sanpy.util.correlation_branches import correlation_branches
from sanpy.base.project_functions import load_project
from sanpy.util.plot import plot_correlations


cmp = "TT"
branch = "neg"

data_path = f"./correlations_west2east/{cmp}"
output_path = f"./acausal_correlations/{cmp}"
fmt = "sac"

min_ntraces = 59

#  EDIT WITH CAUTION
# =======================================================
files = glob(os.path.join(data_path,"*"))

for fpath in files:
    fname = os.path.basename(fpath)
    s1, s2, _ = fname.split("_")

    st = read(fpath, fmt=fmt)
    tr = st[0]

    ntraces = tr.stats.sac.user7
    if ntraces < min_ntraces:
        print(f"skipping {fname}", ntraces)

    tr_branch =  correlation_branches(tr, branch=branch)

    ofile = os.path.basename(fpath)
    ofile = ofile.upper()
    ofile = os.path.join(output_path, ofile)
    tr_branch.write(ofile, fmt="SAC")
