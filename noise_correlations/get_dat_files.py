import matplotlib.pyplot as plt
import numpy as np
import os
from obspy import read, Stream
from glob import glob
from sanpy.util.correlation_branches import correlation_branches
from sanpy.base.project_functions import load_project
from sanpy.util.plot import plot_correlations
from scipy.io import savemat

# THIS OUTPUTS CORRELATION FUNCTIONS (CF)
# IT COPIES THE INPUT TRACE TO BOTH CORRELATION BRANCHES
# THUS, BOTH CF BRANCHES ARE THE SAME

cmp = "TT"
data_path = f"./acausal_correlations/{cmp}"
output_path = f"./acausal_correlations_dat/{cmp}"
staelv_file = "./stations_elevation.txt"

# ---------------------------------------------------
obs_list = [f"ZF.OBS{x:02}" for x in range(1,13)]
obs_list.extend(["ZF.NORTH", "ZF.SOUTH"])

island_list = ["ZF.BREEM", "ZF.QUMAN"]

coast_list = ["ZF.LAVA","ZF.KHUF","SA.EWJHS", "SA.WJHS"]

inland_list = ["KL.WEST","KL.EAST","KL.SOUTH","KL.DEEP"]
# ---------------------------------------------------

#  EDIT WITH CAUTION
# =======================================================
esta = np.loadtxt(staelv_file,usecols=0,dtype="U20")
esta_elv = np.loadtxt(staelv_file,usecols=1)

print(esta)
print(esta_elv)

files = glob(os.path.join(data_path,"*"))

corr_notobs  = []
corr_obs = []
corr_island = []
corr_coast = []
corr_inland = []

for fpath in files:
    fname = os.path.basename(fpath)

    st = read(fpath)
    tr = st[0]

    # stations coordinates
    s1 = tr.stats.sac.kevnm
    s1lon = tr.stats.sac.evlo
    s1lat = tr.stats.sac.evla

    ii = np.argwhere(esta == s1).flatten()
    s1elv = esta_elv[ii[0]]

    s2 = tr.stats.sac.kstnm
    s2lon = tr.stats.sac.stlo
    s2lat = tr.stats.sac.stla

    ii = np.argwhere(esta == s2).flatten()
    s2elv = esta_elv[ii[0]]

    # fill array
    col1 = np.zeros(tr.stats.npts+2)
    col2 = col1.copy()
    col3 = col1.copy()

    col1[0] = s1lon
    col1[1] = s2lon
    col1[2:] = tr.times()

    col2[0] = s1lat
    col2[1] = s2lat
    col2[2:] = tr.data

    col3[0] = s1elv
    col3[1] = s2elv
    col3[2:] = tr.data

    outarray = np.array([col1,col2,col3]).T
    outfile = fname.replace(".SAC",".dat")

    np.savetxt(os.path.join(output_path,outfile), outarray)

    if s1 not in obs_list and s2 not in obs_list:
        corr_notobs.append(outfile)
    elif s1 in obs_list and s2 in obs_list:
        corr_obs.append(outfile)
    elif s1 in island_list or s2 in island_list:
        corr_island.append(outfile)
    elif s1 in coast_list or s2 in coast_list:
        corr_coast.append(outfile)
    elif s1 in inland_list or s2 in inland_list:
        corr_inland.append(outfile)

corr_notobs.sort()
corr_obs.sort()
corr_island.sort()
corr_coast.sort()
corr_inland.sort()

np.savetxt(os.path.join(output_path,"corr_notobs.txt"), corr_notobs, fmt="%s")
np.savetxt(os.path.join(output_path,"corr_obs.txt"), corr_obs, fmt="%s")
np.savetxt(os.path.join(output_path,"corr_island.txt"), corr_island, fmt="%s")
np.savetxt(os.path.join(output_path,"corr_coast.txt"), corr_coast, fmt="%s")
np.savetxt(os.path.join(output_path,"corr_inland.txt"), corr_inland, fmt="%s")
