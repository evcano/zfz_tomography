import matplotlib.pyplot as plt
import numpy as np
import os
from obspy import read, Stream
from obspy.io.sac.sactrace import SACTrace
from glob import glob

cmp = "TT"
data_path = f"../noise_correlations_tomo/stack_tfpws/all/{cmp}"
output_path = f"./correlations_west2east/{cmp}"
fmt = "sac"

#  EDIT WITH CAUTION
# =======================================================
files = glob(os.path.join(data_path, "*"))

for fil in files:
    fname = os.path.basename(fil)
    s1, s2, suffix = fname.split("_")

    # acausal branch contains waves from s1 to s2
    # causal branch contains waves from s2 to s1
    st = read(fil,fmt=fmt)
    tr = st[0]

    if s1 != tr.stats.sac['kevnm']:
        raise Exception

    if s2 != tr.stats.sac['kstnm']:
        raise Exception

    s1lon = tr.stats.sac.evlo
    s1lat = tr.stats.sac.evla

    s2lon = tr.stats.sac.stlo
    s2lat = tr.stats.sac.stla

    # reverse correlation
    # acausal branch  waves from WEST to EAST
    # causal branch waves from EAST to WEST
    if s1lon < s2lon:
        # we do nothing since acausal already has waves from s1 to s2
        west_station = s1
        newfname = fname
    else:
        west_station = s2
        tr.data = tr.data[::-1]

        tr = SACTrace.from_obspy_trace(trace=tr, keep_sac_header=True)
        tr.kevnm = s2
        tr.kstnm = s1

        newfname = f"{s2}_{s1}_{suffix}"

    tr.write(os.path.join(output_path, newfname))
