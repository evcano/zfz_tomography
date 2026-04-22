import os
from sanpy.base.functions import query_pairs, query_virtual_source
from sanpy.base.project_functions import load_project
from sanpy.util.plot import plot_correlations
from sanpy.util.plot import plot_greens


# PARAMETERS
# =========
project_path = '/data/valeroe/red_sea_obs/noise_corr/zf_correlation.pkl'
data_path = "/data/valeroe/red_sea_obs/noise_corr/stacked_corr/all/ZZ"
data_format = 'sac'
data_type = 'correlations'

season = 'all'
bandpass = [1./10., 1./5.]
virtual_source = "ZF.OBS01"
conditions = []
apparent_velocity = False

maxlag = 300.0
global_normalization = False
yaxis = 'dis'
amplitude_only = False

# DONT EDIT BELOW THIS LINE
# =========================
P = load_project(project_path)

# query observations
if conditions:
    pairs = query_pairs(P, conditions)
else:
    pairs = P.pairs_list

if virtual_source:
    pairs = query_virtual_source(pairs, virtual_source)

# plot figures
if data_type == 'correlations':
    plot_correlations(data_path=data_path,
                      cmp="ZZ",
                      data_format=data_format,
                      pairs=pairs,
                      maxtime=maxlag,
                      bandpass=bandpass,
                      global_normalization=global_normalization,
                      yaxis=yaxis,
                      amplitude_only=amplitude_only,
                      apparent_velocity=apparent_velocity)

elif data_type == 'green':
    plot_greens(data_path=os.path.join(P.par['greens_path'], season),
                data_format=data_format,
                pairs=pairs,
                maxtime=maxlag,
                bandpass=[minfq, maxfq],
                global_normalization=global_normalization,
                yaxis=yaxis,
                amplitude_only=amplitude_only,
                apparent_velocity=apparent_velocity)
