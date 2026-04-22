import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from pathlib import Path
from matplotlib import rc

rc("text", usetex=True)

dirs = ["./all_modes", "./mode_0", "./mode_1", "./mode_2"]
labels = ["All modes", "Mode 0", "Mode 1", "Mode 2"]

f1 = 0.0 # 1 = north
f2 = 0.0 # 1 = east
f3 = 1.0 # 1 = down
theta = 1.0 # measured clockwise from the north

fmin = 0.083
fmax = 0.3
 
fig1, axes1 = plt.subplots(4,2,figsize=(8,5), constrained_layout=True)

# loop over modes
for i, dir_ in enumerate(dirs):
    datadir = Path(dir_)
    
    st = read(datadir / "*")
    st.detrend("demean")
    st.detrend("linear")
    st.taper(0.1)
    st.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)
    
    tr_zvf = st.select(channel="ZVF")
    tr_zhf = st.select(channel="ZHF")
    tr_rvf = st.select(channel="RVF")
    tr_rhf = st.select(channel="RHF")
    tr_thf = st.select(channel="THF")
    
    print(tr_zvf[0])
    times = tr_zvf[0].times()
    
    #z = (f1*np.cos(theta) + f2*np.sin(theta)) * tr_zhf[0].data + f3*tr_zvf[0].data
    #r = (f1*np.cos(theta) + f2*np.sin(theta)) * tr_rhf[0].data + f3*tr_rvf[0].data
    #t = (f1*np.sin(theta) - f2*np.cos(theta)) * tr_thf[0].data

    z = tr_zvf[0].data
    t = tr_thf[0].data
  
    if i == 0: 
        zmax = np.max(np.abs(z))
        tmax = np.max(np.abs(t))
    if i == 1:
        zmax2 = np.max(np.abs(z[times<=90]))
        tmax2 = np.max(np.abs(t[times<=90]))
    if i > 1:
        print(zmax2 / np.max(np.abs(z[times<=90])))
        print(tmax2 / np.max(np.abs(t[times<=90])))

    axes1[i,0].plot(times, z, c="k", lw=1.2, alpha=0.9, label=labels[i])
    axes1[i,1].plot(times, t, c="k", lw=1.2, alpha=0.9, label=labels[i])

    axes1[i,0].set_ylim(-zmax, zmax)
    axes1[i,1].set_ylim(-tmax, tmax)

    axes1[i,0].legend()
    axes1[i,1].legend()

    if i != 3:
        axes1[i,0].set_xticklabels([])
        axes1[i,1].set_xticklabels([])
    else:
        axes1[i,0].set_xlabel("Time [s]")
        axes1[i,1].set_xlabel("Time [s]")

    axes1[i,0].set_ylabel("Amplitude")
    axes1[i,1].set_yticklabels([])

axes1[0,0].set_title(f"a) Z Green's function {1/fmax:.2f} - {1/fmin:.2f} s")
axes1[0,1].set_title(f"b) T Green's function {1/fmax:.2f} - {1/fmin:.2f} s")

fig1.savefig("zgf.png",dpi=200)
plt.show()
plt.close()
