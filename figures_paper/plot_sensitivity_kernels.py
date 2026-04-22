import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
from disba import PhaseDispersion
from disba import EigenFunction
from disba import PhaseSensitivity
from disba import GroupSensitivity
from matplotlib import rc


rc("text", usetex=True)


def get_mean_initial_model(mpath,nchains):
    for i in range(1,nchains+1):
        fname = os.path.join(mpath,f"model_{i:03}.dat")
        m = np.loadtxt(fname,usecols=(0,2),skiprows=1)
        if i == 1:
            mean = m.copy()
        else:
            mean[:,1] += m[:,1]

    mean[:,1] /= nchains
    return model


def get_layered_model(model):
    nlay = model.shape[0]

    h = np.zeros(nlay*2)
    vs = np.zeros(h.shape)

    a = 0
    for i in range(0,nlay):
        b = a + 1
        if i == 0:
            h[a] = 0.0
        else:
            h[a] = h[a-1]

        h[b] = h[a] + model[i,0]
        if model[i,2] != 0.0:
            vs[a] = model[i,2]
            vs[b] = model[i,2]
        else:
            vs[a] = np.nan
            vs[b] = np.nan

        a = b + 1

    return h, vs


mfile = "../1d_inversion/model_mean_layered.dat"

# Periods must be sorted starting with low periods
t = np.linspace(1.0, 10, 50)

figname = "1dkernels.png"
figsize = (6.5,3.5)
fs = 9
fs2 = 8
xticks = np.arange(2,8)
xticks2 = np.arange(-1,2.)
yticks = np.arange(0,24,2)

# --------------------------------------------------------------------
velocity_model = np.loadtxt(mfile)
# convert layer thickness to depth, for plotting purposes only
h, vs = get_layered_model(velocity_model)


# SENSITIVITY KERNES
# ---------------------------------------------------
# ps returns a namedtuple (depth, kernel, period, velocity, mode,wave, type, parameter)
ps = PhaseSensitivity(*velocity_model.T)
gs = GroupSensitivity(*velocity_model.T)

periods = np.arange(3,13)

fig, ax = plt.subplots(1,4)

base = cm.get_cmap(name="rainbow")
colist = base(np.linspace(0,1,len(periods)))

for i,per in enumerate(periods):
    parameter = "velocity_s"

    gr = gs(per, mode=0, wave="rayleigh", parameter=parameter)
    pr = ps(per, mode=0, wave="rayleigh", parameter=parameter)

    gl = gs(per, mode=0, wave="love", parameter=parameter)
    pl = ps(per, mode=0, wave="love", parameter=parameter)

    gr.kernel[-1] = gr.kernel[-2]
    pr.kernel[-1] = pr.kernel[-2]
    gl.kernel[-1] = gl.kernel[-2]
    pl.kernel[-1] = pl.kernel[-2]

    ax[0].plot(gr.kernel,gr.depth,label=f"T={gr.period} s",linewidth=0.8,c=colist[i])
    ax[1].plot(pr.kernel,pr.depth,label=f"T={pr.period} s",linewidth=0.8,c=colist[i])

    ax[2].plot(gl.kernel,gl.depth,label=f"T={gl.period} s",linewidth=0.8,c=colist[i])
    ax[3].plot(pl.kernel,pl.depth,label=f"T={pl.period} s",linewidth=0.8,c=colist[i])

ax[0].set_title(r"\textbf{(a)} Rayleigh waves" + "\ngroup-velocity",
                fontsize=fs, loc="left")

ax[1].set_title(r"\textbf{(b)} Rayleigh waves" + "\nphase-velocity",
                fontsize=fs, loc="left")

ax[2].set_title(r"\textbf{(c)} Love waves" + "\ngroup-velocity",
                fontsize=fs, loc="left")

ax[3].set_title(r"\textbf{(d)} Love waves" + "\nphase-velocity",
                fontsize=fs, loc="left")

ax[0].set_ylabel("Depth [km]")

ax[0].set_xlabel("dU/d$V_s$",fontsize=fs)
ax[1].set_xlabel("dc/d$V_s$",fontsize=fs)
ax[2].set_xlabel("dU/d$V_s$",fontsize=fs)
ax[3].set_xlabel("dc/d$V_s$",fontsize=fs)

for i in range(0,4):
    if i > 0:
        ax[i].tick_params(labelleft=False)

    ax[i].set_ylim(yticks[-1],yticks[0])

    ax[i].tick_params(labelsize=fs)

leg = ax[i].legend(fontsize=fs2,
                   framealpha=0.8,
                   edgecolor="k",
                   fancybox=False
                  )


fig.set_size_inches(figsize)
fig.tight_layout()
fig.savefig(figname,dpi=300)
