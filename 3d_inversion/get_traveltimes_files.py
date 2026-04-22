import matplotlib.pyplot as plt
import numpy as np
import os


def read_otimes(fname, paxis, paxis_keep):
    outarr = np.zeros(paxis.size)
    outarr[:] = -1  # not available traveltimes are marked with -1

    X = np.loadtxt(fname)
    obs_per = X[:,0]
    obs_tim = X[:,1]

    if obs_per.size != np.unique(obs_per).size:
        raise Exception(
            "there is more than 1 measurement for the same period"
        )

    for i in range(obs_per.size):
        if obs_per[i] in paxis_keep:
            idx = np.where(paxis == obs_per[i])
            outarr[idx] = obs_tim[i].copy()

    return outarr


stafile = "./stations_coordinates.txt"
datapath = "./observed_times"
outpath = "./mctomo_inputfiles"

cmpts = ["zz", "tt"]
mtypes = ["group", "phase"]

zz_grp_pmin = 3
zz_grp_pmax = 11

zz_pha_pmin = 4
zz_pha_pmax = 10

tt_grp_pmin = 3
tt_grp_pmax = 12

tt_pha_pmin = 3
tt_pha_pmax = 11

# ------------------------------------
net_list = np.loadtxt(stafile, usecols=(0), dtype="U20")
sta_list = np.loadtxt(stafile, usecols=(1), dtype="U20")
sta_codes = []
for i in range(len(net_list)):
    sta_codes.append(f"{net_list[i]}.{sta_list[i]}")

for cmp_ in cmpts:
    cmp_ = cmp_.upper()

    if cmp_ == "ZZ":
        pmin = min(zz_pha_pmin, zz_grp_pmin)
        pmax = max(zz_pha_pmax, zz_grp_pmax)
    elif cmp_ == "TT":
        pmin = min(tt_pha_pmin, tt_grp_pmin)
        pmax = max(tt_pha_pmax, tt_grp_pmax)
    else:
        raise Exception("wrong component")

    paxis = np.arange(pmin, pmax + 1)
    fqaxis = 1./paxis
    nfq = len(fqaxis)

    datapath2 = os.path.join(datapath, cmp_)
    print("reading data from ", datapath2)

    for mt in mtypes:
        mt = mt.upper()

        if mt == "GROUP":
            curveprefix = "GDisp"
        elif mt == "PHASE":
            curveprefix = "CDisp.T"

        if cmp_ == "ZZ" and mt == "GROUP":
            paxis_keep = np.arange(zz_grp_pmin, zz_grp_pmax + 1)
        elif cmp_ == "ZZ" and mt == "PHASE":
            paxis_keep = np.arange(zz_pha_pmin, zz_pha_pmax + 1)
        elif cmp_ == "TT" and mt == "GROUP":
            paxis_keep = np.arange(tt_grp_pmin, tt_grp_pmax + 1)
        elif cmp_ == "TT" and mt == "PHASE":
            paxis_keep = np.arange(tt_pha_pmin, tt_pha_pmax + 1)

        outfile = os.path.join(outpath, f"otimes_{cmp_}_{mt}.dat")

        counter1 = 0
        counter2 = np.zeros(paxis.size)

        with open(outfile, "w") as _file:
            _file.write(f"{nfq}\n")  # number of frequencies

            for fq in fqaxis:
                _file.write(f"{fq:.3f} ")  # frequency

            _file.write("\n")

            for src in sta_codes:
                for rec in sta_codes:
                    fname = f"{curveprefix}.{src}_{rec}_{cmp_}.dat"

                    if not os.path.exists(os.path.join(datapath2,fname)):
                        _file.write("0\n")
                        continue
                    else:
                        otimes = read_otimes(os.path.join(datapath2,fname), paxis, paxis_keep)

                        _file.write("1\n")

                        for otime in otimes:
                            _file.write(f"{otime:.6f} 0.01\n")  # observed traveltimes

                        counter1 += 1
                        counter_tmp = np.ones(otimes.size)
                        counter_tmp[np.where(otimes==-1)] = 0
                        counter2 += counter_tmp


        print(f"{cmp_} {mt}")
        print(f"periods: ", paxis)
        print(f"no files read: ", counter1)
        print(f"no measurements: ", counter2, "\n")
        print(np.sum(counter2))
