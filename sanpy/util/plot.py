import matplotlib.pyplot as plt
import numpy as np
import os
from obspy import read, Stream
from sanpy.util.apparent_velocity import compute_apparent_velocity


def plot_correlations(data_path, cmp, data_format, pairs=None,
                      maxtime=None, gain=1.0, alpha=1.0,
                      bandpass=None, global_normalization=False, yaxis=None,
                      amplitude_only=False, apparent_velocity=False, save=None,
                      branch="both", figsize=(6,9), fs=8):

    # read data
    st = Stream()

    if len(pairs) > 0:
        files = [f"{x}_{cmp}.{data_format}" for x in pairs]

        for f in files:
            stpath = os.path.join(data_path, f)

            if os.path.isfile(stpath):
                st += read(stpath, format=data_format)
    else:
        stpath = os.path.join(data_path, '*')
        st += read(stpath, format=data_format)

    ntr = len(st)

    # filter data
    if bandpass:
        st.detrend("linear")
        st.detrend("demean")
        st.taper(0.1)
        st.filter('bandpass',
                  freqmin=bandpass[0],
                  freqmax=bandpass[1],
                  corners=3,
                  zerophase=True)

    # cut data
    if maxtime:
        maxlag = ((st[0].stats.npts - 1) / 2) * st[0].stats.delta

        for i in range(0, ntr):
            st[i] = st[i].slice(st[i].stats.starttime + maxlag - maxtime,
                                st[i].stats.starttime + maxlag + maxtime)

    maxtime = ((st[0].stats.npts - 1) / 2) * st[0].stats.delta
    lags = np.linspace(-maxtime, maxtime, st[0].stats.npts)


    # normalize data
    if global_normalization:
        st.normalize(global_max=True)
    else:
        st.normalize(global_max=False)

    # sort data according to interstation distance
    distances = []
    for tr in st:
        distances.append(tr.stats.sac.dist)

    idx = np.argsort(np.array(distances))
    distances = np.sort(distances)
    print(distances[0], distances[-1])
    data = np.zeros((ntr, st[0].stats.npts))
    ax2_labels = []

    for i, j in enumerate(idx):
        data[i, :] = st[j].data
        ax2_labels.append(f"{st[j].stats.sac.kevnm}_{st[j].stats.sac.kstnm}")

    if branch == "causal":
        data = data[:,np.where(lags>=0.0)]
        data = data.reshape(data.shape[0], data.shape[2])
        lags = lags[lags>=0.0]
    elif branch == "acausal":
        data = data[:,np.where(lags<=0.0)]
        data = data.reshape(data.shape[0], data.shape[2])
        data = np.fliplr(data)  # time reverse
        lags = lags[lags<=0.0]
        lags = -lags[::-1]  # time reverse

    # setup figure
    plt.rcParams.update({'font.size': fs})

    fig, ax = plt.subplots(figsize=figsize)

    # labels
    if bandpass:
        ax.set_title(
            f'noise correlations {cmp.upper()} {1/bandpass[1]:.3f} - {1/bandpass[0]:.3f} s')
    else:
        ax.set_title(f'noise correlations {cmp.upper()}')

    ax.set_xlabel('lag [s]')

    if yaxis and yaxis == 'dis' and amplitude_only is False:
        ax.set_ylabel('interstation distance [km]')
    else:
        ax.set_ylabel('unitless')

    # plot data
    if amplitude_only:
        ax.imshow(data, extent=[lags[0], lags[-1], 0, ntr-1])
    else:
        offset = 0

        for i in range(0, ntr):
            if "KL" in ax2_labels[i]:
                c = "b"
            else:
                c = "k"

            y = data[i,:] * gain

            if yaxis and yaxis == 'dis':
                y += distances[i]
            else:
                y += offset
                offset = np.max(y)

            ax.plot(lags, y, c=c, lw=0.7, alpha=alpha)
            # station-pair name label
            ax.text(maxtime*0.9, np.mean(y), ax2_labels[i])


    if apparent_velocity:
        if yaxis and yaxis == 'dis':
            for s in apparent_velocity:
                plt.plot(s*np.array(distances), distances, 'r', lw=0.4,
                         alpha=alpha)
                plt.text(s*distances[-1], distances[-1], f"{1/s:.2f} km/s",
                         alpha=alpha)
                if branch == "both":
                    plt.plot(-s*np.array(distances), distances, 'r',
                             lw=0.4,alpha=alpha)
                    plt.text(-s*distances[-1], distances[-1], f"{1/s:.2f} km/s",
                             alpha=alpha)

    print('{} correlations plotted'.format(ntr))

    ax.set_xlim(lags[0], lags[-1])

    if save:
        plt.savefig(save,dpi=300)
    else:
        plt.show()
    return


def plot_greens(data_path, cmp, data_format, pairs=None, maxtime=None,
                bandpass=None, global_normalization=False, yaxis=None,
                amplitude_only=False, apparent_velocity=False):

    # read data
    st = Stream()

    if pairs:
        files = [f"{x}_{cmp}.{data_format}" for x in pairs]

        for f in files:
            stpath = os.path.join(data_path, f)

            if os.path.isfile(stpath):
                st += read(stpath, format=data_format)
    else:
        stpath = os.path.join(data_path, '*')
        st += read(stpath, format=data_format)

    ntr = len(st)

    # filter data and normalize
    if bandpass:
        st.detrend("linear")
        st.detrend("demean")
        st.taper(0.1)
        st.filter('bandpass',
                  freqmin=bandpass[0],
                  freqmax=bandpass[1],
                  corners=2,
                  zerophase=True)

    # cut maximum time
    if maxtime:
        for i in range(0, ntr):
            st[i] = st[i].slice(st[i].stats.starttime,
                                st[i].stats.starttime + maxtime)

    times = st[0].times()

    # normalize data
    if global_normalization:
        st.normalize(global_max=True)
    else:
        st.normalize(global_max=False)

    # sort data according to interstation distance
    distances = []
    for tr in st:
        distances.append(tr.stats.sac.dist)

    idx = np.argsort(np.array(distances))
    distances = np.sort(distances)

    data = np.zeros((ntr, st[0].stats.npts))
    for i, j in enumerate(idx):
        data[i, :] = st[j].data

    # estimate apparent velocity
    if apparent_velocity:
        c1, c2 = compute_apparent_velocity(data, times, distances)

    # setup figure
    fig, ax = plt.subplots()

    ax.set_title(f"Empirical Green's functions {cmp.upper()}")
    ax.set_xlabel('Time [s]')

    if yaxis and yaxis == 'dis' and amplitude_only is False:
        ax.set_ylabel('Interstation distance [km]')
    else:
        ax.set_ylabel('Unitless')

    # plot data
    if amplitude_only:
        ax.imshow(data, extent=[times[0], times[-1], 0, ntr-1])
    else:
        offset = 0

        for i in range(0, ntr):
            if yaxis and yaxis == 'dis':
                data[i, :] += distances[i]
            else:
                data[i, :] += offset
                offset = np.max(data[i, :])

            ax.plot(times, data[i, :], c='k', lw=0.5, alpha=0.5)

    print("{} empirical Green's functions plotted".format(len(idx)))

    if apparent_velocity:
        print('Apparent velocity: {} km/s'.format(1.0/c1))

        if yaxis and yaxis == 'dis':
            plt.plot(c1*np.array(distances)+c2, distances, 'r')

    plt.show()

    return
