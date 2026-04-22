import matplotlib.pyplot as plt
import numpy as np
import obspy as obs
import os
import shutil
from mpi4py import MPI
from glob import glob
from obstools.atacr.utils import rotate_dir
from obspy.signal.invsim import cosine_taper
from scipy import interpolate
from scipy import signal

# -------------------------------------------------
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nproc = comm.Get_size()

data_path = "/data/valeroe/red_sea_obs/data_processed/ZF"
tf_path = "/data/valeroe/red_sea_obs/correct_obs_tilt/transfer_functions"
outdir = "/data/valeroe/red_sea_obs/data_tilt_clock_corrected"

dt = 0.5  # dt of data to correct
dt_tf = 0.5 # dt of data used to get transfer functions
tf_duration = 1800.0  # duration of transfer function in time

tf_taper = None
co_threshold = 0.5
dont_weight_coherence = ["OBS08"]

stations_list = [f"OBS{x:02}" for x in range(1, 13)]
stations_list.extend(['NORTH','SOUTH'])

# EDIT WITH CAUTION
# -------------------------------------------------
sta = stations_list[myrank]
print(myrank,sta)

outdir = os.path.join(outdir,"ZF",sta)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# vertical component data
wfiles = glob(os.path.join(data_path, sta, f"*.*HZ.*.mseed"))
wfiles.sort()

for fileZ in wfiles:
    outfile = os.path.join(outdir, os.path.basename(fileZ))
    shutil.copyfile(src=fileZ, dst=outfile)

    # read waveforms
    file1 = fileZ.replace("HZ","H1")
    file2 = fileZ.replace("HZ","H2")

    if not os.path.isfile(file1) or not os.path.isfile(file2):
        print("missing one component")
        continue
    else:
        st = obs.Stream()
        st += obs.read(fileZ, format="mseed")
        st += obs.read(file1, format="mseed")
        st += obs.read(file2, format="mseed")

    # make sure all traces span the same time range
    st.sort(keys=["starttime"])
    sutc = st[-1].stats.starttime  # latest starttime
    day = sutc.date.isoformat()

    st.sort(keys=["endtime"])
    eutc = st[0].stats.endtime  # earliest endtime

    if sutc > eutc:
        print("trace is too short")
        continue

    st = st.slice(starttime=sutc, endtime=eutc)

    # detrend and taper
    st.detrend("linear")
    st.detrend("demean")
    st.taper(0.1)

    # get tr by component
    trZ = st.select(component="Z")[0]
    tr1 = st.select(component="1")[0]
    tr2 = st.select(component="2")[0]

    # fft traces
    tr_npts = trZ.stats.npts
    tr_fq = np.fft.fftfreq(tr_npts, d=dt)

    ftZ = np.fft.fft(trZ.data, n=tr_npts)
    ft1 = np.fft.fft(tr1.data, n=tr_npts)
    ft2 = np.fft.fft(tr2.data, n=tr_npts)

    # read power and cross spectral densities
    cZZ_file = os.path.join(tf_path, sta, f"cZZ_{day}.npy")

    if not os.path.isfile(cZZ_file):
        print("missing transfer function")
        continue
    else:
        cZZ = np.load(os.path.join(tf_path, sta, f"cZZ_{day}.npy"))
        cHH = np.load(os.path.join(tf_path, sta, f"cHH_{day}.npy"))
        cHZ = np.load(os.path.join(tf_path, sta, f"cHZ_{day}.npy"))
        tilt, _ = np.load(os.path.join(tf_path, sta, f"tilt_coh_{day}.npy"))

    # frequency axis of transfer function
    tf_npts = int(tf_duration/dt_tf)
    tf_fq = np.fft.fftfreq(tf_npts, d=dt_tf)

    # compute coherence and phase
    Co = np.abs(cHZ)**2 / (cHH*cZZ)
    Ph = np.angle(cHZ/cHH)

    # weight the coherence as in Tian et al
    if sta not in dont_weight_coherence:
        Co = Co * np.abs(np.cos(Ph))

    # smooth weighted coherence
    win = signal.windows.hann(20)
    Co = signal.convolve(Co, win, mode="same") / sum(win)

    # compute transfer function
    TF_ZH = np.conj(cHZ)/cHH

    # mute transfer function according to weighted coherence
    TF_ZH[Co < co_threshold] = 0.0

    # adjust transfer function length
    itp_real = interpolate.interp1d(tf_fq, np.real(TF_ZH), fill_value="extrapolate")
    itp_im = interpolate.interp1d(tf_fq, np.imag(TF_ZH),  fill_value="extrapolate")
    TF_ZH_2 = itp_real(tr_fq) + 1j * itp_im(tr_fq)

    # taper transfer function
    if tf_taper:
        ineg = tr_fq < 0
        ipos = tr_fq >= 0

        fq_pos = tr_fq[ipos]
        taper_pos = cosine_taper(npts=fq_pos.size, freqs=fq_pos, flimit=tf_taper)

        taper = np.zeros(tr_fq.shape)
        taper[ipos] = taper_pos
        taper[ineg] = taper_pos[:sum(ineg)][::-1]

        TF_ZH_2 = TF_ZH_2 * taper

    # remove tilt noise from vertical component
    ftH = rotate_dir(ft1, ft2, tilt)
    corrspec = ftZ - TF_ZH_2 * ftH
    corrtime = np.real(np.fft.ifft(corrspec))

    # save denoised trace
    corrtime = obs.Trace(data=corrtime, header=trZ.stats)
    corrtime.data = corrtime.data.astype("float32")
    corrtime.write(outfile,format="mseed")

print("done", sta)
