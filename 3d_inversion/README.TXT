3-D inversion was conducted using a modified version of MCTomo (https://github.com/evcano/MCTomo/tree/zafran).

DIRECTORIES

./acausal_correlations_dat
	Symlink. Contains noise correlations used to obtain dispersion curves.
	Files are in MATLAB format, as read by the EGFAnalysisTimeFreq program.

./disp_cur_TT
	Contains dispersion curves in EGFAnalysisTimeFreq format

./disp_cur_ZZ
	Contains dispersion curves in EGFAnalysisTimeFreq format

./final_chains_2million
	Contains MCTomo inversion files. There is one directory for each of the 20 Markok Chains.
	Read the MCTomo manual to understand the meaning of each subdirectory.

./initial_models
	Contains initial 1-D models used for each Markov Chain.
	They are the same models as in ../1d_inversion/mean_models.
	However, the first line indicates the number of rows and columns as required by MCTomo.

./mctomo_inputfiles
	MCTomo input files. These files were used for all 20 Markov Chains.

./observed_times
	Contains frequency-dependent traveltimes used for MCTomo inversion.

SCRIPTS

./convert_dispersion_curves_to_traveltime.m
	Convert dispersion curves in ./disp_cur_* from velocity to traveltimes

./get_receivers_file.py
	Reads coordinates of receivers and generates sreceivers.dat required by MCTomo.

./get_traveltimes_files.py
	Reads files in ./observed_times and generates otimes_*.dat files required by MCTomo.

./get_water_depth_file.py
	Reads study area bathymetry (downloaded from GEBCO), interpolates to a regular grid,
	applies Gaussian smoothing, and saves the result to waterFile.dat, which is read by
	MCTomo.

./plot_data_noise_hist.py
	Plot histogram of the inverted data-noise hyperparameters.

./plot_mcmc_samples.py
	Plots negative likelihood, misfit, and number of cells for each Markov Chain (Figure 9 of the paper).

./plot_model_3d.py
	Reads estimated model (mean of posterior distribution) in MCTomo format and plots a 3-D figure
	using plotly.

./plot_model_3d_slices.py
	Reads estimated model (mean of posterior distribution) in MCTomo format and plots horizontal and
	vertical slices using plotly.


./plot_model_3d_slices_movie.py
	Reads the estimated model (mean of posterior distribution) in MCTomo format and plots an animation
	showing horizontal and vertical slices using plotly.

FILES

./stations_coordinates.txt
	Coordinates of stations used in this study.
