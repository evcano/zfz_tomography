DIRECTORIES

./bad_windows
	Contains files indicating the date and time of transient signals recorded by each station.
	These files are used to avoid transient signals during the computation of ambient noise correlations.
	
SCRIPTS

./compute_running_max_amplitude.py
	For each station and component, slides a time window and computes the maximum amplitude.
	Used to detect transient signals.
	
./determine_bad_windows.py
	Detects time windows containing transient signals and writes output to ./bad_windows .
	Needs files written by ./compute_running_max_amplitude.py
