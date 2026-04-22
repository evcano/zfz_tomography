SCRIPTS

./estimate_clock_drift.py
	Estimates the clock drift of a station using the procedure described in
	Hannemann, K., Kruger, F., & Dahm, T., 2013. Measuring of clock drift rates and static time offsets of ocean bottom stations by means of ambient noise, Geophysical Journal International, 196(2), 1034–1042.
	
./fix_linear_clock_drift.py
	Uses linear interpolation to fix the clock drift of each station.

FILES

./*_clock.yml
	Input parameters for the estimate_clock_drift.py script.
	
./obs_skew.txt
	Measured clock skew of stations.
	
./obs_estimated_clock_drift.txt
	Estimated clock drift of each station.
 
./obs_drift_estimation_order.txt
	The order in which the clock drift of each station was estimated.
