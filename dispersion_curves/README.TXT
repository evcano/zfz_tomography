NOTES
In this directory, the words "coast, island, and obs" indicate measurements involving ocean-bottom seismometers (OBSs) and land stations, OBSs and island stations, and OBSs and OBSs, respectively.

A modified version of the program EGFAnalysisTimeFreq (originally developed by Huajian Yao, USTC; available at https://github.com/ShaoboYang-USTC/EGFAnalysisTimeFreq )
was used to pick dispersion curves. As no license is provided with the original repository, the modified code is not redistributed. For further information about the modifications,
please contact eduardo.valero.cano@gmail.com 

DIRECTORIES

./average_disp_img
	Contains stacked group- and phase-velocity diagrams used to obtain group-velocity bounds
	and phase-velocity reference curves (see Figures 3 and 4 of the paper).
	There are 3 files per diagram: taxis_*.npy (horizontal axis), vaxis_*.npy (vertical axis) and S_*.npy (diagram).
	
./disp_cur_TT (also ./disp_cur_ZZ)
	Contains dispersion curves.
		./disp_cur_TT/automatically_picked: output of EGFAnalysisTimeFreq.
			First row: station 2 longitude, latitude
			Second row: station 1 longitude, latitude
			Columns: filter central period [s], observed period [s], velocity [km/s], SNR [dB], quality (0 = bad, 1 =  good)
			
		./disp_cur_TT/manually_checked: output of check_disp_curves_manually.py
			Columns: period [s], velocity [km/s], SNR [dB]
			
		./disp_cur_TT/final_curves: output of remove_outlier_measurements.py . These are the dispersion curves used for 1-D and 3-D inversions.
			Columns: period [s], velocity [km/s]

./disp_img_TT (also ./disp_img_ZZ)
	Contains dispersion diagrams shown in Figures 3 and 4 of the paper.
	Plot using ../figures_paper/plot_disp_images.py

./mean_curves
	Mean of all dispersion curves after quality control.
	Columns: period [s], velocity [km/s]
	
./reference_curves
	Contains group-velocity bounds and phase-velocity reference curves used to guide EGFAnalysisTimeFreq.
	*-refcurve.txt: Reference curve. Columns: period [s], velocity [km/s], variance [km/s] (plus/minus 0.5 km/s from reference, not used in practice).
	*-twindow.txt: Velocity bounds. Columns: period [s], middle bound [km/s], lower bound [km/s], upper bound[km/s]
	
SCRIPTS

./check_disp_curves_manually.py
	Code to manually inspect dispersion curves. It also removes low-quality measurements based on SNR and interstation distance.
	
./pick_reference_curves.py
	Code to manually pick group-velocity bounds and phase-velocity reference curves.
	
./remove_outlier_measurements.py
	Removes measurements outside the specified period range and those deviating more than the specified standard deviation.
	
./stack_ftan_images.py
	Stacks dispersion diagrams (group or phase) written by EGFAnalysisTimeFreq. The stacked diagrams are used to pick
	group-velocity bounds and phase-velocity reference curves.
		
FILES

./stations_coordinates.txt
	Station's coordinates (latitude, longitude).
