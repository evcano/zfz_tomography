OBS tilt was corrected using OBStools.

SCRIPTS

./compute_vertical_tilt_transfer_function.py
 	Computes the transfer function between the vertical and horizontal components rotated in the direction
 	of maximum tilt using OBStools.
	
./correct_vertical_tilt.py
	Uses transfer functions computed by compute_vertical_tilt_transfer_function.py to remove tilt
	noise from the horizontal components.
