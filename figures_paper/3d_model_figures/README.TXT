SCRIPTS

./get_3d_model_in_ascii.py
	Reads estimated model in MCTomo format and writes it to an ASCII file.

./get_3d_model_in_numpy.py
	Reads estimated model in MCTomo format and writes it to a numpy file.

./plot_slices_paper.py
	Plots Figure 10 and A2 of the paper.

./plot_study_area.py
	Contains functions used by plot_slices_paper.py

FILES

./average_model_ascii.txt
	Estimated 3-D Vs model in ASCII format.

./average_model.npy
	Estimated 3-D Vs model in numpy format.
	
./dep_mesh.npy
	Depth (in km) of the estimated Vs model. Use this file to plot average_model.npy .

./lat_mesh.npy
	Latitude of the estimated Vs model. Use this file to plot average_model.npy .

./lon_mesh.npy
	Longitude of the estimated Vs model. Use this file to plot average_model.npy .

./std_model_ascii.txt
	Standard deviation of the estimated 3-D Vs model in ASCII format.

./std_model.npy
	Standard deviation of the estimated 3-D Vs model in numpy format.
	Coordinates are indicated by lat_mesh.npy, lon_mesh.npy and dep_mesh.npy .

