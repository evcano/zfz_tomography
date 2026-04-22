1-D inversion was conducted using BayHunter (https://github.com/jenndrei/BayHunter/tree/af06830c266fb9790e4c4c59503c08aed2cf8023).
To consider a water layer during the inversion, add the following lines after <def run_model(self, h, vp, vs, rho, **params):> in BayHunter/src/surf96_modsw.py:

	# water layer                                                                                                     
	h = np.append(1.2, h)
	vp = np.append(1.5, vp)
	vs = np.append(0.0, vs)
	rho = np.append(1.0, rho)
	nlayer = len(h)


DIRECTORIES

./mean_models
	Output of 'get_mean_model_per_chain.py
 
./observed_obs
	Contains observed dispersion curves used for the inversion.

./results/data
	Contains output data files produced by BayHunter during the inversion process. 
	Since they are large files (>100 GB), I did not upload them.
 
SCRIPTS
 
get_best_model_per_chain.py
	For each of the 20 Markov Chains with the highest likelihood, get the model with the minimum misfit.
 
get_mean_model.py
	Obtain the mean (or median and mode as well) model from all 100 Markov Chains.
 
get_mean_model_per_chain.py
	For each of the 20 Markov Chains with the highest likelihood, compute the mean model and store it in "./mean_models". 
	These models serve as initializations for the 20 Markov Chains of the 3-D inversion.

plot_summary.py
	BayHunter script to plot inversion results summary
 
run_inversion.py
	Runs BayHunter inversion

FILES

config.ini
	BayHunter ini file with inversion parameters

model_*.dat
	1-D velocity model (depth [km], vp [km/s], vs[km/s], rho [g/cm3]). It is computed from the final posterior distribution.


