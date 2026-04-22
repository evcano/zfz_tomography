# Tomography of the Zabargad Fracture Zone

This repository includes scripts and results related to:

- E. Valero-Cano, L. Parisi, H.A. Shiddiqi, S. Jónsson, and P.M. Mai, Transdimensional ambient-noise tomography of the Zabargad Fracture Zone, Red Sea (accepted for publication in Geophysical Journal International, 2026).

Download using:

  `git clone https://github.com/evcano/zabargad_fracture_zone_tomography.git`\
  `cd ./zabargad_fracture_zone_tomography`\
  `git lfs pull`

Figures shown in the paper can be reproduced using the scripts in `./figures_paper`. Some figures require data provided by third parties. Please read the `README.TXT` files in each directory for more information.

## Data locations

- Seismic noise cross-correlations: `./noise_correlations/stacked_correlations_sorted_west_to_east`
- Surface wave dispersion curves: `./dispersion_curves/disp_cur_??/final_curves`
- Estimated 1-D shear-wave velocity model: `./1d_inversion/model_mean.dat`
- Estimated 3-D shear-wave velocity model: `./figures_paper/3d_model_figures/average_model_ascii.txt`

## Methods

1-D tomography was conducted using BayHunter (Dreiling et al., 2019):

- https://github.com/jenndrei/BayHunter/tree/af06830c266fb9790e4c4c59503c08aed2cf8023

3-D tomography was conducted using a modified version of MCTomo (Zhang et al., 2018):

- https://github.com/evcano/MCTomo/tree/zafran

## Citation
If you use data, code, or results from this repository, please cite the corresponding publication listed above.
For questions, contact: eduardo.valero.cano@gmail.com

## References

- Dreiling, J., & Tilmann, F. (2019). BayHunter – McMC transdimensional Bayesian inversion of receiver functions and surface wave dispersion. GFZ Data Services.

- Zhang, X., Curtis, A., Galetti, E., & de Ridder, S. (2018). 3-D Monte Carlo surface wave tomography. Geophysical Journal International, 215(3), 1644–1658.
