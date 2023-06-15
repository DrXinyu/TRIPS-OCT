Updated on 15 June 2023.
Fix the bugs 1) ddet computation should be comment. 2) unwrap optix axis is not on.

This is code for image reconstruction using triple input polarization sensitive optical coherence tomography, decribed in a paper titled "Posterior scleral birefringence measured by triple-input polarization-sensitive imaging as a biomarker of myopia progression" published on Nature Biomedical Engineering, https://doi.org/10.1038/s41551-023-01062-w

Run TRIPS_Mueller_reconstruction.m to evaluate the code for TRIPS-OCT Mueller matrix reconstruction.
Run TRIPS_PMD_and_symmetrization.m to evaluate the code for spectral binning to mitigate PMD and recover symmetry.

Detailed funtions in comments.

Data: 
One typical guinea pig retina frame in-vivo
S123.mat is the Stokes vectors of the B-scan to be reconstructed.
DATA LINK: 
https://www.dropbox.com/s/l3me1zog2npejyg/S123.mat?dl=0

bigS123.mat is the Stokes vectors to generate calibration metrics, randomly selected from the volume scan.
DATA LINK:
https://www.dropbox.com/s/rlr5cwscj796szt/bigS123.mat?dl=0


System requirement:
GPU: above GTX1080Ti is recommended. 

Contact:
Xinyu Liu
liux0060@e.ntu.edu.sg

Martin Villiger
mvilliger@mgh.harvard.edu
