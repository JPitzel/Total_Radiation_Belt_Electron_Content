[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11323838.svg)](https://doi.org/10.5281/zenodo.11323838)

# Total Radiation Belt Electron Content

For a defined electron population, the Total Radiation Belt Electron Content (TRBEC) represents the number of electrons present in the Van Allen Radiation Belts. This repository contains the necessary code to compute the Total Radiation Belt Electron Content (TRBEC) using particle flux measurements from the Van Allen Probes Magnetic Electron Ion Spectrometer (MagEIS) instrument, either by defining the electrons as a function of the non-adiabatic invariants of kinetic energy (E), pitch angle (PA), and Roederer L-value (L*), or the adiabatic invariants of the magnetic moment (μ or mu), second adiabatic invariant (K), and L*. A complete technical description of the methodology of this work is available in _Estimating Electron Precipitation using the Total Radiation Belt Electron Content_ (J. Pitzel 2022, https://dx.doi.org/10.11575/PRISM/40159) and is outlined more concisely in a submission (as of May 2024) to the _Journal of Geophysical Research: Space Physics_ as a "Technical Reports: Methods" paper.

<p align="center">
  <img src="https://github.com/JPitzel/Total_Radiation_Belt_Electron_Content/blob/main/images/single_particle_motion_illustration.png" width="750">
</p>

## File Description

The files contain the necessary MATLAB code to compute the TRBEC directly from MagEIS data products or by using the provided intermediate data product.

The 'TRBEC_E_PA.m' and 'TRBEC_mu_K.m' files compute the TRBEC in the appropriate coordinates using the following publicly available data. For non-adiabatic coordinates of E and PA, the input data is the level 3 MagEIS flux measurements available from the Coordinated Data Analysis Web (CDAWeb) (https://cdaweb.gsfc.nasa.gov/) or the Los Alamos National Laboratory (LANL) (https://rbsp-ect.newmexicoconsortium.org/data_pub/). For the adiabatic coordinates of μ and K, the required input data is the adiabatic PSD available from the Radiation Belt Storm Probes Gateway (RBSPGWAY) (https://rbspgway.jhuapl.edu/), as well as the non-adiabatic data for magnetic field information. The provided MATLAB code also requires the MATLAB curve fitting toolbox (for use of the "smooth" function) and NASA's CDF reader for MATLAB (https://cdf.gsfc.nasa.gov/html/matlab_cdf_patch.html). Each TRBEC data point is computed from each Van Allen Probes half orbit. The computation time for each point is on the order of seconds on a moderately equipped machine. Therefore, a secondary function is to output an intermediate data product called the Differential Content (DC) to more efficiently compute different sets of integration bounds. The DC is the intermediate TRBEC integrand after integration over PA or K for all trapped particle trajectories and is most useful for computing the TRBEC of multiple sets of particle populations more quickly.

The 'TRBEC_E_PA_from_DC.m' and 'TRBEC_mu_K_from_DC.m' files use the DC to quickly compute the TRBEC given a set of integration bounds, either in E, L* or μ, L*, and display the results. The DC dataset has been precomputed in E/L* for the entire Van Allen Probes mission, but only from January to March 2013 for the adiabatic bounds as required for the case study in the JGR article.
