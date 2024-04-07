# Total-Radiation-Belt-Electron-Content

MATLAB code used to compute the Total Radiation Belt Electron Content using particle flux measurements from the Van Allen Probes Magnetic Electron Ion Spectrometer (MagEIS) instrument.

## Brief File Description
The TRBEC_E_PA.m\TRBEC_mu_K.m files contain the necessary code to compute the TRBEC following the procedure outlined in Estimating Electron Precipitation using the Total Radiation Belt Electron Content (J. Pitzel 2022, https://dx.doi.org/10.11575/PRISM/40159). For the non-adiabatic bounds using E and PA the input data is the L3 MagEIS data avalible from the Coordinated Data Analysis Web (CDAWeb) or the Los Alamos National Laboratory (LANL). For adiabatic bounds using mu and K the required input data is the adiabatic PSD in adiabatic coordinates avalible from the Radiation Belt Storm Probes Gateway (RBSPGWAY). The code also requires two librarys to run, the first is an spdf reader from NASA (https://cdf.gsfc.nasa.gov/html/matlab_cdf_patch.html) and the second is the IRBEM magnetic field library (https://prbem.github.io/IRBEM/). Per half orbit the computation time of the TRBEC is on the order of seconds so a secondary function is to output an intermediate data product called the differential content (DC). The DC is the TRBEC integrand after integration over PA or K for all trapped partiles and can be used for computing multiple sets of integration bounds quickly.

The TRBEC_E_PA_from_DC.m/TRBEC_mu_K_from_DC.m files use the DC quickly compute the TRBEC given a set of integration bounds in either E,L* or mu,L*. The DC dataset has been precomputed in E/PA for the entire RBSP mission, but only for January-March 2013 for the adiabatic bounds.
