# Bayesian Inference pipeline
This repository contains C++-based simulation modules [Murmu et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract), C-based power spectrum computation tools [Murmu et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract), and Python scripts that integrate these components into an ```emcee```-driven MCMC pipeline for cosmological data analysis using Bayesian inference.

## LIM SIM and New Power spec
The intensity and power spectrum framework has been updated for the pipeline such that it can take specific parameter values as input to generate maps, and the power spectrum is calculated for a specified range of k values. The updated code and its dependencies are demonstrated here for the [LIM simulator](https://github.com/ParthKothari2030/LIM_simulator) and [PowerSpectrum](https://github.com/ParthKothari2030/PowerSpectrum).

## Inference codes
Four modules have been developed to integrate the LIM simulator and Powerspectrum for inference using MCMC.
- `CodeExecutor`:
