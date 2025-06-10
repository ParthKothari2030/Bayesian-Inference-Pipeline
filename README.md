# Bayesian Inference pipeline
This repository contains C++-based simulation modules [Murmu et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract), C-based power spectrum computation tools [Murmu et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract), and Python scripts that integrate these components into an ```emcee```-driven MCMC pipeline for cosmological data analysis using Bayesian inference.

## LIM SIM and New Power spec
The intensity and power spectrum framework has been updated for the pipeline such that it can take specific parameter values as input to generate maps, and the power spectrum is calculated for a specified range of k values. The updated code and its dependencies are demonstrated here for the [LIM simulator](https://github.com/ParthKothari2030/LIM_simulator) and [PowerSpectrum](https://github.com/ParthKothari2030/PowerSpectrum).

## Inference codes
Four modules have been developed to integrate the LIM simulator and Powerspectrum for inference using MCMC.
- `CodeExecutor`: Python subprocess wrapper for executing compiled C/C++ executables.
- `COMAPMeerKAT_Data`: Data loader for COMAP and MeerKAT power spectrum observations used in cosmological parameter estimation.
- `CosmologicalBayesianMethods`: MCMC-ready likelihood calculators for cosmological power spectrum analysis and parameter inference.
- `Main`: ```emcee``` based MCMC code for parameter estimation.
- `TimeCheck`: script to check the time taken for one iteration of map + powerspectrum calculation using Python wrapper.

## Pipeline usage
 A conda environment needs to be created before using this pipeline. ```.yaml``` files contain the dependencies for the creation of the environment called *Parameter_estimation*.
 To install conda on the local machine, refer to Conda [docs](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). 


 Once Conda is installed, use this command to create the environment:
``` bash 
conda env create -f environment.yml
```
Activating the environment:

``` bash
conda activate Parameter_estimation.
```

