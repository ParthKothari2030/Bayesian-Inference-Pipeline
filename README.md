# Bayesian Inference pipeline
This repository contains C++-based simulation modules [Murmu et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract), C-based power spectrum computation tools [Murmu et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract), and Python scripts that integrate these components into an [```emceee```](https://emcee.readthedocs.io/en/stable/)-driven MCMC pipeline for cosmological data analysis using Bayesian inference.

## LIM SIM and New Power spec
The intensity and power spectrum framework has been updated for the pipeline such that it can take specific parameter values as input to generate maps, and the power spectrum is calculated for a specified range of k values. The updated code and its dependencies are demonstrated here for the [LIM simulator](https://github.com/ParthKothari2030/LIM_simulator) and [PowerSpectrum](https://github.com/ParthKothari2030/PowerSpectrum).

## Inference codes
Four modules have been developed to integrate the LIM simulator and Powerspectrum for inference using MCMC.
- `CodeExecutor`: Python subprocess wrapper for executing compiled C/C++ executables.
- `COMAPMeerKAT_Data`: Data loader for COMAP and MeerKAT power spectrum observations used in cosmological parameter estimation.
- `CosmologicalBayesianMethods`: MCMC-ready likelihood calculators for cosmological power spectrum analysis and parameter inference.
- `Main`: [```emceee```](https://emcee.readthedocs.io/en/stable/) based MCMC code for parameter estimation.
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
### Setting up C++ and C execution 

Adjust the ```paths_HI``` and ```paths_CO``` to store the maps in a directory from which the Powerspectrum code can access the map to compute the power spectrum. It is ideal to store the calculated Powerspectrum file, ```.powspec```, to be stored where ```.py``` files are executed, to reduce latency in reading.  

### MCMC
After activating the environment and providing all the paths, you are ready to use the pipeline.
```bash
python3 Main.py
```
and hit Enter. The MCMC parameter estimation is now up and running. 
**NOTE**: Please refer to [```emceee```](https://emcee.readthedocs.io/en/stable/) docs to understand how the ```Main.py``` utilized emceee to perform inference. Also refer to its chain saving process via ```backend``` process. This backend is also enabled in this pipeline.
