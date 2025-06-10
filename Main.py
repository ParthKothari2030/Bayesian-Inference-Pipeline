import numpy as np 
import CodeExecutor as ce
import emcee
import CosmologicalBayesianMethods as cbm
import ComapMeerKAT_Data as cmd
import time


"""
Module name: Main.
Function: MCMC using open source library called "emcee".
Module requirement: Code executor, emcee, numpy, COMAP MeerKAT data, time.
MCMC funcationalities: Independent and Joint analysis.
"""



start = time.time()



def Independent_CO_Likelihood(params):


        alpha_ext, beta_ext = params
        alpha_ext, beta_ext = np.round(alpha_ext, 2) , np.round(beta_ext, 2)

        
        # Priors for Bayesian Inference
        if not (0.85 <= alpha_ext <= 1.49 and -2.94 <= beta_ext <= 2.15): # 3 sigma-priors from the parameter space
                return -np.inf

        else:
                lim_simulator_CO = ce.CodeExecutor("/path/to/CO_LIM_SIM_3.0",
                                        "LIM_simulator",
                                        "param.txt",
                                        "paths_CO.txt",
                                        str(alpha_ext), str(beta_ext))
                lim_simulator_CO.run()
        

                powerspectrum_CO = ce.CodeExecutor("/path/to/New_power_spec_3_CO",
                                        "PowerSpectrum",
                                        "input")
                powerspectrum_CO.run()
                
                # Fetching COMAP data
                COMAPData_class = cmd.COMAPPowerSpectrum(3.0).get_data()

                # Fetches Model Powerspec
                ModernPS_CO = cbm.TheoreticalPowerSpectrum('CO_JointAnalysis_0.56MPc_z3.010.bin.powspec').read_power_spectrum()
                ModernPS_CO = ModernPS_CO[:-2]
                
                # Calculates Likelihood
                Likelihood = cbm.MultiVariateLogGaussianLikelihood(ModernPS_CO[:,1], COMAPData_class[:,1], COMAPData_class[:,2]).gaussian_likelihood()
                
                # print(f'Likelihood: {Likelihood}')
                
                return Likelihood


def Independent_HI_Likelihood(params):
        
        OmegaHI_ext = np.round(params,5).item()

        
        # Priors for Bayesian Inference
        
        if not (0.00032 <=OmegaHI_ext<= 0.00068): #  1 sigma-priors (Rhee et al. 2018) --> z = 0.32
                return -np.inf
        
        # if not (0.00051 <=OmegaHI_ext<= 0.00103): #  1 sigma-priors (Rao et al. 2017) --> z = 0.44
                # return -np.inf

        else:
                lim_simulator_HI = ce.CodeExecutor("/path/to/HI_LIM_SIM_0.33",  # --> change path for different z 
                                        "LIM_simulator",
                                        "param.txt",
                                        "paths_HI.txt",
                                        str(OmegaHI_ext))
                lim_simulator_HI.run()
        

                powerspectrum_HI = ce.CodeExecutor("/path/to/New_power_spec_0.33_HI", # --> change path for different z 
                                        "PowerSpectrum",
                                        "input")
                powerspectrum_HI.run()
                
                # Fetches MeerKAT data 
                MeerData_Class = cmd.MeerKATPowerSpectrum(0.32).get_data()
                # MeerData_Class = cmd.MeerKATPowerSpectrum(0.44).get_data() --> 0.44

                # Fetches Theoretical Powerspectrum
                ModernPS_HI = cbm.TheoreticalPowerSpectrum('HI_0.42Mpc_RSD_z.0.33.bin.powspec').read_power_spectrum()
                ModernPS_HI = ModernPS_HI[1:8:2]

                # ModernPS_HI = cbm.TheoreticalPowerSpectrum('HI_0.42Mpc_RSD_z.0.44.bin.powspec').read_power_spectrum() --> for 0.44
                
                
                # Calculates Likelihood
                Likelihood = cbm.MultiVariateLogGaussianLikelihood(ModernPS_HI[:,1], MeerData_Class[:,1], MeerData_Class[:,2]).gaussian_likelihood()
                
                # print(f'Likelihood: {Likelihood}')
                
                return Likelihood

def JointLikelihood(params):
        
        # Multiple params 
        ext_alpha, ext_OmegaHI  = params

        # rounding off to significant digits
        ext_alpha, ext_OmegaHI = np.round(ext_alpha,2), np.round(ext_OmegaHI, 5)
   
        
        if not (0.85 <= ext_alpha <= 1.49 and 0.00032 <=ext_OmegaHI<= 0.00068):
                return -np.inf

        else:
                # Runs the LIM simulator and PowerSpectrum from the shell
                
                # Constant beta (Greve) -> Subjecct to change!
                Const_Beta = 2.0

                lim_simulator_CO = ce.CodeExecutor("/path/to/CO_LIM_SIM_3.0",
                                        "LIM_simulator",
                                        "param.txt",
                                        "paths_CO.txt",
                                        str(ext_alpha), str(Const_Beta))
                lim_simulator_CO.run()
        

                powerspectrum_CO = ce.CodeExecutor("/path/to/New_power_spec_3_CO",
                                        "PowerSpectrum",
                                        "input")
                powerspectrum_CO.run()
                
                
                lim_simulator_HI = ce.CodeExecutor("/path/to/HI_LIM_SIM_0.33",
                                        "LIM_simulator",
                                        "param.txt",
                                        "paths_HI.txt",
                                        str(ext_OmegaHI))
                lim_simulator_HI.run()
        

                powerspectrum_HI = ce.CodeExecutor("/path/to/New_power_spec_0.33_HI",
                                        "PowerSpectrum",
                                        "input")
                powerspectrum_HI.run()
                
                
                # Fetching COMAP and MeerKAT data

                COMAPData_class = cmd.COMAPPowerSpectrum(3.0).get_data()
                MeerData_Class = cmd.MeerKATPowerSpectrum(0.32).get_data()
                
                # Fetches Theoretical Powerspectrum for CO and HI

                ModernPS_CO = cbm.TheoreticalPowerSpectrum('CO_JointAnalysis_0.56MPc_z3.010.bin.powspec').read_power_spectrum()
                ModernPS_CO = ModernPS_CO[:-2]
                
                ModernPS_HI = cbm.TheoreticalPowerSpectrum('HI_JointAnalysis_0.42Mpc_z.0.33.bin.powspec').read_power_spectrum()
                ModernPS_HI = ModernPS_HI[1:8:2]
                
                
                #Calculates Likelihood
                Joint_Likelihood = cbm.JointLikelihood(
                                    Model_1=ModernPS_CO[:,1],Data_1=COMAPData_class[:,1],Error_1=COMAPData_class[:,2],
                                    Model_2=ModernPS_HI[:,1],Data_2=MeerData_Class[:,1],Error_2=MeerData_Class[:,2]
                                    ).JointLogLikelihood()
                
                # print(f'Likelihood: {Joint_Likelihood}')
                
                return np.array([Joint_Likelihood])


# Number of dimensions and walkers
ndim = 2
nwalkers = 4
nsteps = 6000

# backup using built in method (HDF5 files) 
filename = "JointAnalysis_COz03_HIz032.h5" 
# Hierarchical structure (1 HDF5 file, multiple chain run saves); change name for each new mcmc run
backend = emcee.backends.HDFBackend(filename,name="JointAnalysis_COz03_HIz032_withprior_for_alpha_Omega_keeping_const_beta_2.00")
backend.reset(nwalkers, ndim)



# Parameter space bounds (Initialization of walkers)
alpha_bounds = [0.85, 1.49]
OmegaHI_bounds = [0.00032, 0.00068]

# Initial positions for the walkers at key points in the parameter space
initial_positions = np.array([
    [alpha_bounds[0], OmegaHI_bounds[0]],  # Min alpha, Min Omega_HI
    [alpha_bounds[1], OmegaHI_bounds[1]],  # Max alpha, Max Omega_HI
    [alpha_bounds[0], OmegaHI_bounds[1]],  # Min alpha, Max Omega_HI
    [alpha_bounds[1], OmegaHI_bounds[0]]   # Max alpha, Min Omega_HI
])

print("Initial positions for walkers:", initial_positions)

# Initiating sampler : Change Likelihood function for Different MCMC
sampler = emcee.EnsembleSampler(nwalkers, ndim, JointLikelihood,backend=backend)

# Running sampler
sampler.run_mcmc(initial_positions, nsteps, progress=True)


print(f'time taken: {time.time() - start}')
print("Code run finished")