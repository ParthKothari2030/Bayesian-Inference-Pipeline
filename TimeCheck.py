import Code_executor as ce
import time

starr = time.time()
lim_simulator = ce.CodeExecutor("/media/disk1/parth/MeerKAT_0.33_MCMC/HI_LIM_SIM_0.33",
                                    "LIM_simulator",
                                    "param.txt",
                                    "paths_HI.txt",
                                    str(0.00072)).run()
    
powerspectrum = ce.CodeExecutor("/media/disk1/parth/MeerKAT_0.33_MCMC/New_power_spec_0.33_HI",
                                    "PowerSpectrum",
                                    "input").run()

print(time.time() - starr)