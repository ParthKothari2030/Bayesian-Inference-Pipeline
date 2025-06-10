import numpy as np
import os

"""
Module name: CosmologicalBayesianMethods.
Function: Contains function for Powerspectrum reading and calculation.
Uses: Likelihood calculation for MCMC using Powerspectrum data.
"""


class TheoreticalPowerSpectrum:
    def __init__(self, filepath):
        """
        Initialize the TheoreticalPowerSpectrum class with a file path.

        Args:
            filepath (str): Path to the power spectrum file.
        """
        self.filepath = filepath  # Store the file path as an instance variable

    def read_power_spectrum(self, remove_file=False):
        """ 
        Reads a power spectrum file and optionally removes it after reading.

        Args:
            remove_file (bool): If True, deletes the file after reading.

        Returns:
            np.ndarray or None: The power spectrum data if successful, else None.
        
        Raises:
            FileNotFoundError: If the file does not exist.
            IOError: If the file cannot be read.
        """
        try:
            PowerSpectrum_model = np.genfromtxt(self.filepath)

            if remove_file:
                os.remove(self.filepath)  # Remove the file only if remove_file=True
                
            # # raise an exeception here when the powerspectrum is zero before returning
            # if PowerSpectrum_model[:,1].all() == 0:
            #     raise ValueError('Something wrong')

            return PowerSpectrum_model

        
        except FileNotFoundError as F:
            print(f"Error: The file '{self.filepath}' was not found.")  
            return F

        except IOError:
            print(f"Error: Unable to read the file '{self.filepath}'.")
            return None

        except Exception as e:
            print(f"Unexpected error: {e}")
            return None
        

class MultiVariateLogGaussianLikelihood:
    """
    Computes the log-likelihood assuming a multivariate Gaussian distribution.

    This class calculates the log-likelihood for a dataset given model predictions,
    observed data, and associated errors under the assumption of independent 
    Gaussian errors.
    
    Attributes:
        Model (np.ndarray): A NumPy array containing model predictions.
        Data (np.ndarray): A NumPy array containing observed data.
        Error (np.ndarray): A NumPy array containing uncertainties in the data.
    """

    def __init__(self, Model, Data, Error):
        """
        Initializes the likelihood calculation with Model, Data, and Error.

        Args:
            Model (np.ndarray): Model predictions.
            Data (np.ndarray): Observed data.
            Error (np.ndarray): Measurement uncertainties (standard deviation for each data point).

        Raises:
            TypeError: If inputs are not NumPy arrays.
            ValueError: If Model, Data, and Error do not have the same shape.
        """
        if not all(isinstance(arr, np.ndarray) for arr in (Model, Data, Error)):
            raise TypeError("All inputs (Model, Data, Error) must be NumPy arrays.")
        
        if not (Model.shape == Data.shape == Error.shape):
            raise ValueError("Model, Data, and Error must have the same shape.")     
        
        self.Model = Model
        self.Data = Data
        self.Error = Error
        self.shape = len(Model)

    def gaussian_likelihood(self):
        """
        Computes the log-likelihood assuming a multivariate Gaussian distribution.

        Returns:
            float: The computed log-likelihood value.

        Notes:
            - Assumes errors are independent and represented by a diagonal covariance matrix.
        """
        Model_minus_Data = self.Model - self.Data
        
        SigmaSquared = np.diag(self.Error**2)
        SigmaSquared_inv = np.linalg.inv(SigmaSquared)

        # Compute log-likelihood
        log_likelihood = -np.dot(Model_minus_Data, np.dot (SigmaSquared_inv, Model_minus_Data.reshape(self.shape,1)))
        
        return log_likelihood
    

class JointLikelihood:
    """
    Computes the Joint Log-Likelihood of two surveys (or more in future development).

    This class utilizes the `MultiVariateLogGaussianLikelihood` class to compute
    individual likelihoods and sums them to obtain the joint likelihood.

    Attributes:
        likelihood1 (float): Log-likelihood value for the first dataset.
        likelihood2 (float): Log-likelihood value for the second dataset.
    
    """
    def __init__(self, Model_1, Data_1, Error_1, Model_2, Data_2, Error_2):
        

        """
        Args:
            Model_* (np.ndarray): Model predictions.
            Data_* (np.ndarray): Observed data.
            Error_* (np.ndarray): Measurement uncertainties (standard deviation for each data point).
        
        Notes:
            - This method calls `MultiVariateLogGaussianLikelihood` to compute individual likelihoods.
        """
        self.likelihood1 = MultiVariateLogGaussianLikelihood(Model_1, Data_1, Error_1).gaussian_likelihood()
        self.likelihood2 = MultiVariateLogGaussianLikelihood(Model_2, Data_2, Error_2).gaussian_likelihood()
    
    def JointLogLikelihood(self):
        
        """
        Computes the Joint Log-Likelihood by summing the individual likelihoods.

        Returns:
            float: The computed log-likelihood value
        """
        return self.likelihood1.item() + self.likelihood2.item()

    


