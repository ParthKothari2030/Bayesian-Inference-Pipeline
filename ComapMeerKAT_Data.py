import numpy as np

"""
Module name: COMAP_MeerKAT_Data.
Function: Fetching Observed Powerspectrum of COMAP and MeerKAT.
Uses: Used in MCMC likelihood calculation -> ```P_data ```.
"""


class MeerKATPowerSpectrum:
    """
    A class to store and retrieve MeerKAT power spectrum data for different redshifts.
    """

    def __init__(self, redshift):
        
        """
        Initializes the MeerKATPowerSpectrum class with the requested redshift.

        Args:
            redshift (float): The redshift value for which power spectrum data is required.
        """
        
        if not isinstance(redshift,float):
            raise ValueError("Redshift must be a float type value.")  
        
        
        self.redshift = redshift  # Store redshift as an instance variable

    def get_data(self):
        """
        Returns the MeerKAT power spectrum data for the given redshift.

        Paul et al. 2023 -> Paper for data

        Returns:
            np.ndarray: 2D array with columns [k {Mpc}, P(k) {mk^2 Mpc^3}, Uncertainty {mk^2 Mpc^3}].
        
        Raises:
            ValueError: If the redshift is not supported.
        """
        
        # Scalable data frame if more redshifts are to be added.

        data = {
            0.32: np.array([
                [1.25, 2.26, 0.57],
                [2.04, 0.96, 0.19],
                [3.30, 0.19, 0.09],
                [5.19, 0.33, 0.09]
            ]),
            0.44: np.array([
                [1.01, 4.34, 0.85],
                [1.68, 1.51, 0.27],
                [2.81, 0.96, 0.13],
                [4.31, 0.56, 0.15]
            ])
        }

        # Checks avaibility of redshift PS in data base. 
        if self.redshift in data:
            return data[self.redshift]
        else:
            raise ValueError(f"""Redshift z = {self.redshift} is not supported.
            Values supported as of now are: 0.32, 0.44""")


class COMAPPowerSpectrum:
    """
    A class to store and retrieve COMAP power spectrum data for different redshifts.
    
    Stutzer et al. 2024 -> Paper for data.

    """
    
    def __init__(self, redshift=3.0):
        """
        Initializes the COMAPPowerSpectrum class with the requested redshift.

        Args:
            redshift (float, optional): The redshift value for which power spectrum data is required. Defaults to 3.0.
        """
        self.redshift = redshift

        if not isinstance(redshift,float):
            raise ValueError("Redshift must be a float type value.")  
    
    def get_data(self):
        """
        Returns the COMAP power spectrum data for the given redshift.

        Returns:
            np.ndarray: 2D array with columns [k {Mpc^{-1}}, P(k) {mk^2 Mpc^3}, Uncertainty {mk^2 Mpc^3}].
        
        Raises:
            ValueError: If the redshift is not supported.
        """

        COMAP_Data = {
        3.0: np.array([
            [0.10,  3600.00, 18200.00, 1820.00],  # k, P(k), Uncertainty, k_sigma
            [0.15, 19333.33,  7266.67, 1090.00],
            [0.21,  2809.52,  6047.62, 1270.00],
            [0.30,  3966.67,  4200.00, 1260.00]
            ])
        }

        if self.redshift in COMAP_Data:
            return COMAP_Data[self.redshift]
        else:
            raise ValueError(f"""Redshift z = {self.redshift} is not supported.
            Values supported as of now are: 3.0 """)
