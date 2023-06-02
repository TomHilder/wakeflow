# height.py
# Written by Thomas Hilder

"""
Contains the HeightFunctions class, which allows users to generate functions that return the disk height as a function of radius, using the parameterisation of their choice.
"""

import numpy as np

from typing import Callable

class HeightFunctions():
    """
    Class containing methods to generate functions that allow users to get disk height as a function of radius (most-likely used for the height of the emitting layer for a certain line). 
    """
    
    @staticmethod
    def powerlaw(z_ref: float, r_ref: float, power_index: float) -> Callable:
        """Power-law height profile height = z_ref * (r / r_ref)**power_index.
        
        Parameters
        ----------
        z_ref : disk height at r_ref
        r_ref : reference radius
        power_index : power-law index for the height profile
        
        returns : Callable function that will return the height value at a provided radius, for the power-law profile.
        """
        # return power-law function
        return lambda r : z_ref * np.power(r / r_ref, power_index)
    
    @staticmethod
    def powerlaw_tapered(z_ref: float, r_ref: float, power_index: float, r_taper: float, power_index_taper: float) -> Callable:
        """Power-law height profile height = z_red * (r / r_ref)**power_intex * exp(- (r / r_taper)**power_index_taper).
        
        Parameters
        ----------
        z_ref : disk height at r_ref
        r_ref : reference radius
        power_index : power-law index for the height profile
        r_taper : taper radius
        power_index_taper : power-law index for the exponential taper
        
        returns : Callable function that will return the height value at a provided radius, for the tapered power-law profile.
        """
        # get regular power-law
        powerlaw = HeightFunctions.powerlaw(z_ref, r_ref, power_index)
        # return modification with taper
        return lambda r : powerlaw(r) * np.exp(
            -1 * np.power(
                r / r_taper, power_index_taper
                )
        )
        
    @staticmethod
    def interpolate_tabulated() -> Callable:
        """Height function defined by interpolating over user-specified values. NOT IMPLEMENTED.
        
        Parameters
        ----------
        """
        # TODO: implement method
        return NotImplementedError