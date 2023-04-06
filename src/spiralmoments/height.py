# height.py
# Written by Thomas Hilder

"""
Contains TODO: add contents description
"""

import numpy as np

from typing import Callable

class HeightFunctions():
    
    @staticmethod
    def powerlaw(z_ref: float, r_ref: float, power_index: float) -> Callable:
        # return power-law function
        return lambda r : z_ref * np.power(r / r_ref, power_index)
    
    @staticmethod
    def powerlaw_tapered(z_ref: float, r_ref: float, power_index: float, r_taper: float, power_index_taper: float):
        # get regular power-law
        powerlaw = HeightFunctions.powerlaw(z_ref, r_ref, power_index)
        # return modification with taper
        return lambda r : powerlaw(r) * np.exp(
            -1 * np.power(
                r / r_taper, power_index_taper
                )
        )
        
    @staticmethod
    def tabulated() -> Callable:
        # TODO: implement method
        pass