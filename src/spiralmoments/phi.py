# phi.py
# Written by Thomas Hilder

"""
Contains TODO: add contents description
"""

import numpy as np

from typing import Callable

class PhiFunctions():
    
    @staticmethod
    def log(b, c):
        return lambda r : b * np.log(r) + c
    
    @staticmethod
    def powerlaw(a, c, power_index):
        return lambda r : a * np.power(r, power_index) + c
    
    @staticmethod
    def powerlaw_plus_log(a, b, c, power_index):
        return lambda r : a * np.power(r, power_index) + b * np.log(r) + c
    
    @staticmethod
    def planet_wake_powerlaw_disk():
        pass
    
    @staticmethod
    def planet_wake_numerical():
        pass
    
    @staticmethod
    def m_mode_spiral_numerical():
        pass