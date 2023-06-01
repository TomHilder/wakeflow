# phi.py
# Written by Thomas Hilder

"""
Contains the PhiFunctions class, which allows users to generate functions that return the azimuthal position of a spiral as a function of radius, using the parameterisation of their choice. Includes options to use physical spirals from the Lin-Shu dispersion relation, including the planet wake shape (Ogilvie and Lubow 2002), as well as parameterised spirals such as logarithmic or power-law spirals.
"""

import numpy as np

from scipy.integrate import quad
from typing import Callable
from .lindblad import lindblad_resonances

class PhiFunctions():
    """
    Class containing methods to generate functions that allow users to get spiral shape functions \phi(r) using either simple or physically-motivated parameterisations. In the latter case, \phi(r) can be generated from arbitrary rotation and sound speed curves.
    """
    
    @staticmethod
    def log(b: float, c: float) -> Callable:
        """Logarithmic spiral, defined by phi = b * log(r) + c. It has constant pitch angle arctan(1 / b).
        
        Parameters
        ----------
        b : spiral scaling factor (pitch angle = arctan(1 / b))
        c : constant offset in phi
        
        returns : Callable function that will return the phi value at a provided radius, given the defined spiral shape.
        """
        return lambda r : b * np.log(r) + c
    
    @staticmethod
    def powerlaw(a: float, c: float, power_index: float) -> Callable:
        """Power-law spiral, defined by phi = a * r**power_index + c. The pitch-angle evolves as a power-law. TODO: add form for pitch-angle
        
        Parameters
        ----------
        a : spiral scaling factor
        c : constant offset in phi
        
        returns : Callable function that will return the phi value at a provided radius, given the defined spiral shape.
        """
        return lambda r : a * np.power(r, power_index) + c
    
    @staticmethod
    def powerlaw_plus_log(a: float, b: float, c: float, power_index: float) -> Callable:
        """Composite spiral with both a power-law and a logarithmic component. Defined by phi = a * r**power_index + b * log(r) + c. Corresponds to
        a spiral that evolves with a pitch-angle with both a constant and a power-law component. TODO: add form for pitch-angle
        
        Parameters
        ----------
        a : power-law spiral scaling factor
        b : log spiral scaling factor (pitch angle approaches arctan(1 / b) as r approaches zero)
        c : constant offset in phi
        
        returns : Callable function that will return the phi value at a provided radius, given the defined spiral shape.
        """
        return lambda r : a * np.power(r, power_index) + b * np.log(r) + c
    
    @staticmethod
    def planet_wake_powerlaw_disk(phi_planet: float, r_planet: float, hr_planet: float, index_cs: float, rotation_sign: int) -> Callable:
        """Ogilvie and Lubow (2002) analytical planet wake shape for a power-law disk with Keplerian rotation.
        
        Parameters
        ----------
        phi_planet : azimuthal location of the planet in radians
        r_planet : orbital radius of the planet
        hr_planet : disk aspect ratio H/r at the planet orbital radius r=r_planet. H is the pressure scale height H = cs / Omega.
        index_cs : power law index for the radial power-law parameterisation of the sound speed in the disk cs = cs_planet * (r / r_planet)**-index_cs
        rotation_sign : sign of the rotation. +1 for clockwise and -1 for anticlockwise.
        
        returns : Callable function that will return the phi value at a provided radius, given the defined planet wake shape.
        """
        # scale r by the planet radius
        rr = lambda r : r / r_planet
        # calculate the position of the wake relative to the planet
        phi_wake = lambda r : np.reciprocal(float(hr_planet)) * (
            np.power(rr(r), index_cs - 0.5) / (index_cs - 0.5) -
            np.power(rr(r), index_cs + 1.0) / (index_cs + 1.0) -
            3 / (
                (2 * index_cs - 1.0) * (index_cs + 1.0)
            )
        )
        # return function for the absolute position of the wake
        return lambda r : phi_planet + rotation_sign * np.sign(r - r_planet) * phi_wake(r)
    
    @staticmethod
    def planet_wake_numerical(phi_planet: float, r_planet: float, Omega_func: float, cs_func: float, rotation_sign: int) -> Callable:
        """Planet wake form given in Rafikov (2002), defined in terms of the rotation curve Omega(r) and the sounds speed cs(r), calculated numerically.
        
        Parameters
        ----------
        phi_planet : azimuthal location of the planet in radians
        r_planet : orbital radius of the planet
        Omega_func : function that takes only r as an argument and returns the rotation Omega at r
        cs_func : function that takes only r as an argument and returns the sound speed cs at r
        rotation_sign : sign of the rotation. +1 for clockwise and -1 for anticlockwise.
        
        returns : Callable function that will numerically evaluate and return the phi value at a provided radius, given the defined planet wake shape.
        """
        # get a function for the integrand
        integrand = lambda r : (Omega_func(r) - Omega_func(r_planet)) / cs_func(r)
        # get a function that returns the result of the integral from r_planet to r
        integral = lambda r : quad(integrand, r_planet, r)[0]
        # return function for absolute position of the wake
        return np.vectorize(
            lambda r : phi_planet + rotation_sign * np.sign(r - r_planet) * integral(r)
            )
    
    @staticmethod
    def m_mode_spiral_numerical(phi_0: float, r_0: float, m: int, Omega_func: float, cs_func: float, rotation_sign: int) -> Callable:
        """Individual spiral mode with azimuthal wave number m, defined in terms of the rotation curve Omega(r) and the sounds speed cs(r). Keep in mind that there are evanescent regions interior to the locations of the Lindblad radii. TODO: add a useful reference here.
        
        Parameters
        ----------
        phi_0 : azimuthal launching location for the spiral wave
        r_0 : radial launching location for the spiral wave
        m : azimuthal wave number of the spiral wave
        Omega_func : function that takes only r as an argument and returns the rotation Omega at r
        cs_func : function that takes only r as an argument and returns the sound speed cs at r
        rotation_sign : sign of the rotation. +1 for clockwise and -1 for anticlockwise.
        
        returns : Callable function that will numerically evaluate and return the phi value at a provided radius, given the origin and mode of the wave.
        """
        # get a function to calculate the integrand (we use nan_to_num to ignore wave propagation in evanescent region)
        integrand = lambda r : (
            np.reciprocal(cs_func(r))
        ) * (
            np.nan_to_num(
                np.sqrt((Omega_func(r) - Omega_func(r_0))**2 - (Omega_func(r) / m)**2)
            )
        )
        # get a function that returns the result of the integral from r_0 to r
        integral = lambda r : quad(integrand, r_0, r)[0]
        # return function for absolute position of the spiral
        return np.vectorize(
            lambda r : phi_0 - rotation_sign * integral(r)
            )