# spirals.py
# Written by Thomas Hilder

"""
Contains the Spiral class that allows uses to manipulate, project and rotate different spiral shapes easily.
"""

import numpy as np

from typing     import Tuple, Callable
from .rotations import rotation_matrix

class Spiral():
    
    # coordinates of spiral
    _rad = None
    _phi = None
    _z   = None
    
    def __init__(self, radii: np.ndarray, phi_function: Callable) -> None:
        # TODO: add some input checks
        self._rad = radii
        self._phi = phi_function(radii)
    
    @property
    def rad_phi(self) -> Tuple[np.ndarray, np.ndarray]:
        # return the coordinates
        return self._rad, self._phi
    
    @property
    def xy(self) -> Tuple[np.ndarray, np.ndarray]:
        # get polar coordinates
        rad, phi = self.rad_phi
        # convert to Cartesian
        x = rad * np.cos(phi)
        y = rad * np.sin(phi)
        # return the results
        return x, y
    
    @property
    def height(self) -> np.ndarray:
        # if height not yet calculated raise an exception
        if self._z == None:
            raise Exception("Height coordinate undefined. Please use get_height first.")
        # if height coords have been calculated then return it
        return self._z
    
    def get_height(self, height_func: Callable) -> None:
        # get radial coordinates
        rad, _ = self.rad_phi
        # set heights using height function on radial coordinates
        self._z = height_func(rad)
    
    def perform_rotation(self, position_angle: float, inclination: float, azimuth_offset: float) -> None:
        # get Cartesian coordinates
        x, y = self.xy
        z    = self.height
        # define rotation matrices
        rot_azi = rotation_matrix(azimuth_offset, "z")
        rot_inc = rotation_matrix(   inclination, "x")
        rot_pos = rotation_matrix(position_angle, "z")
        # perform rotations with a loop
        for i in range(x.shape[0]):
            # rotate around z-axis to induce azimuthal offset
            x[i], y[i], z[i] = np.dot(rot_azi, [x[i], y[i], z[i]])
            # rotate around x-axis of sky-plane to match inclination
            x[i], y[i], z[i] = np.dot(rot_inc, [x[i], y[i], z[i]])
            # rotate around z-axis of sky-plane to match position angle
            x[i], y[i], z[i] = np.dot(rot_pos, [x[i], y[i], z[i]])
        # save results to private attributes
        self._x = x
        self._y = y
        self._z = z