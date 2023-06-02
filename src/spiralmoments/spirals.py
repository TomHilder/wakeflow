# spirals.py
# Written by Thomas Hilder

"""
Contains the Spiral class that allows uses to manipulate, project and rotate different spiral shapes easily.
"""

import numpy as np

from typing     import Tuple, Callable
from .rotations import rotation_matrix

class Spiral():
    """
    Class representing spiral shape, to be instantiated by user. Allows for easy 3D projection and rotation to find the line-of-sight shape for a spiral on the surface of a disk.
    """
    
    # polar disk coordinates of spiral
    _rad = None
    _phi = None
    _z_disk   = None
    
    # cartesian sky coordinates of spiral
    _x = None
    _y = None
    _z = None
    
    def __init__(self, radii: np.ndarray, phi_function: Callable) -> None:
        """Instantiate a Spiral object, defined by a function phi(r) and given r values, where phi and r are the azimuthal and radial disk coordinates respectively. 
        
        Parameters:
        -----------
        radii : array containing disk radii at which spiral shape is to be calculated
        phi_function : Callable function that takes disk radius as an argument and returns the respective azimuthal coordinate
        """
        # TODO: add some input checks
        self._rad = radii
        self._phi = phi_function(radii)
    
    @property
    def rad_phi(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get radial and azimuthal disk coordinates of the spiral.
        
        Parameters:
        -----------
        
        returns : Tuple containing numpy.ndarrays, where the first is the radial disk coordinates of the spiral, and the second is the azimuthal.
        """
        # return the coordinates
        return self._rad, self._phi
    
    @property
    def xy(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get the cartesian disk coordinates of the spiral. 
        
        Parameters:
        -----------
        
        returns: Tuple containing numpy.ndarrays, where the first is the disk x coordinates of the spiral, and the second is the y coordinates.
        """
        # get polar coordinates
        rad, phi = self.rad_phi
        # convert to Cartesian
        x = rad * np.cos(phi)
        y = rad * np.sin(phi)
        # return the results
        return x, y
    
    @property
    def height(self) -> np.ndarray:
        """Get height of the spiral on the disk surface. Requires that get_height has been called previously.
        
        Parameters:
        -----------
        
        returns : numpy.ndarray containing the z-coordinates of the spiral in the disk frame.
        """
        # if height not yet calculated raise an exception
        if self._z_disk == None:
            raise Exception("Height coordinate undefined. Please use get_height first.")
        # if height coords have been calculated then return it
        return self._z_disk
    
    def get_height(self, height_func: Callable) -> None:
        """Calculate the height of the spiral on the disk-surface using function that gives the disk height as a function of radius.
        
        Parameters:
        -----------
        height_func : Callable function that takes disk radius as an argument and returns the respective disk height
        
        returns : None
        """
        # get radial coordinates
        rad, _ = self.rad_phi
        # set heights using height function on radial coordinates
        self._z_disk = height_func(rad)
    
    def perform_rotation(self, position_angle: float, inclination: float, azimuth_offset: float) -> None:
        """Rotate the spiral shape to match the position on the sky, given a disk position angle and inclination. Also rotates the disk around its own z-axis by azimuth_offset (to match a planet position or similar).
        
        Parameters:
        -----------
        position_angle : disk position angle in radians, defined counter-clockwise from the positive x-axis direction
        inclincation : disk inclination in radians
        azimuth_offset : angle through which to rotate the spiral around the disk z-axis.
        
        returns : None
        """
        
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
    
    @property
    def sky_coords(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get the sky coordinates of the spiral. Requires that perform_rotation has been called.
        
        Parameters:
        -----------
        
        returns : Sky coordinates of spiral.
        """
        if self._x is not None:
            return self._x, self._y, self._z
        else:
            raise Exception("Sky coordinates not yet calculated. Please call perform_rotation first.")