# rotations.py
# Written by Thomas Hilder

"""
Contains the rotation_matrix function used to perform 3D rotations of vector fields.
"""

import numpy as np

def rotation_matrix(angle: float, axis='x') -> np.ndarray:
    """Function for rotation matrices
    
    Parameters
    ----------
    angle : angle (in radians) around which vector field should be rotated
    axis : axis around which vector field should be rotated, either 'x', 'y', or 'z'
    
    returns : numpy.ndarray with shape (3,3) representing the rotation matrix.
    """
    
    # get phi in [0, 2pi)
    angle = np.mod(angle, 2*np.pi)
    
    # get corresponding matrix for chosen axis
    if axis == "x":
        arr = [
            [1,             0,              0],
            [0, np.cos(angle), -np.sin(angle)],
            [0, np.sin(angle),  np.cos(angle)]
        ]
    elif axis == "y":
        arr = [
            [ np.cos(angle), 0, np.sin(angle)],
            [ 0,             1,             0],
            [-np.sin(angle), 0, np.cos(angle)]
        ]
    elif axis == "z":
        arr = [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle),  np.cos(angle), 0],
            [0,                          0, 1]
        ]
    # if bogus axis choice throw error
    else:
        raise ValueError("axis must be x, y or z")
    # return array
    return np.array(arr)
    