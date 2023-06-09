# utilities.py
# Written by Thomas Hilder

"""
Contains the functions used for smoothing over the edge of the linear box between regimes.
"""

import numpy                as np

# NOTE: contents are intended for internal use and should not be directly accessed by users

def sigmoid(x, x0, a):
    return (1 + np.exp(-1*(x - x0) / a))**-1

def sigmoid_smoothing(f1, f2, x, x0, a):
    return (1 - sigmoid(x, x0, a)) * f1 + sigmoid(x, x0, a) * f2 