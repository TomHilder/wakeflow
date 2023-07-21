# utilities.py
# Written by Thomas Hilder

"""
Contains the functions used for smoothing over the edge of the linear box between regimes.
"""

import numpy                as np

# NOTE: contents are intended for internal use and should not be directly accessed by users

def sigmoid(x, x0, a):
    return (1 + np.exp(-1*(x - x0) / a))**-1

def sigmoid_negx(x, x0, a):
    return (1 + np.exp((x - x0) / a))**-1

def tanh_sigmoid(x, x0, a):
    return 0.5 + 0.5 * np.tanh((x - x0) / a)

def tanh_sigmoid_negx(x, x0, a):
    return 0.5 + 0.5 * np.tanh(-(x - x0) / a)

def sigmoid_smoothing(f1, f2, x, x0, a):
    return (1 - sigmoid(x, x0, a)) * f1 + sigmoid(x, x0, a) * f2 

def sigmoid_smoothing_3_component(f1, f2, f3, x, x0_1, x0_2, a):
    s1_func = tanh_sigmoid
    s2_func = tanh_sigmoid_negx
    # make sigmoids
    sigmoid_1 = s1_func(x, x0_1, a)
    sigmoid_2 = s2_func(x, x0_2, a)
    # perform weighting
    left = (1 - sigmoid_1) * f1
    middle = (sigmoid_1 + sigmoid_2 - 1) * f2
    right = (1 - sigmoid_2) * f3
    # return total
    return left + middle + right