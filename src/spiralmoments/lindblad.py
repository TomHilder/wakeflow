# lindblad.py
# Written by Thomas Hilder

"""
Contains a function to calculate the Lindblad resonance locations
"""

# TODO: add a function that calculates the resonant positions for an arbitrary Omega function (will involve root finding)

import numpy as np

from typing import Tuple

def lindblad_resonances(m: int, r_0: float) -> Tuple[float, float]:
    # convert m to a float to make sure we don't do integer division
    m_float = float(m)
    # inner resonance
    r_inner = np.power(
        1 - np.reciprocal(m_float),
        2. / 3.
    )
    # outer resonance
    r_outer = np.power(
        1 + np.reciprocal(m_float),
        2. / 3.
    )
    # return results scaled by r_0
    return r_inner * r_0, r_outer * r_0