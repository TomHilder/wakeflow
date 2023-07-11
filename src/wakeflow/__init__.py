# __init__.py
# Written by Thomas Hilder

"""
Wakeflow allows users to generate and manipulate semi-analytic models of planet wakes.
"""

# wakeflow package versions
__version__ = "1.4.0"

# give users access to the WakeflowModel class
from .wakeflow import WakeflowModel

# give users access to the velocity perturbations only interface
from .v_perts_only import _VelocityPerturbations

# give users access to the Grid class, realistically they shouldn't need it but it 
# may be desirable by advanced users or the developers
from .grid import _Grid