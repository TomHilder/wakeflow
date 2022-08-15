"""
__init__.py

Written by Thomas Hilder
Last modified 11.08.2022

This file gives users access to classes/functions contained in Wakeflow, and will be modified as functionality increases.
"""

# wakeflow package versionss
__version__ = "1.0.4"

# give users access to the WakeflowModel class
from .wakeflow import WakeflowModel

# give users access to the Grid class
from .grid import Grid