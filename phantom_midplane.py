from phantom_interface import PhantomDump
from setup import run_setup
from grid import Grid
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import ticker, cm


# RUN SETUP
params = run_setup()

# MAKE EMPTY GRID FOR PHANTOM MIDPLANE
g = Grid(params)
g.make_grid()
g.make_empty_disk()

# MAKE PHANTOM DUMP OBJECT
PD = PhantomDump(params, g)
PD.get_polar_grid()

# ADD PHANTOM MID-PLANE EXTRAPOLATED TO 3D TO g
g.add_phantom_midplane(PD)
print(g.info)

# SAVE TO FITS FILE FOR MCFOST
g.write_fits_file()

