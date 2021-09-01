#------------------------------#
# external library imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import ticker, cm

# wakeflow imports
from setup import run_setup
from grid import Grid
from phantom_interface import PhantomDump
from linear_perts import LinearPerts
from non_linear_perts import NonLinearPerts
#------------------------------#



#------------------------------#
# Setup
#------------------------------#

params = run_setup()

#------------------------------#



#------------------------------#
# Unperturbed disk
#------------------------------#

# make empty grid for unperturbed disk
grid_background = Grid(params)
grid_background.make_grid()

# fill grid with Keplerian, power law disk
grid_background.make_keplerian_disk()

#------------------------------#



#------------------------------#
# Linear Perturbations
#------------------------------#

# make empty grid for linear perturbations
grid_lin_perts = Grid(params)
grid_lin_perts.make_grid()
grid_lin_perts.make_empty_disk()

# extract linear perturbations from file
lin_perts = LinearPerts(params)
lin_perts.cut_box_annulus_segment()

# add the linear perturbations onto grid
grid_lin_perts.add_linear_perturbations(lin_perts, grid_background.rho)

#------------------------------#



#------------------------------#
# Non-Linear Perturbations
#------------------------------#

# make empty grid for non-linear perturbations
grid_nonlin_perts = Grid(params)
grid_nonlin_perts.make_grid()
grid_nonlin_perts.make_empty_disk()

# initialise non-linear perturbations
nonlin_perts = NonLinearPerts(params, grid_nonlin_perts)

# extract initial condition from the linear perturbations
nonlin_perts.extract_ICs(lin_perts)

# solve for non-linear perturbations
nonlin_perts.get_non_linear_perts()

# add non-linear perturbations to grid
grid_nonlin_perts.add_non_linear_perturbations(nonlin_perts, grid_background.rho)

#------------------------------#



#------------------------------#
# Merge results
#------------------------------#

grid_background.merge_grids(grid_lin_perts)
grid_background.merge_grids(grid_nonlin_perts)
#grid_background.show_disk2D(0)

grid_lin_perts.merge_grids(grid_nonlin_perts)
grid_lin_perts.show_disk2D(0)

#------------------------------#



#------------------------------#
# Write FITS file
#------------------------------#

grid_background.write_fits_file()

#------------------------------#