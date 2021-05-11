from setup import run_setup
from grid import Grid
from linear_perts import LinearPerts
from non_linear_perts import NonLinearPerts
import numpy as np
import matplotlib.pyplot as plt

# RUN SETUP
params = run_setup()

# GENERATE GRID AND BACKGROUND KEPLERIAN DISK
g1 = Grid(params)
g1.make_grid()
g1.make_keplerian_disk()
print(g1.info)

# EXTRACT LINEAR PERTURBATIONS NEAR PLANET
lin = LinearPerts(params)
lin.cut_box_square()

# CALCULATE NON-LINEAR PERTURBATIONS
nonlin = NonLinearPerts(params, g1)
nonlin.extract_burgers_ICs_sq(lin)
nonlin.get_non_linear_perts()

# PLOT NON-LINEAR PERTS
plt.imshow(nonlin.rho)
plt.show()