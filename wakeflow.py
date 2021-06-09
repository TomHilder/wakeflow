from setup import run_setup
from grid import Grid
from phantom_interface import PhantomDump
from linear_perts import LinearPerts
from non_linear_perts import NonLinearPerts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import ticker, cm

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

# MAKE EMPTY GRID FOR NON-LINEAR PERTURBATIONS
g2 = Grid(params)
g2.make_grid()
g2.make_empty_disk()

# ADD NON-LINEAR PERTURBATIONS TO g2
g2.add_non_linear_perturbations(nonlin)
print(g2.info)

# MERGE KEPLERIAN DISK AND NON-LINEAR PERTURBATIONS
g1.merge_grids(g2)
print(g1.info)

# CHECK MIDPLANES
print("Displaying midplane to check solution")
g1.show_disk2D(0)

"""
# MAKE EMPTY GRID FOR PHANTOM MIDPLANE
gp = Grid(params)
gp.make_grid()
gp.make_empty_disk()

# MAKE PHANTOM DUMP OBJECT
PD = PhantomDump(params, gp)
PD.get_polar_grid()

# ADD PHANTOM MID-PLANE EXTRAPOLATED TO 3D TO g
gp.add_phantom_midplane(PD)
print(gp.info)

g1.show_disk2D(0, save=True, name='analytics')
gp.show_disk2D(0, save=True, name='phantom')
"""

"""
g1.merge_phantom_densities(g)
"""

"""
# PLOT TOTAL DISK V_R

#print(g1.v_phi[0,0,40:])
plt.imshow(np.log10(g1.rho[:,0,:]))
plt.colorbar()
plt.show()

plt.scatter(g1.R[:,40,:].flatten(), np.log10(g1.rho[:, 0, :]).flatten())
plt.show()
"""

"""
plt.show()
plt.imshow(g1.v_r[:,9,:])
plt.colorbar()
plt.show()
plt.imshow(g1.v_r[:,19,:])
plt.colorbar()
plt.show()
plt.imshow(g1.v_r[:,29,:])
plt.colorbar()
plt.show()
plt.imshow(g1.v_r[:,39,:])
plt.colorbar()
plt.show()
plt.imshow(g1.v_r[:,49,:])
plt.colorbar()
plt.show()
"""

"""
_, ax = plt.subplots(subplot_kw=dict(projection='polar'))
myplot = ax.contourf(g1.PHI[:,0,0], g1.R[0,0,:], g1.v_r[:,0,:].transpose(), levels=300, cmap='RdBu')
plt.colorbar(myplot)
plt.show()

_, ax = plt.subplots(subplot_kw=dict(projection='polar'))
myplot = ax.contourf(g1.PHI[:,0,0], g1.R[0,0,:], g1.v_phi[:,0,:].transpose(), levels=300, cmap='RdBu')
plt.colorbar(myplot)
plt.show()
"""


# SAVE TO FITS FILE FOR MCFOST
g1.write_fits_file()