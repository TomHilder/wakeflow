from setup import run_setup
from grid import Grid
from phantom_interface import PhantomDump
from linear_perts import LinearPerts
from non_linear_perts import NonLinearPerts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import ticker, cm
from scipy.interpolate import RectBivariateSpline


def midplane_comparison(phantom_dump, analytics_file_location):

    # ======= DATA IMPORT ======= #

    # rename for convenience
    PD = phantom_dump

    # import analytics cartesian results
    A_vr = np.transpose(np.load(analytics_file_location + '/vr.npy'))
    A_vphi = np.transpose(np.load(analytics_file_location + '/vphi.npy'))
    A_rho = np.transpose(np.load(analytics_file_location + '/density.npy'))

    # scale density for analytics
    A_rho *= np.median(PD.rho_xy) / np.median(A_rho)

    # ======= INTERPOLATION TO SAME GRID AND RESIDUALS ======= #

    # get grid for analytics
    A_x = np.linspace(-500, 500, 1000)
    A_y = np.linspace(-500, 500, 1000)
    #A_X, A_Y = np.meshgrid(A_x, A_y)

    # interpolating functions for analytics
    interp_v_r = RectBivariateSpline(A_y, A_x, A_vr)
    interp_v_phi = RectBivariateSpline(A_y, A_x, A_vphi)
    interp_v_rho = RectBivariateSpline(A_y, A_x, A_rho)

    # new grid for analytics (same as phantom)
    new_x = np.linspace(-500, 500, 600)
    new_y = np.linspace(-500, 500, 600)
    new_X, new_Y = np.meshgrid(new_x, new_y)

    # evaluate interpolation on Grid object grid
    A_int_vr = interp_v_r.ev(new_Y, new_X) + 1E-8
    A_int_vphi = interp_v_phi.ev(new_Y, new_X)
    A_int_rho = interp_v_rho.ev(new_Y, new_X)

    # find residuals
    
    R_vr = PD.vr_xy - A_int_vr
    R_vphi = np.abs(A_int_vphi - PD.vphi_xy)
    R_rho = np.abs(A_int_rho - PD.rho_xy)
    
    """
    R_vr = np.abs((PD.vr_xy) - A_int_vr)
    R_vphi = np.abs((PD.vphi_xy) - A_int_vphi / PD.vphi_xy)
    R_rho = np.abs((PD.rho_xy) - A_int_rho / PD.rho_xy)
    """

    # ======= PLOTTING ======= #
    
    fig, ax = plt.subplots(3, 3, figsize=(10,8), sharex=True, sharey=True, constrained_layout=True)

    # V_R
    pcm_vr = ax[0,0].imshow(A_vr, cmap='RdBu', interpolation='none', extent=[-500,500,500,-500], vmin=-0.5, vmax=0.5)
    ax[0,1].imshow(PD.vr_xy, cmap='RdBu', interpolation='none', extent=[-500,500,500,-500], vmin=-0.5, vmax=0.5)
    ax[0,2].imshow(R_vr, cmap='RdBu', interpolation='none', extent=[-500,500,500,-500], vmin=-0.5, vmax=0.5)

    #V_PHI
    pcm_vphi = ax[1,0].imshow(-1*A_vphi, cmap='inferno', interpolation='none', extent=[-500,500,500,-500], vmin=0, vmax=8)
    ax[1,1].imshow(-1*PD.vphi_xy, cmap='inferno', interpolation='none', extent=[-500,500,500,-500], vmin=0, vmax=8)
    ax[1,2].imshow(R_vphi, cmap='inferno', interpolation='none', extent=[-500,500,500,-500], vmin=0, vmax=8)

    # RHO
    pcm_rho = ax[2,0].imshow(A_rho, cmap='inferno', interpolation='none', extent=[-500,500,500,-500], norm=LogNorm(vmin=1E-19, vmax=1E-14))
    ax[2,1].imshow(PD.rho_xy, cmap='inferno', interpolation='none', extent=[-500,500,500,-500], norm=LogNorm(vmin=1E-19, vmax=1E-14))
    ax[2,2].imshow(R_rho, cmap='inferno', interpolation='none', extent=[-500,500,500,-500], norm=LogNorm(vmin=1E-19, vmax=1E-14))

    # colourbars
    fig.colorbar(pcm_vr, ax=ax[0, :], location='right', label='Radial Velocity [km/s]')
    fig.colorbar(pcm_vphi, ax=ax[1, :], location='right', label='Azimutal Velocity [km/s]')
    fig.colorbar(pcm_rho, ax=ax[2, :], location='right', label='Density [g/cm$^3$]')

    # labels
    ax[0,0].set(ylabel='[au]', aspect='equal', adjustable='box')
    ax[1,0].set(ylabel='[au]', aspect='equal', adjustable='box')
    ax[2,0].set(ylabel='[au]', xlabel='[au]', aspect='equal', adjustable='box')
    ax[2,1].set(xlabel='[au]', aspect='equal', adjustable='box')
    ax[2,2].set(xlabel='[au]', aspect='equal', adjustable='box')

    # titles
    ax[0,0].set_title('Analytical Model')
    ax[0,1].set_title('Simulated Model')
    ax[0,2].set_title('Residuals')

    plt.savefig('midplane_comparisons.pdf')
    #plt.show()



# ========== TESTING ========== #


# RUN SETUP
params = run_setup()

# MAKE EMPTY GRID FOR PHANTOM MIDPLANE
g = Grid(params)
g.make_grid()
g.make_empty_disk()

# MAKE PHANTOM DUMP OBJECT
PD = PhantomDump(params, g)

midplane_comparison(PD, 'phantom_pixelmaps/HD163 analytic cartesian closer')
