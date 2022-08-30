# linear_perts.py
# Written by Thomas Hilder and Francesco Bollati

"""
Contains the LinearPerts class responsible for handling the linear regime of the models.
"""

# type hinting without circular imports
from typing                 import TYPE_CHECKING

import sys, pkg_resources, tarfile
import numpy                    as np
import matplotlib.pyplot        as plt
from scipy.interpolate      import RectBivariateSpline

if TYPE_CHECKING:
    from .model_setup       import _Parameters

# NOTE: contents are intended for internal use and should not be directly accessed by users

# class for storing the results from the linear regime
class _LinearPerts():
    """
    Class to store the results from the linear regime nearby the planet.
    """

    # intitalise object given parameters. locate linear results data and read in
    def __init__(self, parameters : '_Parameters') -> None:

        # grab parameters object
        self.p = parameters

        # get location of linear perturbations data files
        pert_loc = pkg_resources.resource_filename('wakeflow', 'data/linear_perturbations.npy')
        mesh_loc = pkg_resources.resource_filename('wakeflow', 'data/linear_perturbations_mesh.npy')

        # read perturbations from files
        try:
            perts = np.load(pert_loc)
            mesh  = np.load(mesh_loc)

        # in the case files have not been extracted from tarballs yet, extract them
        except FileNotFoundError:

            # open tarballs
            pert_tar = tarfile.open(f"{pert_loc}.tar.gz")
            mesh_tar = tarfile.open(f"{mesh_loc}.tar.gz")

            # extract npy files
            loc = pkg_resources.resource_filename('wakeflow', 'data')
            pert_tar.extractall(loc)
            mesh_tar.extractall(loc)

            # close tarballs
            pert_tar.close()
            mesh_tar.close()

            # read perturbations from files
            perts = np.load(pert_loc)
            mesh  = np.load(mesh_loc)

        # get perturbation arrays
        self.pert_v_r   = perts[0]
        self.pert_v_phi = perts[1]
        self.pert_rho   = perts[2]

        # grid
        self.X = mesh[0]
        self.Y = mesh[1]

        # linear perturbations read in grid
        x = self.X[0,:]
        y = self.Y[:,0]

        # define constants for linear perts
        self.x_box = 2 * self.p.scale_box
        
        # cut square box grid in linear regime
        self.x_cut = x[np.argmin(x < -self.x_box) : np.argmin(x < self.x_box) + 1]
        self.y_cut = y[np.argmin(y < -self.x_box) : np.argmin(y < self.x_box) + 1]

        # test plot
        if False:
            plt.figure()
            plt.contourf(x, y, self.pert_rho, levels=100, cmap='RdBu')
            plt.xlim(self.x_cut[0],self.x_cut[-1])
            plt.ylim(self.y_cut[0],self.y_cut[-1])
            plt.colorbar()
            plt.show()

    # old method of extracting linear perturbations 
    def _cut_box_square(self) -> None:
        """Deprecated. This is only used in the case where you want to plot the linear solution in t, eta space
        """

        # box size (in units of Hill radius), with default scale_box = 1. (note for conversions that self.p.l = 1 Hill radius in cgs)
        box_size = 2*self.p.scale_box
        artificial_y_scale = 6

        # linear perturbations read in grid
        x = self.X[0,:]
        y = self.Y[:,0]

        # cut square box grid in linear regime
        x_cut = x[np.argmin(x < -box_size)                    : np.argmin(x < box_size)                    + 1]
        y_cut = y[np.argmin(y < -artificial_y_scale*box_size) : np.argmin(y < artificial_y_scale*box_size) + 1]

        # find cut indicies 
        x_cut_i1 = np.argmin(x < -box_size)
        x_cut_i2 = np.argmin(x <  box_size) + 1
        y_cut_i1 = np.argmin(y < -artificial_y_scale * box_size)
        y_cut_i2 = np.argmin(y <  artificial_y_scale * box_size) + 1

        # cut perturbation arrays
        cut_v_r   = self.pert_v_r  [y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]
        cut_v_phi = self.pert_v_phi[y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]
        cut_rho   = self.pert_rho  [y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]

        self.cut_rho   = cut_rho
        self.cut_v_r   = cut_v_r
        self.cut_v_phi = cut_v_phi

        # scale to cgs units (and account for rotation direction)
        self.pert_v_r_sq   = cut_v_r   * self.p.c_s_planet * (self.p.m_planet/self.p.m_thermal)
        self.pert_v_phi_sq = cut_v_phi * self.p.c_s_planet * (self.p.m_planet/self.p.m_thermal)
        self.pert_rho_sq   = cut_rho   *                     (self.p.m_planet/self.p.m_thermal)

        self.pert_rho_sq_unscaled = cut_rho # this does not need to account for direction

        # account for rotation direction
        if self.p.a_cw == -1:
            self.pert_v_r_sq   =  np.flipud(self.pert_v_r_sq)
            self.pert_v_phi_sq = -np.flipud(self.pert_v_phi_sq)
            self.pert_rho_sq   =  np.flipud(self.pert_rho_sq)

        # save grids
        self.x_cut = x_cut 
        self.y_cut = y_cut
        self.x_sq  = self.p.l * x_cut + self.p.r_planet
        self.y_sq  = self.p.l * y_cut
        
    # extract linear perturbations and interpolate onto annulus segment grid for global results
    def _cut_box_annulus_segment(self) -> None:
        """Extract the part of the linear solution needed for the model, and interpolate onto appropriate grid
        """

        # box size (in units of Hill radius), note for conversions that self.p.l = 1 Hill radius in cgs
        x_box_size = 2*self.p.scale_box
        y_box_size = 2*self.p.scale_box_ang

        # linear perturbations read in grid
        x = self.X[0,:]
        y = self.Y[:,0]

        # cut square box in linear regime
        x_cut = x[np.argmin(x < -x_box_size) : np.argmin(x < x_box_size) + 1]
        y_cut = y[np.argmin(y < -y_box_size) : np.argmin(y < y_box_size) + 1]

        # annulus segment grid, granularity from square box
        r = np.linspace(
            self.p.r_planet - x_box_size*self.p.l, 
            self.p.r_planet + x_box_size*self.p.l, 
            len(x_cut)
        )
        y_ = np.linspace(
            -y_box_size*self.p.l,
            y_box_size*self.p.l,
            len(y_cut)
        )

        R, Y_ = np.meshgrid(r, y_)

        # for the points on our annulus segment, use transformations to find the x,y values on the original perturbations grid
        # if the user has chosen a linear box with very large angular extent, R and Y will be equal at some point, creating an error
        with np.errstate(all='raise'):
            try:

                # pick x values for the box
                if self.p.box_warp:
                    X_pert_grid = (R - self.p.r_planet) / self.p.l
                else:
                    X_pert_grid = (np.sqrt(R**2 - Y_**2) - self.p.r_planet) / self.p.l

                # get phi values
                PHI = np.arctan2(Y_, np.sqrt(R**2 - Y_**2))

            except FloatingPointError:
                print("Error: Reduce 'scale_box_ang' parameter.")
                print("Exiting")
                sys.exit(1)
        
        # unscale Y
        Y_pert_grid = Y_ / self.p.l

        # cut big perturbations grid to just outside annulus
        self.r_min = self.p.r_planet - x_box_size*self.p.l
        self.r_max = self.p.r_planet + x_box_size*self.p.l

        x_min_global = np.sqrt(self.r_min**2 - (y_box_size*self.p.l)**2)
        x_min_local  = (x_min_global - self.p.r_planet) / self.p.l

        # find cut indicies (remember we need to scale back to units of Hill radius )
        x_cut_i1 = np.argmin(x <  x_min_local)
        x_cut_i2 = np.argmin(x <  x_box_size) + 1
        y_cut_i1 = np.argmin(y < -y_box_size)
        y_cut_i2 = np.argmin(y <  y_box_size) + 1

        # cut grid
        x_int_cut = x[x_cut_i1 : x_cut_i2]
        y_int_cut = y[y_cut_i1 : y_cut_i2]

        # cut perturbation arrays
        cut_v_r   = self.pert_v_r   [y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]
        cut_v_phi = self.pert_v_phi [y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]
        cut_rho   = self.pert_rho   [y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]

        if False:
            plt.contourf(x_int_cut, y_int_cut, cut_rho, cmap="RdBu", vmin=-1, vmax=1, levels=100)
            plt.show()

        # account for rotation direction
        if self.p.a_cw == -1:
            cut_v_r   =  np.flipud(cut_v_r)
            cut_v_phi = -np.flipud(cut_v_phi)
            cut_rho   =  np.flipud(cut_rho)

        # interpolation over cut (cartesian) grid
        interp_v_r   = RectBivariateSpline(y_int_cut, x_int_cut, cut_v_r)
        interp_v_phi = RectBivariateSpline(y_int_cut, x_int_cut, cut_v_phi)
        interp_v_rho = RectBivariateSpline(y_int_cut, x_int_cut, cut_rho)

        # evaluate interpolation over our annulus segment
        self.pert_v_r_ann   = interp_v_r.ev  (Y_pert_grid, X_pert_grid)
        self.pert_v_phi_ann = interp_v_phi.ev(Y_pert_grid, X_pert_grid)
        self.pert_rho_ann   = interp_v_rho.ev(Y_pert_grid, X_pert_grid)

        if False:
            plt.imshow(self.pert_rho_ann, cmap="RdBu", vmin=-1, vmax=1)
            plt.show()

        # scale to cgs units
        self.pert_v_r_ann   *= self.p.c_s_planet*(self.p.m_planet/self.p.m_thermal)
        self.pert_v_phi_ann *= self.p.c_s_planet*(self.p.m_planet/self.p.m_thermal)
        self.pert_rho_ann   *=                   (self.p.m_planet/self.p.m_thermal)

        # scale velocities to km/s
        self.pert_v_r_ann   *= 1e-5
        self.pert_v_phi_ann *= 1e-5

        # save annulus grid
        self.r_ann   = r
        self.y_ann   = y_
        self.R_ann   = R
        self.PHI_ann = PHI

        if False:
            # plotting (for debugging)
            _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            myplot = ax.contourf(PHI, R, self.pert_rho_ann/(self.p.m_planet/self.p.m_thermal), levels=300, vmin=-1, vmax=1, cmap='RdBu')
            ax.set_ylim(0, self.p.r_outer)
            plt.colorbar(myplot)
            plt.show()
