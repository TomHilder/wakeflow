# linear_perts.py
# Written by Thomas Hilder, Daniele Fasano and Francesco Bollati

"""
Contains the LinearPerts class responsible for handling the linear regime of the models.
"""

# type hinting without circular imports
from typing                 import TYPE_CHECKING

import sys, pkg_resources, tarfile
import numpy                    as np
import matplotlib.pyplot        as plt
from scipy.interpolate      import RectBivariateSpline
from .transformations       import _Eta_vector

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
        
        # initialise perturbations properties
        self.info = {
            "Type": None,
            "Grid": None,
            "Size": [0, 0, 0]
            }
        
        # read in perturbations from data files
        mesh, perts = read_perturbation_files(self.p.lin_type)
        
        # get perturbation arrays
        self.pert_v_r   = perts[0]
        self.pert_v_phi = perts[1]
        self.pert_rho   = perts[2]
        
        # define constants for linear perts
        self.x_box_left   = 2 * self.p.scale_box_left
        self.x_box_right  = 2 * self.p.scale_box_right
        self.x_box_top    = 2 * self.p.scale_box_ang_top
        self.x_box_bottom = 2 * self.p.scale_box_ang_bottom
        
        # Using linear perturbations computed following Bollati et al. 2021 with shearing sheet assumption
        # NOTE: this is deprecated and will probably lead to wrong results
        if self.p.lin_type == "shearing_sheet":
            # grid
            self.X = mesh[0]
            self.Y = mesh[1]
            # linear perturbations read in grid
            x = self.X[0,:]
            y = self.Y[:,0]
            # cut square box grid in linear regime
            self.x_cut = x[np.argmin(x < -self.x_box_left)   : np.argmin(x < self.x_box_right) + 1]
            self.y_cut = y[np.argmin(y < -self.x_box_bottom) : np.argmin(y < self.x_box_top) + 1]
            #updating info of linear perts        
            self.info["Type"]    = "shearing_sheet"
            self.info["Grid"]    = "cartesian"
            self.info["Size"][0] = x.shape[0]
            self.info["Size"][1] = y.shape[0]
            self.info["Size"][2] = self.p.n_z
            
        #Using global linear perturbations computed using Miranda's code.         
        elif self.p.lin_type == "global":
            # grid
            self.R   = mesh[0]
            self.PHI = mesh[1]
            #Rescale R with planet radius
            self.R *= self.p.r_planet
            # linear perturbations read in grid
            r   = self.R[0,:]
            phi = self.PHI[:,0]
            #defining cartesian grid
            x = np.linspace(-np.max(r), np.max(r), len(r))
            y = np.linspace(-np.max(r), np.max(r), len(r))
            #creating cartesian mesh
            self.X, self.Y = np.meshgrid(x, y)
            #updating info of linear perts      
            self.info["Type"]    = "global"
            self.info["Grid"]    = "cylindrical"
            self.info["Size"][0] = r.shape[0]
            self.info["Size"][1] = phi.shape[0]
            self.info["Size"][2] = self.p.n_z
            
        #Using perturbations from hydro simulation (not yet implemented)
        elif self.p.lin_type == "simulation":
            raise NotImplementedError("Reading perturbations from a simulation is not yet supported.")

    # old method of extracting linear perturbations 
    def _cut_box_square(self) -> None:
        """Deprecated. This is only used in the case where you want to plot the linear solution in t, eta space
        """

        # box size (in units of Hill radius), with default scale_box = 1. (note for conversions that self.p.l = 1 Hill radius in cgs)
        box_size = 2*self.p.scale_box_left
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
    def _cut_annulus_segment(self) -> None:
        """Extract the part of the linear solution needed for the model, and interpolate onto appropriate grid
        """

        # segment radial size (in units of Hill radius), note for conversions that self.p.l = 1 Hill radius in cgs. 
        r_box_size_left  = self.x_box_left
        r_box_size_right = self.x_box_right
        # segment azimuthal size (in units of \pi). 
        phi_box_size_top    = self.x_box_top / 2
        phi_box_size_bottom = self.x_box_bottom / 2

        #interpolate perturbations on a cylindrical grid and evaluate them in the annulus segment
        if self.info["Type"] == "shearing_sheet" or self.info["Type"] == "simulation":

            # linear perturbations read in grid
            x = self.X[0,:]
            y = self.Y[:,0]

            # cut square box in linear regime
            x_cut = x[np.argmin(x < -r_box_size_left) : np.argmin(x < r_box_size_right) + 1]

            # annulus segment grid, radial granularity from square box, angular granularity fixed for now
            r = np.linspace(
                self.p.r_planet - r_box_size_left*self.p.l, 
                self.p.r_planet + r_box_size_right*self.p.l, 
                len(x_cut)
            )
            
            phi = np.linspace(
                -phi_box_size_bottom*np.pi,
                phi_box_size_top*np.pi,
                1000
            )

            R, PHI = np.meshgrid(r, phi)
            
            # preparing pertubations for interpolation
            v_r_cart   = self.pert_v_r   
            v_phi_cart = self.pert_v_phi 
            rho_cart   = self.pert_rho
            
            #flip perts  if rotation is clockwise
            if self.p.a_cw == -1:
                v_r_cart   =  np.flipud(v_r_cart)
                v_phi_cart = -np.flipud(v_phi_cart)
                rho_cart   =  np.flipud(rho_cart)
                
            # interpolation over global linear (cartesian) grid
            interp_v_r   = RectBivariateSpline(y, x, v_r_cart)
            interp_v_phi = RectBivariateSpline(y, x, v_phi_cart)
            interp_v_rho = RectBivariateSpline(y, x, rho_cart)

            #evaluation of cartesian coordinates corresponding to our polar grid
            X_pert_grid = (R * np.cos(PHI) - self.p.r_planet)/self.p.l
            Y_pert_grid = R * np.sin(PHI)/self.p.l

            # evaluate interpolation over our annulus segment
            self.pert_v_r_ann   = interp_v_r.ev  (Y_pert_grid, X_pert_grid)
            self.pert_v_phi_ann = interp_v_phi.ev(Y_pert_grid, X_pert_grid)
            self.pert_rho_ann   = interp_v_rho.ev(Y_pert_grid, X_pert_grid)

            #plot for debugging
            if False:
                plt.imshow(self.pert_rho_ann, cmap="RdBu", vmin=-1, vmax=1, origin='lower', extent=(r[0], r[-1], phi[0], phi[-1]))
                plt.axis('auto')
                plt.xlabel('R [au]')
                plt.ylabel(r'$\varphi$ [rad]')
                plt.show()
            if False:
                eta = _Eta_vector(self.p.r_planet + r_box_size_left*self.p.l, phi, self.p.r_planet, self.p.hr, self.p.q, self.p.p, self.p.cw_rotation, self.p.m_planet, self.p.m_thermal, self.p.nl_wake)
                plt.plot(phi, self.pert_rho_ann[:,-1])
                ax = plt.gca()
                ax.set_xlabel(r'$\varphi$ [rad]')
                ax.set_ylabel(r'$\sigma$')
                ax2 = ax.twiny()
                ax2.plot(eta, self.pert_rho_ann[:,-1])
                #ax2.plot(eta, self.pert_rho_ann[:,np.argmin(r_<self.p.r_planet + x_box_size_r*self.p.l)])
                ax2.set_xlabel(r'$\eta$')
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
            self.phi_ann = phi
            self.R_ann   = R
            self.PHI_ann = PHI

        #The global solution is computed on a cylindrical grid, no need to interpolate
        elif self.info["Type"] == "global":

            # linear perturbations read in grid
            r   = self.R[0,:]
            phi = self.PHI[:,0]

            #masks for restriction
            r_mask = np.logical_and(r >= self.p.r_planet - r_box_size_left*self.p.l, r <= self.p.r_planet + r_box_size_right*self.p.l)
            phi_mask = np.logical_and(phi >= -phi_box_size_bottom*np.pi, phi <= phi_box_size_top*np.pi)
            #print(r_mask.shape, phi_mask.shape)
            #restricting grid to the annulus segment
            r_ann   = r[r_mask]
            phi_ann = phi[phi_mask]

            R_ann, PHI_ann = np.meshgrid(r_ann, phi_ann)

            # preparing pertubations for annulus segment
            v_r_cyl   = self.pert_v_r.T   
            v_phi_cyl = self.pert_v_phi.T
            rho_cyl   = self.pert_rho.T

            #plot for debugging
            if False:
                plt.imshow(rho_cyl, cmap="RdBu", vmin=-1, vmax=1, origin='lower', extent=(r_ann[0], r_ann[-1], phi_ann[0], phi_ann[-1]))
                plt.axis('auto')
                plt.xlabel('R [au]')
                plt.ylabel(r'$\varphi$ [rad]')
                plt.show()
                plt.imshow(v_r_cyl, cmap="RdBu", vmin=-1, vmax=1, origin='lower', extent=(r_ann[0], r_ann[-1], phi_ann[0], phi_ann[-1]))
                plt.axis('auto')
                plt.xlabel('R [au]')
                plt.ylabel(r'$\varphi$ [rad]')
                plt.show()
                plt.imshow(v_phi_cyl, cmap="RdBu", vmin=-1, vmax=1, origin='lower', extent=(r_ann[0], r_ann[-1], phi_ann[0], phi_ann[-1]))
                plt.axis('auto')
                plt.xlabel('R [au]')
                plt.ylabel(r'$\varphi$ [rad]')
                plt.show()
            #print(v_r_cyl.shape)

            #flip perts  if rotation is clockwise
            if self.p.a_cw == -1:
                v_r_cyl   =  np.flipud(v_r_cyl)
                v_phi_cyl = -np.flipud(v_phi_cyl)
                rho_cyl   =  np.flipud(rho_cyl)

            #restricting perturbations to the annulus segment
            self.pert_v_r_ann   = v_r_cyl[:,r_mask][phi_mask,:]
            self.pert_v_phi_ann = v_phi_cyl[:,r_mask][phi_mask,:]
            self.pert_rho_ann   = rho_cyl[:,r_mask][phi_mask,:]

            #plot for debugging
            if False:
                plt.imshow(self.pert_rho_ann, cmap="RdBu", vmin=-1, vmax=1, origin='lower', extent=(r_ann[0], r_ann[-1], phi_ann[0], phi_ann[-1]))
                plt.axis('auto')
                plt.xlabel('R [au]')
                plt.ylabel(r'$\varphi$ [rad]')
                plt.show()
                plt.imshow(self.pert_v_r_ann, cmap="RdBu", vmin=-1, vmax=1, origin='lower', extent=(r_ann[0], r_ann[-1], phi_ann[0], phi_ann[-1]))
                plt.axis('auto')
                plt.xlabel('R [au]')
                plt.ylabel(r'$\varphi$ [rad]')
                plt.show()
                plt.imshow(self.pert_v_phi_ann, cmap="RdBu", vmin=-1, vmax=1, origin='lower', extent=(r_ann[0], r_ann[-1], phi_ann[0], phi_ann[-1]))
                plt.axis('auto')
                plt.xlabel('R [au]')
                plt.ylabel(r'$\varphi$ [rad]')
                plt.show()
            if False:
                eta = _Eta_vector(self.p.r_planet + r_box_size_left*self.p.l, phi_ann, self.p.r_planet, self.p.hr, self.p.q, self.p.p, self.p.cw_rotation, self.p.m_planet, self.p.m_thermal, self.p.nl_wake)
                plt.plot(phi_ann, self.pert_rho_ann[:,-1])
                ax = plt.gca()
                ax.set_xlabel(r'$\varphi$ [rad]')
                ax.set_ylabel(r'$\sigma$')
                ax2 = ax.twiny()
                ax2.plot(eta, self.pert_rho_ann[:,-1])
                #ax2.plot(eta, self.pert_rho_ann[:,np.argmin(r_<self.p.r_planet + x_box_size_r*self.p.l)])
                ax2.set_xlabel(r'$\eta$')
                plt.show()

            # scale to cgs units
            self.pert_v_r_ann   *= self.p.c_s_planet*(self.p.m_planet/self.p.m_thermal)
            self.pert_v_phi_ann *= self.p.c_s_planet*(self.p.m_planet/self.p.m_thermal)
            self.pert_rho_ann   *=                   (self.p.m_planet/self.p.m_thermal)

            # scale velocities to km/s
            self.pert_v_r_ann   *= 1e-5
            self.pert_v_phi_ann *= 1e-5

            # save annulus grid
            self.r_ann   = r_ann
            self.phi_ann = phi_ann
            self.R_ann   = R_ann
            self.PHI_ann = PHI_ann



def read_perturbation_files(lin_type="global"):
    # get location of linear perturbations data files
    if lin_type == "global":
        pert_loc = pkg_resources.resource_filename('wakeflow', 'data/global_linear_perturbations.npy')
        mesh_loc = pkg_resources.resource_filename('wakeflow', 'data/global_linear_perturbations_mesh.npy')
    elif lin_type == "shearing_sheet":
        pert_loc = pkg_resources.resource_filename('wakeflow', 'data/linear_perturbations.npy')
        mesh_loc = pkg_resources.resource_filename('wakeflow', 'data/linear_perturbations_mesh.npy')
    elif lin_type == "simulation":
        raise NotImplementedError("Reading the perturbations from simulation data is not supported yet.")
    else:
        raise ValueError("lin_type must be either 'global' or 'shearing_sheet'")
    # try to read perturbations from files
    try:
        perts = np.load(pert_loc)
        mesh  = np.load(mesh_loc)
    # if files are note unzipped yet, do that
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
    return mesh, perts