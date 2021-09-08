import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
import pymcfost
from mcfost_interface import read_mcfost_grid_data
from matplotlib.colors import LogNorm
import astropy.io.fits as fits

# NOTE: grid dimensions are (x,z,y) or (phi,z,r)

class Grid:
    def __init__(self, parameters):

        # grab parameters object
        self.p = parameters

        # initialise grid properties
        self.info = {
            "Type": None,
            "Log_r": None,
            "Size": [0, 0, 0],
            "Contains": "Empty"
            }

        # initialise arrays
        self.v_r = None 
        self.v_phi = None 
        self.rho = None

    def make_grid(self):

        print(f"Constructing {self.p.grid_type} Grid ")

        # define disk height (not used for mcfost grid)
        self.height = self.p.hr * (self.p.r_outer / self.p.r_ref)**(0.5 - self.p.q) * self.p.r_outer

        # make cartesian grid
        if self.p.grid_type == "cartesian":
            self.make_cartesian_grid()

        # make cylindrical grid
        elif self.p.grid_type == "cylindrical":
            self.make_cylindrical_grid()

        # make mcfost grid
        elif self.p.grid_type == "mcfost":
            self.make_mcfost_grid()

    def make_cartesian_grid(self):
        
        # make grid from specifications in parameter file
        self.x = np.linspace(-self.p.r_outer, self.p.r_outer, self.p.n_x)
        self.y = np.linspace(-self.p.r_outer, self.p.r_outer, self.p.n_y)
        self.z_xy = np.linspace(0, self.height, self.p.n_z)

        self.X, self.Z_xy, self.Y  = np.meshgrid(self.x, self.z_xy, self.y, indexing='ij')

        # update grid info
        self.info["Type"] = "cartesian"
        self.info["Size"][0] = self.x.shape[0]
        self.info["Size"][1] = self.z_xy.shape[0]
        self.info["Size"][2] = self.y.shape[0]

    def make_cylindrical_grid(self):
         
        # make grid from specifications in parameter file
        if self.p.r_log:
            self.r = np.geomspace(self.p.r_inner, self.p.r_outer, self.p.n_r)
        else:
            self.r = np.linspace(self.p.r_inner, self.p.r_outer, self.p.n_r)    
        self.phi = np.linspace(0, 2*np.pi, self.p.n_phi)
        self.z = np.linspace(0, self.height, self.p.n_z)

        self.PHI, self.Z, self.R = np.meshgrid(self.phi, self.z, self.r, indexing='ij')

        # update grid info
        self.info["Type"] = "cylindrical"
        self.info["Size"][0] = self.phi.shape[0]
        self.info["Size"][1] = self.z.shape[0]
        self.info["Size"][2] = self.r.shape[0]
        self.info["Log_r"] = self.p.r_log

    def make_mcfost_grid(self):

        # generate mcfost grid data from specifications in parameter file
        self.R, self.Z = read_mcfost_grid_data(self.p)
        self.r = self.R[0, 0, :]
        self.phi = np.linspace(0, 2*np.pi, self.p.n_phi)

        PHI_temp = np.repeat(self.phi[:, np.newaxis], self.p.n_r, axis=1)
        self.PHI = np.repeat(PHI_temp[:, np.newaxis, :], self.p.n_z, axis=1)

        # note that we do not have a self.z because it depends on radius, as the grid is flared

        # update grid info
        self.info["Type"] = "mcfost"
        self.info["Size"][0] = self.phi.shape[0]
        self.info["Size"][1] = self.Z.shape[1]
        self.info["Size"][2] = self.r.shape[0]
        self.info["Log_r"] = "mcfost"

    def get_r_phi_coords(self):

        # need the equivalent (r,phi) coordinates when using cartesian grid
        self.R_xy = np.sqrt(self.X**2 + self.Y**2)
        self.PHI_xy = np.arctan2(self.Y, self.X)

    def make_keplerian_disk(self):

        print("Making Keplerian disk ")

        # get radii
        if self.info["Type"] == "cartesian":
            self.get_r_phi_coords()
            r = self.R_xy
            z = self.Z_xy
        else:
            r = self.R
            z = self.Z

        p = self.p.p
        q = self.p.q
        hr = self.p.hr
        r_ref = self.p.r_ref
        rho_ref = self.p.rho_ref
        c_s_0 = self.p.c_s_0
        G = self.p.G_const
        m_star = self.p.m_star
        m_sol = self.p.m_solar
        au = self.p.au

        # Keplerian velocities
        self.v_r = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))
        self.v_kep = 1e-5 * np.sqrt(self.p.G_const * self.p.m_star * self.p.m_solar / (r * self.p.au))
        self.v_phi = np.copy(self.v_kep) * self.p.a_cw

        # pressure and gravity height dependent correction
        corr = np.sqrt(
            -1* (p + q) * hr**2 * (r / r_ref)**(1 - 2*q) + (1 - q) + q*r / np.sqrt(r**2 + z**2)
        )    

        # pressure correction for midplane
        """corr  = np.sqrt( 1 - (2 * q + p) * hr**2 * (r / r_ref)**(1 - 2*q) )"""

        # perform correction
        self.v_phi *= corr

        # unperturbed density profle (full)
        """self.rho = (
            rho_ref * (r / r_ref)**-p * np.exp(
                ((self.p.G_const*self.p.m_star*self.p.m_solar) / (hr**2 * (r / r_ref)**(1 - 2*q) * self.v_kep**2)) \
                    * (1 / np.sqrt((r*self.p.au)**2 + (z*self.p.au)**2) + 1 / (r*self.p.au))
            )
        )"""
        c_s = c_s_0 * (r/r_ref)**-q

        self.rho = (
            rho_ref * (r/r_ref)**-p * np.exp(
                (G*m_star*m_sol / c_s**2) * (1/np.sqrt((r*au)**2 + (z*au)**2) - 1/(r*au))
            )
        )

        # unperturbed density profile (thin disk approx.)
        #self.rho = rho_ref * (r / r_ref)**-p * np.exp( - 0.5 * (z / r_ref)**2 * hr**-2 * (r / r_ref)**(2 * q - 3))

        # update grid info
        self.info["Contains"] = "Keplerian velocity field, unperturbed density"

    def make_empty_disk(self):

        # make grids full of zeros
        self.v_r = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))
        self.v_phi = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))
        self.rho = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))

        # update grid info
        self.info["Contains"] = "Zeros"

    def add_linear_perturbations(self, LinearPerts, rho_background):

        # box size (in units of Hill radius), note for conversions that self.p.l = 1 Hill radius in cgs
        box_size = 2*self.p.scale_box
        min_phi = -self.p.scale_box_ang*np.arcsin(box_size*self.p.l / self.p.r_planet)
        max_phi = self.p.scale_box_ang*np.arcsin(box_size*self.p.l / self.p.r_planet)
        min_r = self.p.r_planet - box_size*self.p.l
        max_r = self.p.r_planet + box_size*self.p.l

        # find (phi, r) grid for either Cartesian or Cylindrical global grid
        if self.info["Type"] == "cartesian":
            self.get_r_phi_coords()
            R, PHI = self.R_xy, self.PHI_xy
        else:
            R, PHI = self.R, self.PHI

        # new PHI grid to use (-pi,pi) instead of (0,2pi), where values are swapped in place, also ditch z coordinate
        # also construct a mask that contains 0 outside linear annulus and 1 inside
        PHI_new = np.zeros((PHI.shape[0],PHI.shape[2]))
        R_new = R[:,0,:]
        linear_mask = np.zeros((PHI.shape[0],PHI.shape[2]))

        for i in range(PHI.shape[0]):
            for j in range(PHI.shape[2]):

                # transforming phi coordinate in place
                if PHI[i,0,j] > np.pi:
                    PHI_new[i,j] = PHI[i,0,j] - 2*np.pi
                else:
                    PHI_new[i,j] = PHI[i,0,j]

                # constructing mask
                if PHI_new[i,j] > min_phi and PHI_new[i,j] < max_phi and R_new[i,j] > min_r and R_new[i,j] < max_r:
                    linear_mask[i,j] = 1
                else:
                    linear_mask[i,j] = 0
                

        # get linear solution           
        lp = LinearPerts

        # assemble interpolation functions over linear perts grid
        interp_v_r = RectBivariateSpline(lp.phi_ann, lp.r_ann, lp.pert_v_r_ann)
        interp_v_phi = RectBivariateSpline(lp.phi_ann, lp.r_ann, lp.pert_v_phi_ann)
        interp_rho = RectBivariateSpline(lp.phi_ann, lp.r_ann, lp.pert_rho_ann)

        # evaluate interpolations on global grid
        global_lin_v_r = interp_v_r.ev(PHI_new, R_new)
        global_lin_v_phi = interp_v_phi.ev(PHI_new, R_new)
        global_lin_rho = interp_rho.ev(PHI_new, R_new)

        # apply mask to only get solution in valid domain
        global_lin_v_r = global_lin_v_r * linear_mask
        global_lin_v_phi = global_lin_v_phi * linear_mask
        global_lin_rho = global_lin_rho * linear_mask

        # Add velocities, identical at all heights for now
        self.v_r += global_lin_v_r[:, np.newaxis, :]
        self.v_phi += global_lin_v_phi[:, np.newaxis, :]

        # Add density, scaling by background density
        self.rho = global_lin_rho[:, np.newaxis, :] * rho_background

        # plot for debugging
        """
        _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        myplot = ax.contourf(PHI[:,0,:], R[:,0,:], global_lin_v_r*linear_mask, levels=300, vmin=-2E4, vmax=2E4, cmap='RdBu')
        ax.set_ylim(0, self.p.r_outer)
        plt.colorbar(myplot)
        plt.show()
        """

        # update grid info
        self.info["Contains"] = "linear perturbations"

    def add_non_linear_perturbations(self, NonLinearPerts, rho_background):

        nonlin = NonLinearPerts

        # Add velocities, identical at all heights for now
        nl_vr = nonlin.vr[:, np.newaxis, :]
        nl_vphi = nonlin.vphi[:, np.newaxis, :]

        self.v_r += nl_vr
        self.v_phi += nl_vphi

        # Define scale height array
        if self.info["Type"] == "cartesian":

            self.get_r_phi_coords()
            r = self.R_xy
            z = self.Z_xy

            size_ones = (len(self.x), self.p.n_z, len(self.y))

        else:

            r = self.R 
            z = self.Z

            size_ones = (self.p.n_phi, self.p.n_z, self.p.n_r)

        self.H = self.p.h_ref * (r / self.p.r_ref)**self.p.beta

        # Add density, scaling vertically as per thin disk assumption
        #nl_rho = (np.ones(size_ones) + nonlin.rho[:, np.newaxis, :]) * self.p.rho_ref * (r/self.p.r_ref)**(-self.p.p)
        #self.rho = nl_rho * np.exp(-0.5 * (z / self.H)**2)

        self.rho = nonlin.rho[:, np.newaxis, :] * rho_background

        # update grid info
        self.info["Contains"] = "non-linear perturbations"

    def add_phantom_midplane(self, PhantomDump):

        PD = PhantomDump

        # Add velocities, identical at all heights
        ph_vr = PD.vr[:, np.newaxis, :]
        ph_vphi = PD.vphi[:, np.newaxis, :]

        self.v_r += ph_vr #/(1 + self.Z)
        self.v_phi += ph_vphi #/(1 + self.Z)

        # Define scale height array
        self.H = self.p.h_ref * (self.R / self.p.r_ref)**self.p.beta

        # Add density, scaling vertically as per thin disk assumption
        ph_rho = PD.rho[:, np.newaxis, :]
        self.rho = ph_rho * np.exp(-0.5 * (self.Z / self.H)**2)

        # update grid info
        self.info["Contains"] = "phantom mid-plane extrapolated to 3D"

    def merge_grids(self, grid_to_merge):

        g = grid_to_merge

        if type(grid_to_merge) is not Grid:
            print("Error: merge_grids must be given a Grid object. Merging failed.")
            return False
        
        if self.p != g.p:
            print("Error: merge_grids requires both grids must have same parameters. Merging failed.")
            return False

        # merge data arrays
        self.v_r += g.v_r
        self.v_phi += g.v_phi
        self.rho += g.rho

        # update info
        self.info["Contains"] += " AND " + g.info["Contains"]

    def merge_phantom_densities(self, grid_to_merge):

        g = grid_to_merge

        if type(grid_to_merge) is not Grid:
            print("Must be given a Grid object")
            return False
        
        if self.p != g.p:
            print("Both grids must have same type parameters")
            return False

        # merge data arrays
        self.rho = g.rho

    def show_disk2D(self, z_slice, save=False, name='disk2Dplot'):

        # plot v_r
        if self.info["Type"] == "cartesian":
            plt.contourf(self.X[:,0,0], self.Y[0,0,:], self.v_r[:,z_slice,:], levels=300, cmap="RdBu")
            plt.colorbar()
        else:
            _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            myplot = ax.contourf(self.PHI[:,0,:], self.R[:,0,:], self.v_r[:,z_slice,:], levels=300, cmap='RdBu')
            ax.set_ylim(0, self.p.r_outer)
            plt.colorbar(myplot)
        plt.title(r"$v_r$")
        if save:
            plt.savefig(f'{self.p.system}/{self.p.name}/{name}_vr_z{z_slice}.pdf')
        plt.show()

        # plot v_phi
        if self.info["Type"] == "cartesian":
            plt.contourf(self.X[:,0,0], self.Y[0,0,:], self.v_phi[:,z_slice,:], levels=300)
            plt.colorbar()
        else:
            _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            myplot = ax.contourf(self.PHI[:,0,:], self.R[:,0,:], self.v_phi[:,z_slice,:], levels=300, cmap='RdBu')
            ax.set_ylim(0, self.p.r_outer)
            plt.colorbar(myplot)
        plt.title(r"$v_{\phi}$")
        if save:
            plt.savefig(f'{self.p.system}/{self.p.name}/{name}_vphi_z{z_slice}.pdf')
        plt.show()

        # plot rho
        if self.info["Type"] == "cartesian":
            plt.contourf(self.X[:,0,0], self.Y[0,0,:], self.rho[:,z_slice,:], levels=300)
            plt.colorbar()
        else:
            _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            myplot = ax.contourf(self.PHI[:,0,:], self.R[:,0,:], self.rho[:,z_slice,:], levels=300)
            ax.set_ylim(0, self.p.r_outer)
            plt.colorbar(myplot)
        plt.title(r"$\rho$")
        if save:
            plt.savefig(f'{self.p.system}/{self.p.name}/{name}_rho_z{z_slice}.pdf')
        plt.show()
    
    def write_fits_file(self):

        if self.info["Type"] == "cartesian":

            print("Cannot write FITS file when using Cartesian grid")

        else:

            # create empty array for vertical velocities
            self.v_z = np.zeros((self.p.n_phi, self.p.n_z, self.p.n_r))

            # create master velocities array
            velocities = 1e3 * np.array([self.v_r, self.v_phi, self.v_z])

            # setup HDUs for fits file
            primary_hdu = fits.PrimaryHDU(np.abs(self.rho))
            second_hdu = fits.ImageHDU(np.abs(self.rho)) # * np.exp(-1 * (self.Z / self.H)**2)
            tertiary_hdu = fits.ImageHDU(velocities)

            # set header properties for mcfost
            primary_hdu.header['hierarch read_gas_velocity'] = 2
            primary_hdu.header['hierarch gas_dust_ratio'] = 100
            primary_hdu.header['hierarch read_gas_density'] = 1
            primary_hdu.header['read_n_a'] = 0

            # sink particles properties in header (code modifications needed for anything other than 2)
            primary_hdu.header['hierarch N_sink'] = 2

            # star properties in header
            primary_hdu.header['hierarch M_star'] = self.p.m_star        # star mass in M_sol
            primary_hdu.header['hierarch x_star'] = 0.                   # star at origin
            primary_hdu.header['hierarch y_star'] = 0.
            primary_hdu.header['hierarch z_star'] = 0.
            primary_hdu.header['hierarch vx_star'] = 0.                  # star not moving
            primary_hdu.header['hierarch vy_star'] = 0.
            primary_hdu.header['hierarch vz_star'] = 0.

            # planet properties in header
            primary_hdu.header['hierarch M_planet'] = self.p.m_planet    # planet mass in M_jup
            primary_hdu.header['hierarch x_planet'] = self.p.r_planet    # planet x position simply orbital radius [AU]
            primary_hdu.header['hierarch y_planet'] = 0.
            primary_hdu.header['hierarch z_planet'] = 0.
            primary_hdu.header['hierarch vx_planet'] = 0.                # planet velocity just keplerian rotation [m/s] (positive for anticlockwise, which is cw=-1)
            primary_hdu.header['hierarch vy_planet'] = np.sqrt(
                self.p.G_const*self.p.m_star*self.p.m_solar/(self.p.r_planet*self.p.au)
                ) / 100 * self.p.a_cw
            primary_hdu.header['hierarch vz_planet'] = 0.

            # setup fits file
            hdul = fits.HDUList([primary_hdu, second_hdu, tertiary_hdu])

            # Write a fits file for mcfost
            print("Writing a fits file for MCFOST")
            fitsname = "wakeflow_model.fits"
            hdul.writeto(f"{self.p.system}/{self.p.name}/mcfost_output/{fitsname}", overwrite=True)