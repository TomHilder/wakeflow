import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
import pymcfost
from mcfost_interface import read_mcfost_grid_data
from matplotlib.colors import LogNorm
import astropy.io.fits as fits

# NOTE: grid dimensions are (x,z,y) or (phi,z,r)

# TODO: Add mapping for linear perturbations into grid object arrays, for both square box and annulus segment
# TODO: Add vertical scaling for adding non-linear perturbations

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
        self.z = np.linspace(0, self.height, self.p.n_z)

        self.X, self.Z, self.Y = np.meshgrid(self.x, self.z, self.y, indexing='ij')

        # update grid info
        self.info["Type"] = "cartesian"
        self.info["Size"][0] = self.x.shape[0]
        self.info["Size"][1] = self.z.shape[0]
        self.info["Size"][2] = self.y.shape[0]

    def make_cylindrical_grid(self):
         
        # make grid from specifications in parameter file
        if self.p.r_log:
            self.r = np.geomspace(self.p.r_inner, self.p.r_outer, self.p.n_r)
        else:
            self.r = np.linspace(self.p.r_inner, self.p.r_outer, self.p.n_r)    
        self.phi = np.linspace(0, 2*np.pi, self.p.n_phi)
        self.z = np.linspace(0, self.height, self.p.n_z)

        self.R, self.Z, self.PHI = np.meshgrid(self.phi, self.z, self.r, indexing='ij')

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

        # often need the equivalent (r,phi) coordinates when using cartesian grid
        self.R_xy = np.sqrt(self.X**2 + self.Y**2)
        self.PHI_xy = np.arctan(self.Y / self.X)

    def make_keplerian_disk(self):

        print("Making Keplerian disk ")

        # get radii
        if self.info["Type"] == "cartesian":
            self.get_r_phi_coords()
            r = self.R_xy
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

        print(self.v_phi)

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

    def add_linear_perturbations(self, LinearPerts):

        pass
        """
        # find (phi, r) grid for either Cartesian or Cylindrical global grid
        if self.info["Type"] == "cartesian":
            self.get_r_phi_coords()
            R, PHI = self.R_xy, self.PHI_xy
        elif self.p.grid_type == "cylindrical":
            R, PHI = self.R, self.PHI
             

        lin = LinearPerts

        # assemble interpolation functions over linear perts grid
        interp_v_r = RectBivariateSpline(lin.phi, lin.r, lin.pert_v_r) #should do transpose of pert_v_r and swap phi and r

        # find part of grid inside linear domain
        #r_cut_i1 = np.argmin()

        # evaluate relevent part of grid using interpolation function and map onto grid

        # extrapolate along z-axis (copy velocities, scale densities)
        """

    def add_non_linear_perturbations(self, NonLinearPerts):

        nonlin = NonLinearPerts

        # Add velocities, identical at all heights for now
        nl_vr = nonlin.vr[:, np.newaxis, :]
        nl_vphi = nonlin.vphi[:, np.newaxis, :]

        self.v_r += nl_vr #/(1 + self.Z)
        self.v_phi += nl_vphi #/(1 + self.Z)

        # Define scale height array
        H = self.p.h_ref * (self.R / self.p.r_ref)**self.p.beta

        # Add density, scaling vertically as per thin disk assumption
        nl_rho = (np.ones((self.p.n_phi, self.p.n_z, self.p.n_r)) + nonlin.rho[:, np.newaxis, :]) * self.p.rho_ref * (self.R/self.p.r_ref)**(-self.p.p)
        self.rho = nl_rho * np.exp(-0.5 * (self.Z / H)**2)

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
        H = self.p.h_ref * (self.R / self.p.r_ref)**self.p.beta

        # Add density, scaling vertically as per thin disk assumption
        ph_rho = PD.rho[:, np.newaxis, :]
        self.rho = ph_rho * np.exp(-0.5 * (self.Z / H)**2)

        # update grid info
        self.info["Contains"] = "phantom mid-plane extrapolated to 3D"

    def merge_grids(self, grid_to_merge):

        g = grid_to_merge

        if type(grid_to_merge) is not Grid:
            print("Must be given a Grid object")
            return False
        
        if self.p != g.p:
            print("Both grids must have same type parameters")
            return False

        # merge data arrays
        self.v_r += g.v_r
        self.v_phi += g.v_phi
        self.rho += g.rho

        # update info
        self.info["Contains"] += " AND " + g.info["Contains"]

    def show_disk2D(self, z_slice):

        # plot v_r
        plt.imshow(self.v_r[:,z_slice,:])
        plt.colorbar()
        plt.title(r"$v_r$")
        plt.show()

        # plot v_phi
        plt.imshow(self.v_phi[:,z_slice,:])
        plt.colorbar()
        plt.title(r"$v_{\phi}$")
        plt.show()

        # plot rho
        plt.imshow(self.rho[:,z_slice,:])
        plt.colorbar()
        plt.title(r"$\rho$")
        plt.show()
    
    def write_fits_file(self):

        # create empty array for vertical velocities
        self.v_z = np.zeros((self.p.n_phi, self.p.n_z, self.p.n_r))

        # create master velocities array
        velocities = 1e3 * np.array([self.v_r, self.v_phi, self.v_z])

        print(np.shape(self.v_r), np.shape(self.v_phi), np.shape(self.v_z))
        print(np.shape(self.rho))
        print(np.shape(velocities))

        # setup HDUs for fits file
        primary_hdu = fits.PrimaryHDU(np.abs(self.rho))
        second_hdu = fits.ImageHDU(np.abs(self.rho))
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
        fitsname = "analytic_disk.fits"
        hdul.writeto(f"{self.p.system}/{self.p.name}/mcfost_output/{fitsname}", overwrite=True)