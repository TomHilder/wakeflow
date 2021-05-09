import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
import pymcfost
from mcfost_interface import read_mcfost_grid_data

# important note: grid dimensions are (x,z,y) or (r,z,phi)

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

        self.R, self.Z, self.PHI = np.meshgrid(self.r, self.z, self.phi, indexing='ij')

        # update grid info
        self.info["Type"] = "cylindrical"
        self.info["Size"][0] = self.r.shape[0]
        self.info["Size"][1] = self.z.shape[0]
        self.info["Size"][2] = self.phi.shape[0]
        self.info["Log_r"] = self.p.r_log

    def make_mcfost_grid(self):

        # generate mcfost grid data from specifications in parameter file
        self.R, self.Z = read_mcfost_grid_data(self.p)
        self.r = self.R[:, 0, 0]
        self.phi = np.linspace(0, 2*np.pi, self.p.n_phi)

        PHI_temp = np.repeat(self.phi[np.newaxis, :], self.p.n_r, axis=0)
        self.PHI = np.repeat(PHI_temp[:, np.newaxis, :], self.p.n_z, axis=1)

        # note that we do not have a self.z because it depends on radius, as the grid is flared

        # update grid info
        self.info["Type"] = "mcfost"
        self.info["Size"][0] = self.r.shape[0]
        self.info["Size"][1] = self.Z.shape[1]
        self.info["Size"][2] = self.phi.shape[0]
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

        # Keplerian velocities
        self.v_r = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))
        self.v_kep = np.sqrt(self.p.G_const * self.p.m_star * self.p.m_solar / (r * self.p.au))
        self.v_phi = np.copy(self.v_kep)

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

        # unperturbed density profile (thin disk approx.)
        self.rho = rho_ref * (r / r_ref)**-p * np.exp( - 0.5 * (z / r_ref)**2 * hr**-2 * (r / r_ref)**(2 * q - 3))

        # update grid info
        self.info["Contains"] = "Keplerian velocity field, unperturbed density"

    def make_empty_disk(self):

        # make grids full of zeros
        self.v_r = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))
        self.v_phi = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))
        self.rho = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))

        # update grid info
        self.info["Contains"] = "Zeros"

    def add_linear_perturbations(self, linear_perts):

        # find (phi, r) grid for either Cartesian or Cylindrical global grid
        if self.info["Type"] == "cartesian":
            self.get_r_phi_coords()
            R, PHI = self.R_xy, self.PHI_xy
        elif self.p.grid_type == "cylindrical":
            R, PHI = self.R, self.PHI
             

        lin = linear_perts

        # assemble interpolation functions over linear perts grid
        interp_v_r = RectBivariateSpline(lin.phi, lin.r, lin.pert_v_r) #should do transpose of pert_v_r and swap phi and r

        # find part of grid inside linear domain
        #r_cut_i1 = np.argmin()

        # evaluate relevent part of grid using interpolation function and map onto grid

        # extrapolate along z-axis (copy velocities, scale densities)

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