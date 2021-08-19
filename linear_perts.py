import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import ticker, cm

class LinearPerts():
    def __init__(self, parameters):

        print("Reading in Linear Perturbations")

        # grab parameters object
        self.p = parameters

        # read perturbations from files
        perts = np.load("linear_perturbations.npy")
        mesh = np.load("linear_perturbations_mesh.npy")

        # get perturbation arrays
        self.pert_v_r = perts[0]
        self.pert_v_phi = perts[1]
        self.pert_rho = perts[2]
        self.X = mesh[0]
        self.Y = mesh[1]

        # define constants for linear perts
        self.x_box = 2*self.p.scale_box
        # box size
        # cut off

    def cut_box_square(self):

        print("Extracting Linear Perturbations in the Vicinity of the Planet")

        # box size (in units of Hill radius), with default scale_box = 1. (note for conversions that self.p.l = 1 Hill radius in cgs)
        box_size = 2*self.p.scale_box

        # linear perturbations read in grid
        x = self.X[0,:]
        y = self.Y[:,0]

        # cut square box grid in linear regime
        x_cut = x[np.argmin(x < -box_size) : np.argmin(x < box_size) + 1]
        y_cut = y[np.argmin(y < -box_size) : np.argmin(y < box_size) + 1]

        # find cut indicies 
        x_cut_i1 = np.argmin(x < -box_size)
        x_cut_i2 = np.argmin(x < box_size) + 1
        y_cut_i1 = np.argmin(y < -box_size)
        y_cut_i2 = np.argmin(y < box_size) + 1

        # cut perturbation arrays
        cut_v_r = self.pert_v_r[y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]
        cut_v_phi = self.pert_v_phi[y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]
        cut_rho = self.pert_rho[y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]

        self.cut_rho = cut_rho

        # scale to cgs units (and account for rotation direction)
        self.pert_v_r_sq = cut_v_r * self.p.c_s_planet*(self.p.m_planet/self.p.m_thermal)
        self.pert_v_phi_sq = cut_v_phi * self.p.c_s_planet*(self.p.m_planet/self.p.m_thermal)
        self.pert_rho_sq  = cut_rho * (self.p.m_planet/self.p.m_thermal)

        self.pert_rho_sq_unscaled = cut_rho # this does not need to account for direction

        # account for rotation direction
        if self.p.a_cw == -1:
            self.pert_v_r_sq = np.flipud(self.pert_v_r_sq)
            self.pert_v_phi_sq = -1*np.flipud(self.pert_v_phi_sq)
            self.pert_rho_sq = np.flipud(self.pert_rho_sq)

        # save grids
        self.x_cut = x_cut 
        self.y_cut = y_cut
        self.x_sq = self.p.l * x_cut + self.p.r_planet
        self.y_sq = self.p.l * y_cut

        # plotting (for debugging)
        _, ax = plt.subplots()
        myplot = ax.contourf(self.x_sq, self.y_sq, self.pert_rho_sq, levels=300, vmin=-4, vmax=4, cmap='RdBu')
        plt.colorbar(myplot)
        plt.show()

        _, ax = plt.subplots()
        myplot = ax.contourf(self.x_sq, self.y_sq, self.pert_v_r_sq, levels=300, cmap='RdBu')
        plt.colorbar(myplot)
        plt.show()

        _, ax = plt.subplots()
        myplot = ax.contourf(self.x_sq, self.y_sq, self.pert_v_phi_sq, levels=300, cmap='RdBu')
        plt.colorbar(myplot)
        plt.show()

    """
    def add_to_global_grid(self, global_grid):

        # grab grid object
        g = global_grid

        # interpolate over square box
        interp_v_r = RectBivariateSpline(self.y_sq, self.x_sq, self.pert_v_r_sq)
        interp_v_phi = RectBivariateSpline(self.y_sq, self.x_sq, self.pert_v_phi_sq)
        interp_rho = RectBivariateSpline(self.y_sq, self.x_sq, self.pert_rho_sq)

        # convert global grid to Cartesian coordinates
        global_R, global_PHI = np.meshgrid(g.r, g.phi)
        global_X = global_R * np.cos(global_PHI)
        global_Y = global_R * np.sin(global_PHI)

        # evaluate interpolation functions over global grid
        global_v_r_box = interp_v_r.ev(global_Y, global_X)
        global_v_phi_box = interp_v_phi.ev(global_Y, global_X)
        global_rho_box = interp_rho.ev(global_Y, global_X)

        # construct global grid of zeros
        zeros = np.zeros(global_v_r_box.shape)

        # update global grid of zeros to have ones on coordinates inside linear box

        # multiply each interpolated global grid by zeros/ones matrix

        # add interpolation results to grid object

        # plot for debugging
        _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        myplot = ax.contourf(global_PHI, global_R, global_rho_box, levels=300)
        plt.colorbar(myplot)
        plt.show()
    """

    def cut_box_annulus_segment(self):

        print("Extracting Linear Perturbations in the Vicinity of the Planet")

        # box size (in units of Hill radius), note for conversions that self.p.l = 1 Hill radius in cgs
        box_size = 2*self.p.scale_box

        # linear perturbations read in grid
        x = self.X[0,:]
        y = self.Y[:,0]

        # cut square box in linear regime
        x_cut = x[np.argmin(x < -box_size) : np.argmin(x < box_size) + 1]
        y_cut = y[np.argmin(y < -box_size) : np.argmin(y < box_size) + 1]

        # annulus segment grid, granularity from square box
        r = np.linspace(
            self.p.r_planet - box_size*self.p.l, 
            self.p.r_planet + box_size*self.p.l, 
            len(x_cut)
        )
        phi = np.linspace(
            -self.p.scale_box_ang*np.arcsin(box_size*self.p.l / self.p.r_planet), 
            self.p.scale_box_ang*np.arcsin(box_size*self.p.l / self.p.r_planet), 
            len(y_cut)
        )
        R, PHI = np.meshgrid(r, phi)

        # for the points on our annulus segment, use transformations to find the x,y values on the original perturbations grid
        X_pert_grid = (R * np.cos(PHI) - self.p.r_planet) / self.p.l
        Y_pert_grid = (R * np.sin(PHI)) / self.p.l

        # cut big perturbations grid to just outside annulus
        self.r_min = self.p.r_planet - box_size*self.p.l
        self.r_max = self.p.r_planet + box_size*self.p.l
        self.theta_max =  np.arcsin(box_size*self.p.l / self.p.r_planet)

        # find cut indicies (remember we need to scale back to units of Hill radius )
        x_cut_i1 = np.argmin(x < (self.r_min * np.cos(self.theta_max) - self.p.r_planet) / self.p.l)
        x_cut_i2 = np.argmin(x < self.r_max / self.p.l) + 1
        y_cut_i1 = np.argmin(y < (-1 * self.r_max * np.sin(self.theta_max)) / self.p.l)
        y_cut_i2 = np.argmin(y < (self.r_max * np.sin(self.theta_max)) / self.p.l) + 1

        # cut grid
        x_int_cut = x[x_cut_i1 : x_cut_i2]
        y_int_cut = y[y_cut_i1 : y_cut_i2]

        # cut perturbation arrays
        cut_v_r = self.pert_v_r[y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]
        cut_v_phi = self.pert_v_phi[y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]
        cut_rho = self.pert_rho[y_cut_i1:y_cut_i2, x_cut_i1:x_cut_i2]

        # interpolation over cut (cartesian) grid
        interp_v_r = RectBivariateSpline(y_int_cut, x_int_cut, cut_v_r)
        interp_v_phi = RectBivariateSpline(y_int_cut, x_int_cut, cut_v_phi)
        interp_v_rho = RectBivariateSpline(y_int_cut, x_int_cut, cut_rho)

        # evaluate interpolation over our annulus segment
        self.pert_v_r_ann = interp_v_r.ev(Y_pert_grid, X_pert_grid)
        self.pert_v_phi_ann = interp_v_phi.ev(Y_pert_grid, X_pert_grid)
        self.pert_rho_ann = interp_v_rho.ev(Y_pert_grid, X_pert_grid)

        # scale to cgs units (and account for rotation direction)
        self.pert_v_r_ann *= self.p.c_s_planet*(self.p.m_thermal/self.p.m_star)
        self.pert_v_phi_ann *= self.p.c_s_planet*(self.p.m_thermal/self.p.m_star)*self.p.a_cw
        self.pert_rho_ann *= (self.p.m_planet/self.p.m_thermal)

        # save annulus grid
        self.r_ann = r
        self.phi_ann = phi
        self.R_ann = R
        self.PHI_ann = PHI

        """
        # plotting (for debugging)
        _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        myplot = ax.contourf(PHI, R, self.pert_rho_ann, levels=300, vmin=0.01, vmax=4)
        ax.set_ylim(0, self.p.r_outer)
        plt.colorbar(myplot)
        plt.show()
        """

