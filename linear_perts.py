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

        # linear perturbations read in grid
        x = self.X[0,:]
        y = self.Y[:,0]

        # define constants for linear perts
        self.x_box = 2*self.p.scale_box
        
        # cut square box grid in linear regime
        self.x_cut = x[np.argmin(x < -self.x_box) : np.argmin(x < self.x_box) + 1]
        self.y_cut = y[np.argmin(y < -self.x_box) : np.argmin(y < self.x_box) + 1]
        

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
        X_pert_grid = (R - self.p.r_planet) / self.p.l
        Y_pert_grid = (self.p.r_planet * PHI) / self.p.l


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

        # account for rotation direction
        if self.p.a_cw == -1:
            cut_v_r = np.flipud(cut_v_r)
            cut_v_phi = -1*np.flipud(cut_v_phi)
            cut_rho = np.flipud(cut_rho)

        # interpolation over cut (cartesian) grid
        interp_v_r = RectBivariateSpline(y_int_cut, x_int_cut, cut_v_r)
        interp_v_phi = RectBivariateSpline(y_int_cut, x_int_cut, cut_v_phi)
        interp_v_rho = RectBivariateSpline(y_int_cut, x_int_cut, cut_rho)

        # evaluate interpolation over our annulus segment
        self.pert_v_r_ann = interp_v_r.ev(Y_pert_grid, X_pert_grid)
        self.pert_v_phi_ann = interp_v_phi.ev(Y_pert_grid, X_pert_grid)
        self.pert_rho_ann = interp_v_rho.ev(Y_pert_grid, X_pert_grid)

        # scale to cgs units
        self.pert_v_r_ann *= self.p.c_s_planet*(self.p.m_planet/self.p.m_thermal)
        self.pert_v_phi_ann *= self.p.c_s_planet*(self.p.m_planet/self.p.m_thermal)
        self.pert_rho_ann *= (self.p.m_planet/self.p.m_thermal)

        # scale velocities to km/s
        self.pert_v_r_ann *= 1e-5
        self.pert_v_phi_ann *= 1e-5

        # save annulus grid
        self.r_ann = r
        self.phi_ann = phi
        self.R_ann = R
        self.PHI_ann = PHI

        """
        # plotting (for debugging)
        _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        myplot = ax.contourf(PHI, R, self.pert_v_r_ann, levels=300, vmin=-2E4, vmax=2E4, cmap='RdBu')
        ax.set_ylim(0, self.p.r_outer)
        plt.colorbar(myplot)
        plt.show()
        """

