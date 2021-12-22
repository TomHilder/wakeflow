import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d
import matplotlib.pyplot as plt
import os, sys, subprocess

class PhantomDump:
    def __init__(self, parameters=None, Grid=None, vr='phantom/hd163/vr.pix', vphi='phantom/hd163/vphi.pix', rho='phantom/hd163/rho.pix'):

        # grab parameters object
        self.p = parameters
        self.radius = self.p.r_outer

        # should be handed an empty grid with the correct dimensions and grid setup used in the run
        self.g = Grid

        # read in midplane arrays
        self.vr_xy = np.loadtxt(vr)
        self.vphi_xy = np.loadtxt(vphi)
        self.rho_xy = np.loadtxt(rho)

        # get phantom grid
        length = self.vr_xy.shape[0]
        self.x_ph = np.linspace(-self.radius, self.radius, length)
        self.y_ph = np.linspace(-self.radius, self.radius, length)
        self.X_ph, self.Y_ph = np.meshgrid(self.x_ph, self.y_ph)

        # test plot
        plt.contourf(self.X_ph, self.Y_ph, self.rho_xy, levels=200)
        plt.colorbar()
        plt.axis("scaled")
        plt.show()

        plt.contourf(self.X_ph, self.Y_ph, self.vphi_xy, levels=200)
        plt.colorbar()
        plt.axis("scaled")
        plt.show()

        plt.contourf(self.X_ph, self.Y_ph, self.vr_xy, levels=200)
        plt.colorbar()
        plt.axis("scaled")
        plt.show()

    def extract_dens_perts(self, width=1):

        R = np.sqrt(self.X_ph**2 + self.Y_ph**2)
        PHI = np.arctan2(self.Y_ph, self.X_ph)

        # calculate the mean with a lambda
        w = width / 2
        f = lambda r : self.rho_xy[(R >= r-w) & (R < r+w) & (PHI < 7*np.pi/4) & (PHI > np.pi/4)].mean()
        r  = np.linspace(R.min(), R.max(), num=int(self.radius))
        mean = np.vectorize(f)(r)
        mean_interp_func = interp1d(r, mean, kind='linear')

        # test plot
        if True:
            fig,ax=plt.subplots()
            ax.plot(r,mean)
            ax.set_yscale('log')
            plt.show()

        # subtract azimuthal average
        print("Subtracting azimuthal average of density profile")
        averaged = mean_interp_func(R)
        self.rho_xy_pert = (self.rho_xy - averaged) / averaged
        print("dens MIN MAX = ", self.rho_xy_pert.min(), self.rho_xy_pert.max())

        # calculate the mean with a lambda
        w = width / 2
        f = lambda r : self.vphi_xy[(R >= r-w) & (R < r+w) & (PHI < 7*np.pi/4) & (PHI > np.pi/4)].mean()
        r  = np.linspace(R.min(), R.max(), num=int(self.radius))
        mean = np.vectorize(f)(r)
        mean_interp_func = interp1d(r, mean, kind='linear')

        # test plot
        if True:
            fig,ax=plt.subplots()
            ax.plot(r,mean)
            ax.set_yscale('log')
            plt.show()

        # subtract azimuthal average
        print("Subtracting azimuthal average of density profile")
        averaged = mean_interp_func(R)
        self.vphi_xy_pert = (self.vphi_xy - averaged) / averaged
        print("vphi MIN MAX = ", self.vphi_xy_pert.min(), self.vphi_xy_pert.max())

    def get_polar_grid(self):

        # get polar grid of grid object
        self.R = self.g.R[:, 0, :]
        self.PHI = self.g.PHI[:, 0, :]

        # find equivalent Cartesian points
        X = self.R * np.cos(self.PHI)
        Y = self.R * np.sin(self.PHI)

        # interpolate over phantom midplane
        interp_v_r = RectBivariateSpline(self.y_ph, self.x_ph, self.vr_xy)
        interp_v_phi = RectBivariateSpline(self.y_ph, self.x_ph, self.vphi_xy)
        interp_v_rho = RectBivariateSpline(self.y_ph, self.x_ph, self.rho_xy)

        # evaluate interpolation on Grid object grid
        self.vr = interp_v_r.ev(Y, X)
        self.vphi = interp_v_phi.ev(Y, X)
        self.rho = interp_v_rho.ev(Y, X)

        # testing plots
        """
        plt.imshow(self.vr)
        plt.show()

        # test plot
        _, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(8,8))
        myplot = ax.contourf(self.PHI, self.R, self.vr, levels=300)
        plt.colorbar(myplot)
        plt.show()
        """

