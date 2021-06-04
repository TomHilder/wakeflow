import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
import os, sys, subprocess

class PhantomDump:
    def __init__(self, parameters, Grid, vr='phantom_pixelmaps/hd163/vr.pix', vphi='phantom_pixelmaps/hd163/vphi.pix', rho='phantom_pixelmaps/hd163/rho.pix',):

        # grab parameters object
        self.p = parameters

        # should be handed an empty grid with the correct dimensions and grid setup used in the run
        self.g = Grid

        # read in midplane arrays
        self.vr_xy = np.flipud(np.loadtxt(vr).transpose())
        self.vphi_xy = -1*np.flipud(np.loadtxt(vphi).transpose())
        self.rho_xy = np.flipud(np.loadtxt(rho).transpose())

        # get phantom grid
        length = self.vr_xy.shape[0]
        self.x_ph = np.linspace(-500, 500, length)
        self.y_ph = np.linspace(-500, 500, length)
        self.X_ph, self.Y_ph = np.meshgrid(self.x_ph, self.y_ph)

        # test plot
        """
        plt.figure(figsize=(8,8))
        plt.contourf(self.x_ph, self.y_ph, self.vr_xy.transpose(), levels=1000)
        plt.xlabel('[au]')
        plt.ylabel('[au')
        plt.show()
        """

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

