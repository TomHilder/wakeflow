# grid.py
# Written by Thomas Hilder

"""
Contains the Grid class on which all Wakeflow models are run/stored.
"""

import numpy                as np
import matplotlib.pyplot    as plt
import astropy.io.fits      as fits
import shutil               as sh
from scipy.interpolate  import RectBivariateSpline
from matplotlib.colors  import LogNorm
from .mcfost_interface  import _read_mcfost_grid_data
from .model_setup       import _Parameters
from .linear_perts      import _LinearPerts
from .non_linear_perts  import _NonLinearPerts

# NOTE: contents are intended mostly for internal use and should not really be accessed by users, but I have still provided
# documentation in the case that advanced users want to mess around with any of it.

# a class to store numpy arrays of results together, so that nothing gets confused. Also stores information and parameters
class _Grid:
    """
    Grid object to generate and store Wakeflow results on. 

    Note that all 3D arrays in the Grid have dimensions (x,z,y) or (phi,z,r).

    Attributes
    ----------
    p : Parameters
        Parameters object (model_setup.py) containing the options used for the model.
    info : dict
        Contains the type, size and contents of the grid.
    height: ndarray
        Scale height of the disk.
    v_r : ndarray
        Contains radial velocities.
    v_phi : ndarray
        Contains azimuthal velocities, positive is counterclockwise.
    rho : ndarray
        Contains densities.
    x : ndarray
        1D. Contains x coordinates of Cartesian grid.
    y : ndarray
        1D. Contains y coordinates of Cartesian grid.
    z_xy : ndarray
        1D. Contains z coordinates of Cartesian grid.
    X : ndarray
        3D. Contains x coordinates of Cartesian grid points.
    Y : ndarray
        3D. Contains y coordinates of Cartesian grid points.
    Z_xy : ndarray
        3D. Contains z coordinates of Cartesian grid points.
    R_xy : ndarray
        3D. Contains radial coordinates of Cartesian grid points.
    PHI_xy : ndarray
        3D. Contains azimuthal coordinates of Cartesian grid points.
    r : ndarray
        1D. Contains radial coordinates of Cylindrical or Mcfost grid.
    phi : ndarray
        1D. Contains azimuthal coordinates of Cylindrical or Mcfost grid.
    z : ndarray
        1D. Contains z coordinates of Cylindrical grid.
    R : ndarray
        3D. Contains radial coordinates of Cylindrical or Mcfost grid points.
    PHI : ndarray
        3D. Contains azimuthal coordinates of Cylindrical or Mcfost grid points.
    Z : ndarray
        3D. Contains z coordinates of Cylindrical or Mcfost grid points.
    """

    # instantiate the class
    def __init__(self, parameters : _Parameters) -> None:
        """
        Parameters
        ----------
        parameters : _Parameters
            _Parameters object (model_setup.py) containing parameters to be used for the Grid.
        """

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
        self.v_r   = None 
        self.v_phi = None 
        self.rho   = None

    # setup the grids with chosen geometry
    def _make_grid(self) -> None:
        """Generates grid according to parameters in self.p and define disk scale height
        """

        # define disk height (not used for mcfost grid)
        self.height = self.p.hr * (self.p.r_outer / self.p.r_ref)**(0.5 - self.p.q) * self.p.r_outer

        # make cartesian grid
        if self.p.grid_type == "cartesian":
            self._make_cartesian_grid()

        # make cylindrical grid
        elif self.p.grid_type == "cylindrical":
            self._make_cylindrical_grid()

        # make mcfost grid
        elif self.p.grid_type == "mcfost":
            self._make_mcfost_grid()

    # setup Cartesian grid
    def _make_cartesian_grid(self) -> None:
        """Generates a Cartesian grid using parameters in self.p
        """
        
        # make grid from specified parameters
        self.x    = np.linspace(-self.p.r_outer, self.p.r_outer, self.p.n_x)
        self.y    = np.linspace(-self.p.r_outer, self.p.r_outer, self.p.n_y)
        self.z_xy = np.linspace(0, self.height, self.p.n_z)

        self.X, self.Z_xy, self.Y  = np.meshgrid(self.x, self.z_xy, self.y, indexing='ij')

        # update grid info
        self.info["Type"]    = "cartesian"
        self.info["Size"][0] = self.x.shape   [0]
        self.info["Size"][1] = self.z_xy.shape[0]
        self.info["Size"][2] = self.y.shape   [0]

    # setup Cylindrical grid
    def _make_cylindrical_grid(self) -> None:
        """Generates a Cylindrical grid using parameters in self.p
        """
         
        # make grid from specifications in parameter file
        if self.p.r_log:
            self.r = np.geomspace(self.p.r_inner, self.p.r_outer, self.p.n_r)
        else:
            self.r = np.linspace (self.p.r_inner, self.p.r_outer, self.p.n_r)    
        #self.phi = np.linspace(0, 2*np.pi, self.p.n_phi)
        self.phi = np.linspace(-np.pi, np.pi, self.p.n_phi)
        self.z   = np.linspace(0, self.height, self.p.n_z)

        self.PHI, self.Z, self.R = np.meshgrid(self.phi, self.z, self.r, indexing='ij')

        # update grid info
        self.info["Type"]    = "cylindrical"
        self.info["Size"][0] = self.phi.shape[0]
        self.info["Size"][1] = self.z.shape  [0]
        self.info["Size"][2] = self.r.shape  [0]
        self.info["Log_r"]   = self.p.r_log

    # setup Mcfost grid by calling Mcfost to generate the geometry, and then read in that geometry
    def _make_mcfost_grid(self) -> None:
        """Generates an Mcfost grid using parameters in self.p
        """

        # generate mcfost grid data from specifications in parameter file
        self.R, self.Z = _read_mcfost_grid_data(self.p)
        self.r         = self.R[0, 0, :]
        self.phi       = np.linspace(0, 2*np.pi, self.p.n_phi)

        PHI_temp = np.repeat(self.phi[:, np.newaxis],    self.p.n_r, axis=1)
        self.PHI = np.repeat(PHI_temp[:, np.newaxis, :], self.p.n_z, axis=1)

        # note that we do not have a self.z because it depends on radius, as the grid is flared

        # update grid info
        self.info["Type"]    = "mcfost"
        self.info["Size"][0] = self.phi.shape[0]
        self.info["Size"][1] = self.Z.shape  [1]
        self.info["Size"][2] = self.r.shape  [0]
        self.info["Log_r"]   = "mcfost"

    # we need (r,phi) points for the Cartesian grid points in order to do transformations
    def _get_r_phi_coords(self) -> None:
        """Find the (r,phi) coordinates corresponding to each (x,y) point
        """

        if self.info["Type"] != "cartesian":
            raise Exception("Can only be used on Cartesian grid.")

        # need the equivalent (r,phi) coordinates when using cartesian grid
        self.R_xy   = np.sqrt   (self.X**2 + self.Y**2)
        self.PHI_xy = np.arctan2(self.Y,     self.X)

    # if the user desires dimensionless results, needs to be used last after all grids are merged
    def _remove_dimensions(self, scale_dens : bool = False) -> None:
        """Convert results to be dimensionless by appropriate scaling; should be used last

        Parameters
        ----------
        scale_dens : bool, optional
            If true, the density is scaled by the rho_ref value in the parameters. Should be left as default for
            scaling perturbations (delta rho / rho) and set to True otherwise.
        """

        grid_length = self.p.U_len / self.p.au

        if self.p.grid_type == "cartesian":
            self.x      /= grid_length
            self.y      /= grid_length
            self.z_xy   /= grid_length
            self.X      /= grid_length
            self.Y      /= grid_length
            self.Z_xy   /= grid_length

        elif self.p.grid_type == "cylindrical":
            self.r      /= grid_length
            self.phi    /= 1.0
            self.z      /= grid_length
            self.R      /= grid_length
            self.PHI    /= 1.0
            self.Z      /= grid_length

        else:
            self.r      /= grid_length
            self.phi    /= 1.0
            self.R      /= grid_length
            self.PHI    /= 1.0
            self.Z      /= grid_length
            
        self.v_r   /= self.p.U_vel
        self.v_phi /= self.p.U_vel

        if scale_dens:
            self.rho /= self.p.rho_ref

    # results for clockwise rotation are the same under a reflection about the y-axis, or phi-axis, 
    # as well as a reverse of sign in v_phi
    def _flip_results(self) -> None:
        """Go from counterclockwise to clockwise rotation; should be used last
        """

        if self.p.grid_type == "cartesian":
            self.y *= -1
            self.Y *= -1
        
        else:
            self.phi *= -1
            self.PHI *= -1

        self.v_phi *= -1

    # setup the unperturbed background disk onto which the perturbations will be mapped
    def _make_keplerian_disk(self) -> None:
        """Generate unperturbed accretion disk in vertical hydrostatic equilibrium, according to parameters in self.p
        """

        # get radii
        if self.info["Type"] == "cartesian":
            self._get_r_phi_coords()
            r = self.R_xy
            z = self.Z_xy
        else:
            r = self.R
            z = self.Z

        # rename for readability in long equations
        p       = self.p.dens_p
        q       = self.p.q
        hr      = self.p.hr
        r_ref   = self.p.r_ref
        r_c     = self.p.r_c
        rho_ref = self.p.rho_ref
        c_s_0   = self.p.c_s_0
        G       = self.p.G_const
        m_star  = self.p.m_star
        m_sol   = self.p.m_solar
        au      = self.p.au

        # Keplerian velocities
        self.v_r   = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))
        self.v_kep = 1e-5 * np.sqrt(G * m_star * m_sol / (r * au))
        self.v_phi = np.copy(self.v_kep) * self.p.a_cw

        # get h/r(r)
        hrf = hr * (r / r_ref) ** (0.5 - q)

        # pressure and gravity height dependent correction
        corr = np.sqrt(
            (-1* (p + 2 * q) * hrf**2) + (1 - 2 * q) + (2 * q * r / np.sqrt(r**2 + z**2))
        )

        # perform correction
        self.v_phi *= corr

        # sound speed profile
        c_s = c_s_0 * (r / r_ref)**-q

        # unperturbed density profile (full) (no taper)
        if r_c == 0: 
            self.rho = (
                rho_ref * (r / r_ref)**-p * np.exp(
                    (G * m_star * m_sol / c_s**2) * (1 / np.sqrt((r * au)**2 + (z * au)**2) - 1 / (r * au))
                )
            )
        else: # with taper
            self.rho = (
                rho_ref * (r / r_ref)**-p * np.exp(-(r / r_c)**(2 - p)) * np.exp(
                    (G * m_star * m_sol / c_s**2) * (1 / np.sqrt((r * au)**2 + (z * au)**2) - 1 / (r * au))
                )
            )

        # unperturbed density profile (thin disk approx.)
        #self.rho = rho_ref * (r / r_ref)**-p * np.exp( - 0.5 * (z / r_ref)**2 * hr**-2 * (r / r_ref)**(2 * q - 3))

        # update grid info
        self.info["Contains"] = "Keplerian velocity field, unperturbed density"

    # fill the grids with zeros, used for storing perturbations on separate grids before mapping onto background
    def _make_empty_disk(self) -> None:
        """fill v_r, v_phi and rho with zeros
        """

        # make grids full of zeros
        self.v_r   = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))
        self.v_phi = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))
        self.rho   = np.zeros((self.info["Size"][0], self.info["Size"][1], self.info["Size"][2]))

        # update grid info
        self.info["Contains"] = "Zeros"

    # add results from the linear regime nearby the planet onto the grid
    def _add_linear_perturbations(self, LinearPerts : _LinearPerts, rho_background : "_Grid.rho") -> None:
        """Add results from the linear regime, stored in LinearPerts object, to the grid

        Parameters
        ----------
        LinearPerts : LinearPerts
            LinearPerts object containing results from the linear regime.
        rho_background : Grid.rho
            unperturbed density from Grid object where make_keplerian_disk has been used.
        """

        # box size (in units of Hill radius), note for conversions that self.p.l = 1 Hill radius in cgs
        x_box_size = 2 * self.p.scale_box
        y_box_size = 2 * self.p.scale_box_ang

        min_r = self.p.r_planet - x_box_size * self.p.l
        max_r = self.p.r_planet + x_box_size * self.p.l

        min_y = -y_box_size * self.p.l
        max_y =  y_box_size * self.p.l

        max_phi =  np.pi / 2
        min_phi = -np.pi / 2

        # find (phi, r) grid for either Cartesian or Cylindrical global grid
        if self.info["Type"] == "cartesian":
            self._get_r_phi_coords()
            R, PHI = self.R_xy, self.PHI_xy
        else:
            R, PHI = self.R, self.PHI

        # new PHI grid to use (-pi,pi) instead of (0,2pi), where values are swapped in place, also ditch z coordinate
        # also construct a mask that contains 0 outside linear annulus and 1 inside
        PHI_new     = np.zeros((PHI.shape[0],PHI.shape[2]))
        Y_new       = np.zeros((PHI.shape[0],PHI.shape[2]))
        linear_mask = np.zeros((PHI.shape[0],PHI.shape[2]))

        R_new = R[:,0,:]

        for i in range(PHI.shape[0]):
            for j in range(PHI.shape[2]):

                # transforming phi coordinate in place
                if PHI[i,0,j] > np.pi:
                    PHI_new[i,j] = PHI[i,0,j] - 2*np.pi
                else:
                    PHI_new[i,j] = PHI[i,0,j]

                Y_new[i,j] = R[i,0,j] * np.sin(PHI[i,0,j])

                # constructing mask
                if PHI_new[i,j] > min_phi and PHI_new[i,j] < max_phi \
                    and Y_new[i,j] > min_y and Y_new[i,j] < max_y \
                    and R_new[i,j] > min_r and R_new[i,j] < max_r:
                    linear_mask[i,j] = 1
                else:
                    linear_mask[i,j] = 0

        # get linear solution           
        lp = LinearPerts

        # assemble interpolation functions over linear perts grid
        interp_v_r   = RectBivariateSpline(lp.y_ann, lp.r_ann, lp.pert_v_r_ann)
        interp_v_phi = RectBivariateSpline(lp.y_ann, lp.r_ann, lp.pert_v_phi_ann)
        interp_rho   = RectBivariateSpline(lp.y_ann, lp.r_ann, lp.pert_rho_ann)

        # evaluate interpolations on global grid
        global_lin_v_r   = interp_v_r.ev  (Y_new, R_new)
        global_lin_v_phi = interp_v_phi.ev(Y_new, R_new)
        global_lin_rho   = interp_rho.ev  (Y_new, R_new)

        # apply mask to only get solution in valid domain
        global_lin_v_r   = global_lin_v_r   * linear_mask
        global_lin_v_phi = global_lin_v_phi * linear_mask
        global_lin_rho   = global_lin_rho   * linear_mask

        # Add velocities, identical at all heights for now
        self.v_r   += global_lin_v_r  [:, np.newaxis, :]
        self.v_phi += global_lin_v_phi[:, np.newaxis, :]

        # Add density, scaling by background density
        self.rho += global_lin_rho[:, np.newaxis, :] * rho_background

        # plot for debugging
        #_, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        #myplot = ax.contourf(PHI[:,0,:], R[:,0,:], global_lin_v_r, levels=300, cmap='RdBu')
        #ax.set_ylim(0, self.p.r_outer)
        #plt.colorbar(myplot)
        #plt.show()

        # update grid info
        self.info["Contains"] = "linear perturbations"

    # add results from non-linear wake propagation onto the grid
    def _add_non_linear_perturbations(self, NonLinearPerts : _NonLinearPerts, rho_background : "_Grid.rho") -> None:
        """Add results from the non-linear regime, stored in NonLinearPerts object, to the grid

        Parameters
        ----------
        NonLinearPerts : NonLinearPerts
            NonLinearPerts object containing results from the non-linear regime.
        rho_background : Grid.rho
            unperturbed density from Grid object where make_keplerian_disk has been used.
        """

        nonlin = NonLinearPerts

        # Add velocities, identical at all heights for now
        nl_vr   = nonlin.vr  [:, np.newaxis, :]
        nl_vphi = nonlin.vphi[:, np.newaxis, :]

        self.v_r   += nl_vr
        self.v_phi += nl_vphi

        # Define scale height array
        if self.info["Type"] == "cartesian":

            self._get_r_phi_coords()
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

        self.rho += nonlin.rho[:, np.newaxis, :] * rho_background

        # update grid info
        self.info["Contains"] = "non-linear perturbations"

#    def add_phantom_midplane(self, PhantomDump):
#
#        PD = PhantomDump
#
#        # Add velocities, identical at all heights
#        ph_vr   = PD.vr  [:, np.newaxis, :]
#        ph_vphi = PD.vphi[:, np.newaxis, :]
#
#        self.v_r   += ph_vr #/(1 + self.Z)
#        self.v_phi += ph_vphi #/(1 + self.Z)
#
#        # Define scale height array
#        self.H = self.p.h_ref * (self.R / self.p.r_ref)**self.p.beta
#
#        # Add density, scaling vertically as per thin disk assumption
#        ph_rho   = PD.rho[:, np.newaxis, :]
#        self.rho = ph_rho * np.exp(-0.5 * (self.Z / self.H)**2)
#
#        # update grid info
#        self.info["Contains"] = "phantom mid-plane extrapolated to 3D"

    # take multiple grid objects and merge them if the geometry is the same, simply by adding the components
    def _merge_grids(self, grid_to_merge : "_Grid") -> None:
        """Merge v_r, v_phi and rho from two compatible Grid objects (ie. with the same _Parameters)
        """

        g = grid_to_merge

        if type(grid_to_merge) is not _Grid:
            raise Exception("merge_grids must be given a Grid object. Merging failed.")
        
        if self.p != g.p:
            raise Exception("merge_grids requires both grids must have same parameters. Merging failed.")

        # merge data arrays
        self.v_r   += g.v_r
        self.v_phi += g.v_phi
        self.rho   += g.rho

        # update info
        self.info["Contains"] += " AND " + g.info["Contains"]

#    def merge_phantom_densities(self, grid_to_merge):
#
#        g = grid_to_merge
#
#        if type(grid_to_merge) is not Grid:
#            print("Must be given a Grid object")
#            return False
#        
#        if self.p != g.p:
#            print("Both grids must have same type parameters")
#            return False
#
#        # merge data arrays
#        self.rho = g.rho

    # create plots of a constant z slice, mostly used for debugging but also shows midplane results to user
    # if they set show_midplane_plots=True
    def _show_disk2D(
        self, 
        z_slice :    int, 
        show :      bool = False, 
        save :      bool = False, 
        dimless :   bool = False, 
        vphi_lim : float = 1.0, 
        cw :        bool = False
    ) -> None:
        """Generate plots of the disk at constant z value with index z_slice

        For example, mid-plane plots are generated with z_slice=0.

        Parameters
        ----------
        z_slice : int
            Index in Grid of z_slice to plot.
        show : bool, optional
            True to show plots to user using matplotlib.pyplot.show(). Default is False.
        save : bool, optional
            True to save plots in directory with results. Default is False.
        dimless: bool, optional
            Set to True if using dimensionless results to get accurately labelled plots.
        vphi_lim : float, optional
            Modifier for the limits on the v_phi plot.
        cw : bool, optional
            True for clockwise rotation and False for anticlockwise. Since anticlockwise is default, False is default.
        """

        # parameters
        savedir      = f"{self.p.system}/{self.p.name}/{self.p.m_planet}Mj"
        contour_lvls = 300

        # find maximum contours
        vr_max   = np.max(self.v_r  [:,z_slice,:])
        vphi_max = np.max(self.v_phi[:,z_slice,:]) * vphi_lim
        if cw == True:
            vphi_max = np.min(self.v_phi[:,z_slice,:]) * vphi_lim
        rho_max  = np.percentile(self.rho[:,z_slice,:], 99.9)

        # Cartesian grid plots
        if self.info["Type"] == "cartesian":

            # v_r
            plt.close("all")
            fig, ax = plt.subplots(dpi=150)
            c       = ax.pcolormesh(self.x, self.y, np.transpose(self.v_r[:,z_slice,:]), vmin=-vr_max, vmax=vr_max, cmap='RdBu', rasterized=True)
            ax.axis('scaled')
            ax.set_title(r"$\delta v_r$")
            if not dimless:
                ax.set_xlabel(r"$x \, [\mathrm{AU}]$")
                ax.set_ylabel(r"$y \, [\mathrm{AU}]$")
                fig.colorbar(c, extend="both", label="$[\mathrm{km / s}]$")
            else:
                ax.set_xlabel(r"$x$")
                ax.set_ylabel(r"$y$")
                fig.colorbar(c, extend="both")
            if save:
                plt.savefig(f'{savedir}/vr_z{z_slice}.pdf')
            if show:
                plt.show()

            # v_phi
            plt.close("all")
            fig, ax = plt.subplots(dpi=150)
            c       = ax.pcolormesh(self.x, self.y, np.transpose(self.v_phi[:,z_slice,:]), vmin=-vphi_max, vmax=vphi_max, cmap='RdBu', rasterized=True)
            ax.axis('scaled')
            ax.set_title(r"$\delta v_{\phi}$")
            if not dimless:
                ax.set_xlabel(r"$x \, [\mathrm{AU}]$")
                ax.set_ylabel(r"$y \, [\mathrm{AU}]$")
                fig.colorbar(c, extend="both", label="$[\mathrm{km / s}]$")
            else:
                ax.set_xlabel(r"$x$")
                ax.set_ylabel(r"$y$")
                fig.colorbar(c, extend="both")
            if save:
                plt.savefig(f'{savedir}/vphi_z{z_slice}.pdf')
            if show:
                plt.show()

            # rho
            plt.close("all")
            fig, ax = plt.subplots(dpi=150)
            c       = ax.pcolormesh(self.x, self.y, np.transpose(self.rho[:,z_slice,:]), vmin=-rho_max, vmax=rho_max, cmap='RdBu', rasterized=True)
            ax.axis('scaled')
            ax.set_title(r"$\delta \rho \, / \rho$")
            if not dimless:
                ax.set_xlabel(r"$x \, [\mathrm{AU}]$")
                ax.set_ylabel(r"$y \, [\mathrm{AU}]$")
            else:
                ax.set_xlabel(r"$x$")
                ax.set_ylabel(r"$y$")
            fig.colorbar(c, extend="both")
            if save:
                plt.savefig(f'{savedir}/rho_z{z_slice}.pdf')
            if show:
                plt.show()

        # Polar grid plots
        else:

            # v_r
            plt.close("all")
            fig, ax = plt.subplots(dpi=150)
            c       = ax.pcolormesh(self.phi, self.r, np.transpose(self.v_r[:,z_slice,:]), vmin=-vr_max, vmax=vr_max, cmap='RdBu', rasterized=True)
            ax.set_title(r"$\delta v_r$")
            if not dimless:
                ax.set_xlabel(r"$\phi \, [\mathrm{rad}]$")
                ax.set_ylabel(r"$r \, \, [\mathrm{AU}]$")
                fig.colorbar(c, extend="both", label="$[\mathrm{km / s}]$")
            else:
                ax.set_xlabel(r"$\phi \, [\mathrm{rad}]$")
                ax.set_ylabel(r"$r$")
                fig.colorbar(c, extend="both")
            if save:
                plt.savefig(f'{savedir}/vr_z{z_slice}.pdf')
            if show:
                plt.show()

            # v_phi
            plt.close("all")
            fig, ax = plt.subplots(dpi=150)
            c       = ax.pcolormesh(self.phi, self.r, np.transpose(self.v_phi[:,z_slice,:]), vmin=-vphi_max, vmax=vphi_max, cmap='RdBu', rasterized=True)
            ax.set_title(r"$\delta v_\phi$")
            if not dimless:
                ax.set_xlabel(r"$\phi \, [\mathrm{rad}]$")
                ax.set_ylabel(r"$r \, \, [\mathrm{AU}]$")
                fig.colorbar(c, extend="both", label="$[\mathrm{km / s}]$")
            else:
                ax.set_xlabel(r"$\phi \, [\mathrm{rad}]$")
                ax.set_ylabel(r"$r$")
                fig.colorbar(c, extend="both")
            if save:
                plt.savefig(f'{savedir}/vphi_z{z_slice}.pdf')
            if show:
                plt.show()

            # rho
            plt.close("all")
            fig, ax = plt.subplots(dpi=150)
            c       = ax.pcolormesh(self.phi, self.r, np.transpose(self.rho[:,z_slice,:]), vmin=-rho_max, vmax=rho_max, cmap='RdBu', rasterized=True)
            ax.set_title(r"$\delta \rho \, / \rho$")
            if not dimless:
                ax.set_xlabel(r"$\phi \, [\mathrm{rad}]$")
                ax.set_ylabel(r"$r \, \, [\mathrm{AU}]$")
            else:
                ax.set_xlabel(r"$\phi \, [\mathrm{rad}]$")
                ax.set_ylabel(r"$r$")
            fig.colorbar(c, extend="both")
            if save:
                plt.savefig(f'{savedir}/rho_z{z_slice}.pdf')
            if show:
                plt.show()
        
        plt.close("all")
    
    # create a fits file from the results formatted specifically for use in MCFOST in combination with
    # the .para file written by pymcfost through wakeflow
    def _write_fits_file(self) -> None:
        """Write a .FITS file from the results compatible with being run in MCFOST
        """

        if self.info["Type"] == "cartesian":

            print("WARNING: Cannot write FITS file when using Cartesian grid.")

        else:

            # create empty array for vertical velocities
            self.v_z = np.zeros((self.p.n_phi, self.p.n_z, self.p.n_r))

            # create master velocities array
            velocities = 1e3 * np.array([self.v_r, self.v_phi, self.v_z])

            # setup HDUs for fits file
            primary_hdu  = fits.PrimaryHDU(np.abs(self.rho))
            second_hdu   = fits.ImageHDU  (np.abs(self.rho)) # * np.exp(-1 * (self.Z / self.H)**2)
            tertiary_hdu = fits.ImageHDU  (velocities)

            # set header properties for mcfost
            primary_hdu.header['hierarch read_gas_velocity'] = 2
            primary_hdu.header['hierarch gas_dust_ratio']    = 100
            primary_hdu.header['hierarch read_gas_density']  = 1
            primary_hdu.header['read_n_a']                   = 0

            # sink particles properties in header (code modifications needed for anything other than 2)
            primary_hdu.header['hierarch N_sink'] = 2

            # star properties in header
            primary_hdu.header['hierarch M_star']  = self.p.m_star       # star mass in M_sol
            primary_hdu.header['hierarch x_star']  = 0.                  # star at origin
            primary_hdu.header['hierarch y_star']  = 0.
            primary_hdu.header['hierarch z_star']  = 0.
            primary_hdu.header['hierarch vx_star'] = 0.                  # star not moving
            primary_hdu.header['hierarch vy_star'] = 0.
            primary_hdu.header['hierarch vz_star'] = 0.

            # planet properties in header
            primary_hdu.header['hierarch M_planet']  = self.p.m_planet    # planet mass in M_jup
            primary_hdu.header['hierarch x_planet']  = self.p.r_planet    # planet x position simply orbital radius [AU]
            primary_hdu.header['hierarch y_planet']  = 0.
            primary_hdu.header['hierarch z_planet']  = 0.
            primary_hdu.header['hierarch vx_planet'] = 0.                # planet velocity just keplerian rotation [m/s] (positive for anticlockwise, which is cw=-1)
            primary_hdu.header['hierarch vy_planet'] = np.sqrt(
                self.p.G_const*self.p.m_star*self.p.m_solar/(self.p.r_planet*self.p.au)
                ) / 100 * self.p.a_cw
            primary_hdu.header['hierarch vz_planet'] = 0.

            # setup fits file
            hdul = fits.HDUList([primary_hdu, second_hdu, tertiary_hdu])

            # Write a fits file for mcfost
            fitsname = "wakeflow_model.fits"
            hdul.writeto(f"{self.p.system}/{self.p.name}/{self.p.m_planet}Mj/{fitsname}", overwrite=True)
            print("Saved FITS file formatted for MCFOST.")

            # Copy mcfost parameter file here too
            sh.copy(
                f"{self.p.system}/{self.p.name}/mcfost/mcfost_{self.p.name}.para",
                f"{self.p.system}/{self.p.name}/{self.p.m_planet}Mj/mcfost.para"
            )
    
    # simply saves the grid and densities/velocities to .npy files so that they may be read later
    def _save_results(self, label : str, printed : str) -> None:
        """ Save results and grid to .npy files in results directory

        Parameters
        ----------
        label : str
            Prepended to saved file names (used in wakeflow.py to give "delta" for perturbations and "total" for perturbations 
            applied on top of an unperturbed disk).
        printed : str
            Message to be printed to the user in format f"{printed} saved to {savedir}".
        """

        savedir = f"{self.p.system}/{self.p.name}/{self.p.m_planet}Mj"

        # save grid:
        if self.p.grid_type == "cartesian":
            np.save(f"{savedir}/X.npy", self.X)
            np.save(f"{savedir}/Z.npy", self.Z_xy)
            np.save(f"{savedir}/Y.npy", self.Y)
        else:
            np.save(f"{savedir}/{label}_PHI.npy", self.PHI)
            np.save(f"{savedir}/{label}_Z.npy", self.Z)
            np.save(f"{savedir}/{label}_R.npy", self.R)

        # save results:
        np.save(f"{savedir}/{label}_v_r.npy", self.v_r)
        np.save(f"{savedir}/{label}_v_phi.npy", self.v_phi)
        np.save(f"{savedir}/{label}_rho.npy", self.rho)

        print(f"{printed} saved to {savedir}")

        