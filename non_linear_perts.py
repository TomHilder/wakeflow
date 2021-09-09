import math
import numpy as np
from scipy.interpolate import griddata
from copy import copy
import matplotlib.pyplot as plt
from burgers import solve_burgers
from transformations import phi_wake, Eta, mod2pi, t, t_integral, t_integrand, g, Lambda_fu, Lambda_fv, get_chi, get_dens_vel, xy_to_etat, plot_r_t

# TODO: Still using Francesco's code without modifications besides variable names and arguments, try to optimise/clean up further

class NonLinearPerts():
    def __init__(self, parameters, Grid):
        
        # grab parameters object
        self.p = parameters

        # should be handed an empty grid with the correct dimensions and grid setup used in the run
        self.g = Grid

    def extract_ICs(self, LinearPerts):

        print('  * Extracting Burgers initial condition from linear density perturbation ...')
        
        # grab linear perturbations object
        lp = LinearPerts

        # grid
        x = lp.X[0,:]
        y = lp.Y[:,0]

        # set maximum eta - semi-width of the support of the azimuthal profile of the linear density perturbation
        eta_max = 25 

        # find the index in the x grid corresponding to the edge of the box
        index = np.argmin(np.abs(x - lp.x_box))

        # extract profile of constant x along edge of box for initial condition
        profile = lp.pert_rho[:,index] / np.sqrt(np.abs(lp.x_box))

        # restrict y range --- I'm not sure why this is necessary but it is. Larger values will cause issues
        y_max = 30
        profile_rest = profile[(y > -y_max) & (y < y_max)]
        y_rest = y[(y > -y_max) & (y < y_max)]

        # ## find eta points for IC profile using full transformation
        # find local cart. coords in real units
        x_IC = self.p.l * np.repeat(lp.x_box, len(y_rest))
        y_IC = self.p.l * y_rest

        # find corresponding global polar coords
        r_IC = x_IC + self.p.r_planet
        phi_IC = y_IC / self.p.r_planet

        # initialise arrays for corresponding t,eta points
        eta_IC = np.zeros(len(y_rest))
        t_IC = np.zeros(len(y_rest))

        # perform transformation
        for i in range(len(y_rest)):
            eta_IC[i] = Eta(r_IC[i], phi_IC[i], self.p.r_planet, self.p.hr_planet, self.p.q, -1)
            t_IC[i] = t(r_IC[i], self.p.r_planet, self.p.hr_planet, self.p.q, self.p.p)

        # restrict eta range using eta_max
        self.eta = eta_IC[(eta_IC > -eta_max) & (eta_IC < eta_max)]
        self.profile = profile_rest[(eta_IC > -eta_max) & (eta_IC < eta_max)]

        # set t0
        self.t0 = t_IC[0]

        # old, approximate IC extraction procedure
        if False:
            y_match = -np.sign(lp.x_box) * 0.5 * lp.x_box**2
            y_cut = y[(y - y_match > -eta_max) & (y - y_match < eta_max)]
            self.profile = profile[(y - y_match > -eta_max) & (y - y_match < eta_max)]
            self.eta = y_cut - y_match*np.ones(len(y_cut))
            self.t0 = t(self.p.r_planet + self.p.l * lp.x_box, self.p.r_planet, self.p.hr, self.p.q, self.p.p)

        # set eta_tilde:
        for i in range(len(self.eta)):
            if self.profile[i] == 0 and self.eta[i] > -10 and self.eta[i] < 0:
                zero = self.eta[i]
            elif i!= (len(self.eta) - 1) and self.profile[i] * self.profile[i + 1] < 0 and self.eta[i] > -10 and self.eta[i] < 0:
                zero = 0.5 * (self.eta[i] + self.eta[i + 1])
        self.eta_tilde = -zero

        # set C:
        deta = self.eta[1] - self.eta[0]
        profile0 = self.profile[self.eta < -self.eta_tilde]
        C0 = -np.trapz(profile0, dx = deta)
        self.C = (self.p.gamma + 1) * (self.p.m_planet / self.p.m_thermal) * C0 / 2**(3/4)

        print('     eta_tilde = ', self.eta_tilde)
        print('     C0 = ', C0)
        print('     t0 = ', self.t0)

        beta_p = self.p.m_planet / self.p.m_thermal

        # ======================================
        # === put linear solution into t,eta ===

        # you will need to run cut_box_square on the linear perts object before doing this.
        # it is okay to run it as well as the annulus cut, they don't overwrite each other

        if self.p.show_teta_debug_plots:

            # run square box cut
            lp.cut_box_square()

            # plot r, t
            plot_r_t(self.p)

            # Local Cartesian grid which the linear solution was calculated on, meshgrid version
            X, Y = np.meshgrid(lp.x_cut,lp.y_cut)

            # grab positive x side of linear solution
            x_len = len(lp.x_cut) // 2

            # restrict box to area considered
            X_rest = X[:,x_len:index]
            Y_rest = Y[:,x_len:index]

            # Find Chi for linear solution, only taking positive values of x
            # The "index" value also makes it so we do not go outside the linear box
            linear_profile = lp.cut_rho[:,x_len:index] / np.sqrt(np.abs(X_rest))

            # Find etas for the box using IC method
            linear_eta = Y_rest + np.sign(X_rest) * 0.5 * X_rest**2

            # Find etas using proper transformations
            #hp = self.p.hr * self.p.r_planet
            hp = self.p.hr_planet * self.p.r_planet
            x_glob = 2 * hp * X_rest / 3.
            y_glob = 2 * hp * Y_rest / 3.
            r_glob = x_glob + self.p.r_planet
            phi_glob = y_glob / self.p.r_planet
            eta_lin = np.zeros(linear_profile.shape)
            t_lin = np.zeros(linear_profile.shape)
            for i in range(eta_lin.shape[0]):
                #print(str(i))
                for j in range(eta_lin.shape[1]):
                    eta_lin[i,j] = Eta(r_glob[i,j], phi_glob[i,j], self.p.r_planet, self.p.hr_planet, self.p.q, -self.p.a_cw)
                    t_lin[i,j] = t(r_glob[i,j], self.p.r_planet, self.p.hr_planet, self.p.q, self.p.p)

            # Plot Chi vs eta when using eta transformation used for IC cut out
            """
            for i in range(0, len(linear_eta[0,:]), 10):
                plt.plot(linear_eta[:,i], linear_profile[:,i])
                plt.plot(eta_lin[:,i], linear_profile[:,i])
                plt.show()
            """
            plt.plot(linear_eta[:,-1], linear_profile[:,-1], label="Approximate $\eta$ transformation")
            plt.plot(eta_lin[:,-1], linear_profile[:,-1], label="Full $\eta$ transformation")
            plt.plot(self.eta, self.profile, ls='--', label="Actual IC used")
            plt.legend(loc="lower left")
            plt.xlim(-10,10)
            plt.xlabel('$\eta$')
            plt.ylabel('$\chi$')
            plt.savefig('ICs.pdf')
            plt.show()

            plt.scatter(t_lin, eta_lin)
            plt.show()

            #eta_lin = linear_eta    ### TURN THIS ON AND OFF TO CHANGE ETA TRANSFORMATION --> COMMENT OUT TO USE INTEGRAL TRANSFORM

            # get grid regular in (t,eta) to interpolate onto
            dt_lin = np.diff(np.sort(t_lin.flatten())).mean()       # find time step interval to use on regular grid from t_lin
                                                                    # takes minimum separation between points in original t_lin grid
            t_lin_reg = np.arange(0, self.t0, dt_lin)
            eta_lin_reg = copy(self.eta)
            T_lin_reg, ETA_lin_reg = np.meshgrid(t_lin_reg, eta_lin_reg)

            # interpolate solution over irregular grid onto regular grid
            print("Interpolating Linear solution into (t, eta) space")
            lin_solution = griddata(
                (t_lin.flatten(), eta_lin.flatten()), 
                linear_profile.flatten(),
                (T_lin_reg, ETA_lin_reg),
                method='linear'
            )

            #plt.plot(linear_eta[:,-1], linear_profile[:,-1])       # approx
            plt.plot(eta_lin[:,-1], linear_profile[:,-1])           # integral trans
            plt.plot(eta_lin_reg, lin_solution[:,-1])               # interpolated integral trans
            #plt.plot(self.eta, self.profile, c='k')                # IC
            plt.show()

            # plotting (for debugging)
            _, ax = plt.subplots()
            myplot = ax.contourf(t_lin_reg, eta_lin_reg, np.nan_to_num(lin_solution), levels=np.arange(-4,4,0.1), cmap='RdBu')
            plt.colorbar(myplot)
            plt.show()

            # stick this into the solutions array, update later code so that solutions array builds on this one
            self.linear_solution = np.nan_to_num(lin_solution)
            self.linear_t = t_lin_reg

        else:
            
            self.linear_solution = 0
            self.linear_t = 0


    def get_non_linear_perts(self):

        beta_p = self.p.m_planet / self.p.m_thermal

        print('  * Solving Burgers equation ...')
        time, eta, solution, eta_inner, solution_inner  = solve_burgers(
            self.eta, self.profile, self.p.gamma, beta_p, self.C, self.p.CFL, self.eta_tilde, self.t0, self.linear_solution, self.linear_t, self.p.show_teta_debug_plots
        )

        print('  * Computing nonlinear perturbations ...')

        tf = time[-1]

        eta_tilde = self.eta_tilde
        C = self.C 
        t0 = self.t0
        Rp = self.p.r_planet
        x_match = 2*self.p.scale_box 
        l = self.p.l
        cw = -self.p.a_cw
        hr = self.p.hr_planet
        q = self.p.q

        gamma = self.p.gamma
        csp = self.p.c_s_planet
        p = self.p.p

        if self.g.info["Type"] == "cartesian":

            x = self.g.x 
            y = self.g.y

            dnl = np.zeros((len(x),len(y)))
            unl = np.zeros((len(x),len(y)))
            vnl = np.zeros((len(x),len(y)))

            for i in range(len(x)):
                for j in range(len(y)):
                    xx = x[i]
                    yy = y[j]
                    rr = np.sqrt(xx**2 + yy**2)
                    pphi = np.arctan2(yy,xx)

                    Chi = get_chi(pphi, rr, time, eta, eta_inner, eta_tilde, C, solution, solution_inner, t0, tf, Rp, x_match, l, cw, hr, q, p)
                    dnl[i,j], unl[i,j], vnl[i,j] = get_dens_vel(rr, Chi, gamma, Rp, cw, csp, hr, q, p) # COMPUTE DENSITY AND VELOCITY PERTURBATIONS

        else:

            dnl = np.zeros((self.p.n_phi,self.p.n_r))
            unl = np.zeros((self.p.n_phi,self.p.n_r))
            vnl = np.zeros((self.p.n_phi,self.p.n_r))

            r = self.g.r
            phi = self.g.phi

            for i in range(self.p.n_r):
                rr = r[i]
                for j in range(self.p.n_phi):
                    pphi = phi[j]

                    Chi = get_chi(pphi, rr, time, eta, eta_inner, eta_tilde, C, solution, solution_inner, t0, tf, Rp, x_match, l, cw, hr, q, p)
                    dnl[j,i], unl[j,i], vnl[j,i] = get_dens_vel(rr, Chi, gamma, Rp, cw, csp, hr, q, p) # COMPUTE DENSITY AND VELOCITY PERTURBATIONS

        self.rho = dnl
        self.vr = 1e-5*unl
        self.vphi = 1e-5*vnl