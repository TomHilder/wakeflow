import numpy as np
from burgers import solve_burgers
from transformations import phi_wake, Eta, mod2pi, t, t_integral, t_integrand, g, Lambda_fu, Lambda_fv, get_chi, get_dens_vel

# TODO: Still using Francesco's code without modifications besides variable names and arguments, try to optimise/clean up further

class NonLinearPerts():
    def __init__(self, parameters, Grid):
        
        # grab parameters object
        self.p = parameters

        # should be handed an empty grid with the correct dimensions and grid setup used in the run
        self.g = Grid

    def extract_burgers_ICs_sq(self, LinearPerts):

        print('  * Extracting Burgers initial condition from linear density perturbation ...')
        
        # grab linear perturbations object
        lp = LinearPerts

        # grid
        x = lp.X[0,:]
        y = lp.Y[:,0]

        # set maximum eta - semi-width of the support of the azimuthal profile of the linear density perturbation
        eta_max = 25 

        index = np.argmin(np.abs(x - lp.x_box))
        profile = lp.pert_rho[:,index] / np.sqrt(np.abs(lp.x_box))
        y_match = -np.sign(lp.x_box) * 0.5 * lp.x_box**2
        y_cut = y[(y - y_match > -eta_max) & (y - y_match < eta_max)]
        self.eta = y_cut - y_match*np.ones(len(y_cut))
        self.profile = profile[(y - y_match > -eta_max) & (y - y_match < eta_max)]

        # set t0 (t correponding to the x_match)
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

    def extract_burgers_ICs_annulus_segment(self, LinearPerts):
        pass

    def get_non_linear_perts(self):

        beta_p = self.p.m_planet / self.p.m_thermal

        print('  * Solving Burgers equation ...')
        time, eta, solution, eta_inner, solution_inner  = solve_burgers(self.eta, self.profile, self.p.gamma, beta_p, self.C, self.p.CFL, self.eta_tilde)

        print('  * Computing nonlinear perturbations ...')

        tf = time[-1]

        dnl = np.zeros((self.p.n_phi,self.p.n_r))
        unl = np.zeros((self.p.n_phi,self.p.n_r))
        vnl = np.zeros((self.p.n_phi,self.p.n_r))

        r = self.g.r
        phi = self.g.phi

        eta_tilde = self.eta_tilde
        C = self.C 
        t0 = self.t0
        Rp = self.p.r_planet
        x_match = 2*self.p.scale_box 
        l = self.p.l
        cw = -self.p.a_cw
        hr = self.p.hr
        q = self.p.q

        gamma = self.p.gamma
        csp = self.p.c_s_planet
        p = self.p.p

        for i in range(self.p.n_r):
            rr = r[i]
            for j in range(self.p.n_phi):
                pphi = phi[j]

                Chi = get_chi(pphi, rr, time, eta, eta_inner, eta_tilde, C, solution, solution_inner, t0, tf, Rp, x_match, l, cw, hr, q, p)
                dnl[j,i], unl[j,i], vnl[j,i] = get_dens_vel(rr, Chi, gamma, Rp, cw, csp, hr, q, p) # COMPUTE DENSITY AND VELOCITY PERTURBATIONS

        self.rho = dnl
        self.vr = unl
        self.vphi = vnl

        print("-- done --")