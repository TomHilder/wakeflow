import numpy                    as np
import matplotlib.pyplot        as plt
from scipy.integrate        import quad
from scipy.interpolate      import RectBivariateSpline

def phi_wake(r, Rp, hr, q, cw): 
    """Eq. (4) Bollati et al. 2021
    """
    rr = r / Rp
    return -cw * np.sign(r - Rp) * (1 / hr) * (rr**(q - 0.5) / (q - 0.5) - rr**(q + 1) / (q + 1) - 3 / ((2 * q - 1) * (q + 1)))

def Eta(r, phi, Rp, hr, q, cw):
    """Eq. (14) Bollati et al. 2021
    """
    coeff    = 1.5 / hr
    phi_w    = mod2pi(phi_wake(r, Rp, hr, q, cw))
    deltaphi = phi - phi_w

    if deltaphi > np.pi:
        deltaphi = deltaphi - 2 * np.pi

    elif deltaphi < -np.pi:
        deltaphi = deltaphi + 2 * np.pi 

    return coeff * deltaphi

def mod2pi(phi):

    if phi >= 0:
        return phi % (2 * np.pi) 

    else:
        if np.abs(phi) < np.pi * 2:
            return phi + 2 * np.pi    

        else:
            phi   = -phi
            resto =  phi % (2 * np.pi)    
            
        return -resto + 2 * np.pi

def t_integrand(x, q, p):
    rho = 5 * q + p
    w   = rho / 2 - 11 / 4
    return np.abs(1 - x**(1.5))**(1.5) * x**w

def t_integral(up, q, p):
    return  quad(t_integrand, 1, up, args=(q,p))[0]

def t(r, Rp, hr, q, p):
    """Equation (43) Rafikov 2002    (Eq. 13 Bollati et al. 2021)
    """
    module_integral = np.abs(t_integral(r / Rp, q, p))
    coeff = 3 * hr**(-5 / 2) / (2**(5 / 4))
    return coeff * module_integral

def g(r, Rp, hr, q, p):
    """Equation (12) Bollati et al. 2021
    """
    coeff = 2**0.25 * hr**0.5
    term1 = (r / Rp)**(0.5 * (1 - p - 3 * q))
    term2 = np.abs((r / Rp)**(-1.5) - 1)**(-0.5)
    return coeff * term1 * term2

def Lambda_fu(r, Rp, csp, hr, gamma, q, p):
    """Eq. (28) Bollati et al. 2021
    """
    coeff = 2**0.75 * csp * hr**(-0.5) / (gamma + 1)
    term1 = np.abs((r / Rp)**(-1.5) - 1)**0.5
    term2 = (r / Rp)**(0.5 * (p + q - 1))
    return coeff * term1 * term2

def Lambda_fv(r, Rp, csp, hr, gamma, q, p):
    """Eq. (29) Bollati et al. 2021
    """
    coeff = 2**0.75 * csp * hr**0.5 / (gamma + 1)
    term1 = np.abs((r / Rp)**(-1.5) - 1)**(-0.5)
    term2 = (r / Rp)**(0.5 * (p - q - 3))
    return coeff * term1 * term2

def get_chi(
    pphi, 
    rr, 
    time_outer,
    time_inner, 
    eta_outer, 
    eta_inner, 
    eta_tilde_outer, 
    eta_tilde_inner,
    C_outer, 
    C_inner,
    solution_outer, 
    solution_inner, 
    t0_outer, 
    t0_inner,
    tf_outer,
    tf_inner, 
    Rp, 
    x_match, 
    l, 
    cw, 
    hr, 
    q, 
    p
):
    # COMPUTATION OF Chi

    if (rr > (Rp - x_match*l)) and (rr <(Rp + x_match*l)): # exclude points inside annulus from linear regime
        return 0.

    # change coordinates of the grid point (rr,pphi) to (t1,eta1)
    t1 = t(rr, Rp, hr, q, p)
    eta1 = Eta(rr, pphi, Rp, hr, q, cw)

    # If the point is in the outer disk, use the outer wake solution
    if (rr - Rp) > 0:

        # use numerical solution before the profile develops N-wave
        if t1 < (tf_outer + t0_outer):   
            index_t = np.argmax((t0_outer + time_outer) > t1)
            grid_t = [t0_outer + time_outer[index_t-1], t0_outer+time_outer[index_t]]

            if eta1 > eta_outer[-1] or eta1 < eta_outer[0]:
                Chi = 0
            else:
                index_eta = np.argmax(eta_outer > eta1)
                grid_eta  = np.array([eta_outer[index_eta-1], eta_outer[index_eta]])
                grid_solution = np.array([
                    [solution_outer[index_eta-1,index_t-1], solution_outer[index_eta-1,index_t]],
                    [solution_outer[index_eta,index_t-1],   solution_outer[index_eta,index_t]]
                    ])

                inter = RectBivariateSpline(grid_eta, grid_t, grid_solution, kx=1, ky=1)
                Chi   = inter(eta1, t1)

        # for large t use the analytical N-wave shape (Eq. 17 Bollati et al. 2021)
        else: 

            extr_left  = +cw * np.sign(rr - Rp) * eta_tilde_outer - np.sqrt(2 * C_outer * (t1 - t0_outer))
            extr_right = +cw * np.sign(rr - Rp) * eta_tilde_outer + np.sqrt(2 * C_outer * (t1 - t0_outer))

            if eta1 > extr_left and eta1 < extr_right:
                Chi = (-cw * np.sign(rr - Rp) * eta1 + eta_tilde_outer) / (t1 - t0_outer)  #eq.(29) nonlinear.pdf
            else:
                Chi = 0

    # If the point is in the inner disk, use the inner wake solution
    else:

        # use numerical solution before the profile develops N-wave
        if t1 < (tf_inner + t0_inner):   
            index_t = np.argmax((t0_inner + time_inner) > t1)
            grid_t = [t0_inner + time_inner[index_t-1], t0_inner+time_inner[index_t]]

            if eta1 > eta_inner[-1] or eta1 < eta_inner[0]:
                Chi = 0
            else:
                index_eta = np.argmax(eta_inner > eta1)
                grid_eta  = np.array([eta_inner[index_eta-1], eta_inner[index_eta]])
                grid_solution = np.array([
                    [solution_inner[index_eta-1,index_t-1], solution_inner[index_eta-1,index_t]],
                    [solution_inner[index_eta,index_t-1],   solution_inner[index_eta,index_t]]
                    ])

                inter = RectBivariateSpline(grid_eta, grid_t, grid_solution, kx=1, ky=1)
                Chi   = -inter(eta1, t1)

        # for large t use the analytical N-wave shape (Eq. 17 Bollati et al. 2021)
        else: 

            extr_left  = +cw * np.sign(rr - Rp) * eta_tilde_inner - np.sqrt(2 * C_inner * (t1 - t0_inner))
            extr_right = +cw * np.sign(rr - Rp) * eta_tilde_inner + np.sqrt(2 * C_inner * (t1 - t0_inner))

            if eta1 > extr_left and eta1 < extr_right:
                Chi = (-cw * np.sign(rr - Rp) * eta1 + eta_tilde_inner) / (t1 - t0_inner)  #eq.(29) nonlinear.pdf
            else:
                Chi = 0
    
    return Chi

def get_dens_vel(rr, Chi, gamma ,Rp, cw, csp, hr, q, p):
    g1  = g(rr, Rp, hr, q, p)
    dnl = Chi * 2 / (g1 * (gamma + 1))     # Eq. (11) Bollati et al. 2021

    Lfu = Lambda_fu(rr, Rp, csp, hr, gamma, q, p)
    Lfv = Lambda_fv(rr, Rp, csp, hr, gamma, q, p)
    unl = np.sign(rr - Rp) * Lfu * Chi           # Eq. (23) Bollati et al. 2021
    vnl = np.sign(rr - Rp) * Lfv * Chi * (-cw) # Eq. (24) Bollati et al. 2021 (the sign of v is reversed if we change cw)

    return dnl, unl, vnl

def plot_r_t(params):
    r = np.linspace(params.r_planet, params.r_outer, 1000)
    times = []
    for rad in r:
        times.append(t(rad, params.r_planet, params.hr_planet, params.q, params.p))

    plt.plot(r, times)
    plt.xlabel("r")
    plt.ylabel('t')
    plt.show()