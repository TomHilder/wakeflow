import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import RectBivariateSpline

# TODO: Speed up functions with Numba, reorganise, optimise


def phi_wake(r, Rp, hr, q, cw): # Eq. (4) Bollati et al. 2021

   rr = r / Rp
   return -cw * np.sign(r-Rp)*(1/hr)*( rr**(q-0.5)/(q-0.5) - rr**(q + 1)/(q+1)-3/((2*q-1)*(q+1)) )

def Eta(r, phi, Rp, hr, q, cw): # Eq. (14) Bollati et al. 2021

   coeff = 1.5/hr
   phi_w = mod2pi(phi_wake(r, Rp, hr, q, cw))
   deltaphi = phi - phi_w

   if deltaphi > np.pi:
       deltaphi = deltaphi -2 * np.pi
   elif deltaphi < -np.pi:
       deltaphi = deltaphi + 2 * np.pi

   return coeff * deltaphi

def mod2pi(phi):
   if phi>=0:
      return phi%(2*np.pi)
   else:
      if np.abs(phi)<np.pi*2:
         return phi+2*np.pi
      else:
         phi = -phi
         resto = phi%(2*np.pi)
      return -resto+2*np.pi

# Equation (43) Rafikov 2002    (Eq. 13 Bollati et al. 2021)

def t_integrand(x, q, p):
   rho = 5*q + p
   w = rho/2 - 11/4
   return np.abs( 1 - x**(1.5) )**(1.5) * x**w

def t_integral(up, q, p):
   return  quad(t_integrand, 1, up, args=(q,p))[0]

def t(r, Rp, hr, q, p):
   module_integral = np.abs( t_integral(r/Rp, q, p) )
   coeff = 3*hr**(-5/2)/(2**(5/4))
   return coeff*module_integral

def xy_to_etat(x_prime, y_prime, Rp, hr, q, p, cw):

    # get r,phi
    """
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)
    """

    hp = hr * Rp

    x = 2 * hp * x_prime / 3.
    y = 2 * hp * y_prime / 3.

    r = x + Rp
    phi = y / Rp


    # I want to do this code:
    """
    t_ = t(r, Rp, hr, q, p))
    eta_ = Eta(r, phi, Rp, hr, q, p)
    """

    # but instead I gotta do this:
    t_ = np.zeros(r.shape)
    eta_ = np.zeros(r.shape)
    for i in range(0, r.shape[0]):
        for j in range(0, r.shape[1]):
            t_[i,j] = t(r[i,j], Rp, hr, q, p)
            eta_[i,j] = Eta(r[i,j], phi[i,j], Rp, hr, q, p)

    return t_, eta_

# Equation (12) Bollati et al. 2021

def g(r, Rp, hr, q, p):

    coeff = 2**0.25*hr**0.5
    term1 = (r/Rp)**(0.5*(1-p-3*q))
    term2 = np.abs( (r/Rp)**(-1.5)-1 )**(-0.5)
    return coeff * term1 * term2

def Lambda_fu(r, Rp, csp, hr, gamma, q, p):       # Eq. (28) Bollati et al. 2021

    coeff = 2**0.75*csp*hr**(-0.5)/(gamma+1)
    term1 = np.abs( (r/Rp)**(-1.5)-1 )**0.5
    term2 = (r/Rp)**(0.5*(p+q-1))
    return coeff * term1 * term2

def Lambda_fv(r, Rp, csp, hr, gamma, q, p):      # Eq. (29) Bollati et al. 2021

    coeff = 2**0.75*csp*hr**0.5/(gamma+1)
    term1 = np.abs( (r/Rp)**(-1.5)-1 )**(-0.5)
    term2 = (r/Rp)**(0.5*(p-q-3))
    return coeff * term1 * term2

def get_chi(pphi, rr, time, eta, eta_inner, eta_tilde, C, solution, solution_inner, t0, tf, Rp, x_match, l, cw, hr, q, p):
    # COMPUTATION OF Chi


    if (rr > (Rp - x_match*l)) and (rr <(Rp + x_match*l)): # exclude points inside annulus from linear regime
        return 0.


    # change coordinates of the grid point (rr,pphi) to (t1,eta1)
    t1 = t(rr, Rp, hr, q, p)
    eta1 = Eta(rr, pphi, Rp, hr, q, cw)

    if t1 < (tf + t0):   # use numerical solution before the profile develops N-wave
        index_t = np.argmax( (t0 + time) > t1 )
        #print(index_t)
        grid_t = [t0+time[index_t-1], t0+time[index_t]]

        # the density (Chi) azimuthal profile is flipped along the azimuthal direction
        # both passing from r > Rp to r < Rp and from cw = -1 to cw = +1:

        if cw*(rr-Rp) < 0:
            if eta1 > eta[-1] or eta1 < eta[0]:
                Chi = 0
            else:
                index_eta = np.argmax( eta > eta1 )
                grid_eta = np.array([eta[index_eta-1],eta[index_eta]])
                grid_solution = np.array([
                    [solution[index_eta-1,index_t-1], solution[index_eta-1,index_t]],
                    [solution[index_eta,index_t-1], solution[index_eta,index_t]]
                    ])

                inter = RectBivariateSpline(grid_eta, grid_t, grid_solution, kx=1, ky=1)
                Chi = inter(eta1, t1)

        else:
            if eta1 > eta_inner[-1] or eta1 < eta_inner[0]:
                Chi = 0
            else:
                index_eta = np.argmax( eta_inner > eta1)
                grid_eta = np.array([eta[index_eta-1],eta[index_eta]])
                grid_eta = np.array([eta[index_eta-1],eta[index_eta]])
                grid_solution = np.array([
                    [solution_inner[index_eta-1,index_t-1], solution_inner[index_eta-1,index_t]],
                    [solution_inner[index_eta,index_t-1], solution_inner[index_eta,index_t]]
                    ])

                inter = RectBivariateSpline(grid_eta, grid_t, grid_solution, kx=1, ky=1)
                Chi = inter(eta1, t1)

    else: # for large t use the analytical N-wave shape (Eq. 17 Bollati et al. 2021)

        extr_left = +cw * np.sign(rr-Rp) * eta_tilde - np.sqrt(2*C*(t1-t0))
        extr_right = +cw * np.sign(rr-Rp) * eta_tilde + np.sqrt(2*C*(t1-t0))

        if eta1 > extr_left and eta1 < extr_right:
            Chi = (-cw * np.sign(rr-Rp) * eta1 + eta_tilde) / (t1-t0)  #eq.(29) nonlinear.pdf
        else:
            Chi = 0

    return Chi

def get_dens_vel(rr, Chi, gamma ,Rp, cw, csp, hr, q, p):
    g1 = g(rr, Rp, hr, q, p)
    dnl = Chi * 2 / (g1*(gamma + 1))     # Eq. (11) Bollati et al. 2021

    Lfu = Lambda_fu(rr, Rp, csp, hr, gamma, q, p)
    Lfv = Lambda_fv(rr, Rp, csp, hr, gamma, q, p)
    unl = np.sign(rr-Rp) * Lfu * Chi           # Eq. (23) Bollati et al. 2021
    vnl = np.sign(rr-Rp) * Lfv * Chi * (-cw) # Eq. (24) Bollati et al. 2021 (the sign of v is reversed if we change cw)

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