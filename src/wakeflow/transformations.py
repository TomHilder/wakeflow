# transformations.py
# Written by Thomas Hilder and Francesco Bollati

"""
Contains functions for mapping the results between (t,eta,chi) coordinates and physical coordinates, see Bollati et al. 2021 for details.
"""

import numpy                    as np
import matplotlib.pyplot        as plt
from scipy.integrate        import quad, odeint
from scipy.interpolate      import RectBivariateSpline

# NOTE: contents are intended for internal use and should not be directly accessed by users

# TODO: update code to be faster and more efficient, most of the code in this file is still leftover from the old version (Analytical Kinks)

# wake shape
def _phi_wake(r, Rp, hr, q, cw): 
    """Eq. (4) Bollati et al. 2021
    """
    rr = r / Rp
    return -cw * np.sign(r - Rp) * (1 / hr) * (rr**(q - 0.5) / (q - 0.5) - rr**(q + 1) / (q + 1) - 3 / ((2 * q - 1) * (q + 1)))

# eta coordinate transformation
def _Eta(r, phi, Rp, hr, q, cw):
    """Eq. (14) Bollati et al. 2021
    """
    coeff    = 1.5 / hr
    phi_w    = _mod2pi(_phi_wake(r, Rp, hr, q, cw))
    deltaphi = phi - phi_w

    if deltaphi > np.pi:
        deltaphi = deltaphi - 2 * np.pi

    elif deltaphi < -np.pi:
        deltaphi = deltaphi + 2 * np.pi 

    return coeff * deltaphi

def _Eta_vector(r, phi, Rp, hr, q, cw):
    """Eq. (14) Bollati et al. 2021

    Vectorised version of _Eta.
    Simply replaces the modular arithmetic which involved if statements with
    modulus operators and constant offsets.
    """
    coeff    = 1.5 / hr
    phi_w    = _phi_wake(r, Rp, hr, q, cw) % (2*np.pi)
    deltaphi = (phi - phi_w + np.pi) % (2*np.pi) - np.pi

    return coeff * deltaphi

# modulo 2pi
def _mod2pi(phi):
    if phi >= 0:
        return phi % (2 * np.pi) 

    else:
        if np.abs(phi) < np.pi * 2:
            return phi + 2 * np.pi    

        else:
            phi   = -phi
            resto =  phi % (2 * np.pi)    
            
        return -resto + 2 * np.pi


# integrand of the t coordinate transformation
def _t_integrand(x, q, p):
    rho = 5 * q + p
    w   = rho / 2 - 11 / 4
    return np.abs(1 - x**(1.5))**(1.5) * x**w

# integral for t coordinate transformation
def _t_integral(up, q, p):
    return  quad(_t_integrand, 1, up, args=(q,p))[0]

def _t_vector(rr, Rp, hr, q, p):
    """Equation (43) Rafikov 2002    (Eq. 13 Bollati et al. 2021)
    This is a vectorised version of _t.
    _t computes an integral where the integrand is independent of the radius r.
    The radius r only changes the end point of the integral. Instead of using quad
    for each end point (which is not as easily vectorisable), we can utilise an ODE solver to 
    obtain the result at every end point without doing redundant work.
    """

    # First we will flatten the array, keeping track of the original shape.
    # We will also obtain the index of each element of the sorted array.
    # this will allow us to "unsort" the array at the end.
    shape = rr.shape
    integral_bounds = (rr / Rp).flatten()
    integral_bounds_sorted_indices = np.argsort(integral_bounds)

    # We must divide the integral endpoints into values above one and below one,
    # since the lower integral bound is one, and the ode solver only accepts time
    # values in ascending order.
    # This is effectively dividing the integrals between ones which integrate backwards,
    # and ones which integrate forwards.
    integral_bounds_mask  = integral_bounds >= 1.0
    num_larger = np.sum(integral_bounds_mask)

    # We will allocate the arrays to store the values >= 1 and < 1, prepending
    # the initial condition, 1.0. Then we will sort the array since this is required
    # by the ode solver.
    larger_vals = np.zeros( (num_larger + 1,) , dtype = integral_bounds.dtype)
    larger_vals[0] = 1.0
    larger_vals[1:] = integral_bounds[integral_bounds_mask]
    larger_vals.sort()

    smaller_vals = np.zeros( (integral_bounds.size - num_larger + 1,),
        dtype = integral_bounds.dtype)
    smaller_vals[0] = 1.0
    smaller_vals[1:] = integral_bounds[~integral_bounds_mask]
    smaller_vals.sort()
    smaller_vals = np.flip(smaller_vals)

    # Now we can evaluate the integral using the ode solver.
    rho = 5 * q + p
    w   = rho / 2 - 11 / 4

    odefun = lambda _,x: np.abs(1 - x**(1.5))**(1.5) * x**w
    larger_t = odeint(odefun, np.array([0]), larger_vals)

    # In the case of the lower ones, we are making a linear coordinate 
    # transformation x -> 1 - x to make the values ascending and starting at
    # 0 rather than -1 to avoid internal problems in the solver.
    odefun = lambda _,x: np.abs(1 - (1-x)**(1.5))**(1.5) * (1-x)**w
    smaller_t = odeint(odefun, np.array([0]), 1-smaller_vals)

    # combine the results by excluding the prepended initial condition.
    sorted_results = np.concatenate((
        np.flip(smaller_t.flatten()[1:]), 
        larger_t.flatten()[1:]
    ))

    # "unsort", reshape and apply final calculations on the results.
    module_integral = np.zeros_like(sorted_results)
    module_integral[integral_bounds_sorted_indices] = np.abs(sorted_results)

    coeff = 3 * hr**(-5 / 2) / (2**(5 / 4))
    return np.reshape(coeff * module_integral, shape)


# t coordinate transformation
def _t(r, Rp, hr, q, p):
    """Equation (43) Rafikov 2002    (Eq. 13 Bollati et al. 2021)
    """
    module_integral = np.abs(_t_integral(r / Rp, q, p))
    coeff = 3 * hr**(-5 / 2) / (2**(5 / 4))
    return coeff * module_integral

# g(r) quantity calculation
def _g(r, Rp, hr, q, p):
    """Equation (12) Bollati et al. 2021
    """
    coeff = 2**0.25 * hr**0.5
    term1 = (r / Rp)**(0.5 * (1 - p - 3 * q))
    term2 = np.abs((r / Rp)**(-1.5) - 1)**(-0.5)
    return coeff * term1 * term2

# needed to get radial velocities
def _Lambda_fu(r, Rp, csp, hr, gamma, q, p):
    """Eq. (28) Bollati et al. 2021
    """
    coeff = 2**0.75 * csp * hr**(-0.5) / (gamma + 1)
    term1 = np.abs((r / Rp)**(-1.5) - 1)**0.5
    term2 = (r / Rp)**(0.5 * (p + q - 1))
    return coeff * term1 * term2

## needed to get azimuthal velocities
def _Lambda_fv(r, Rp, csp, hr, gamma, q, p):
    """Eq. (29) Bollati et al. 2021
    """
    coeff = 2**0.75 * csp * hr**0.5 / (gamma + 1)
    term1 = np.abs((r / Rp)**(-1.5) - 1)**(-0.5)
    term2 = (r / Rp)**(0.5 * (p - q - 3))
    return coeff * term1 * term2

# find chi for a particular grid point, either from Burger's eqn solution or self-similar solution
def _get_chi(
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
    p,
    t1
):
    # COMPUTATION OF Chi

    if (rr > (Rp - x_match*l)) and (rr <(Rp + x_match*l)): # exclude points inside annulus from linear regime
        return 0.

    # change coordinates of the grid point (rr,pphi) to (t1,eta1)
    #t1_orig = _t(rr, Rp, hr, q, p)

    #if np.abs(t1-t1_orig) / t1_orig > 1e-2:
    #    print(t1-t1_orig)


    eta1 = _Eta(rr, pphi, Rp, hr, q, cw)

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


def _get_chi_vector(
    pphi, 
    rr,
    tt,
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
    x_match_l,
    x_match_r,
    l, 
    cw, 
    hr, 
    q, 
    p,
): 
    """
    This is a vectorised version of the previous _get_chi function.
    Key changes are:
     - Using masks rather than logical statements to replace
        if statements.
     - Using one RectBivariateSpline interpolation function for all
        points needing interpolation, rather than recreating a function
        for every point.
    """
 # COMPUTATION OF Chi

    Chi = np.zeros_like(rr)

    eta_array = _Eta_vector(rr, pphi, Rp, hr, q, cw)

    # Inner and outer masks will account for the annulus directly.
    outer_mask = rr - Rp >= x_match_r * l
    inner_mask = rr - Rp <= -x_match_l * l

    before_N_wave_mask = np.logical_or(
        tt < (tf_outer + t0_outer),
        tt < (tf_inner + t0_inner)
        )

    # Outer masks.

    """
    if (rr - Rp) >= x_match*l:
        if t1 >= (tf_outer + t0_outer):
            if eta1 > extr_left and eta1 < extr_right:
    """
    m = np.logical_and(outer_mask,
        np.logical_and(~before_N_wave_mask,
        np.abs( eta_array - cw * np.sign(rr - Rp) * eta_tilde_outer) > np.sqrt(2 * C_outer * np.abs(tt - t0_outer))
        ))
    Chi[m] = 0.0

    """
    if (rr - Rp) >= x_match*l:
        if t1 >= (tf_outer + t0_outer):
            if eta1 <= extr_left or eta1 >= extr_right:
    """
    m = np.logical_and(outer_mask,
        np.logical_and(~before_N_wave_mask,
        np.abs( eta_array - cw * np.sign(rr - Rp) * eta_tilde_outer) <= np.sqrt(2 * C_outer * np.abs(tt - t0_outer))
        ))
    Chi[m] = (-cw * np.sign(rr[m] - Rp) * eta_array[m] + eta_tilde_outer) / (tt[m] - t0_outer)


    """
    if (rr - Rp) >= x_match*l:
        if t1 < (tf_outer + t0_outer):
            if eta1 > eta_outer[-1] or eta1 < eta_outer[0]:
    """
    interp_outer = RectBivariateSpline(eta_outer, t0_outer + time_outer, solution_outer, kx=1, ky=1)
    m = np.logical_and(outer_mask, 
        np.logical_and(before_N_wave_mask,
        np.logical_and(eta_outer[0] < eta_array,
            eta_array < eta_outer[-1]
        )))
    Chi[m] = interp_outer(eta_array[m], tt[m], grid=False)

    # Inner masks.
    
    """
    if (rr - Rp) <= -x_match*l:
        if t1 >= (tf_inner + t0_inner):
            if eta1 > extr_left and eta1 < extr_right:
    """
    m = np.logical_and(inner_mask,
        np.logical_and(~before_N_wave_mask,
        np.abs( eta_array - cw * np.sign(rr - Rp) * eta_tilde_inner) >= np.sqrt(2 * C_inner * np.abs(tt - t0_inner))
        ))
    Chi[m] = 0.0

    """
    if (rr - Rp) <= -x_match*l:
        if t1 >= (tf_inner + t0_inner):
            if eta1 <= extr_left or eta1 >= extr_right:
    """
    m = np.logical_and(inner_mask,
        np.logical_and(~before_N_wave_mask,
        np.abs( eta_array - cw * np.sign(rr - Rp) * eta_tilde_inner) < np.sqrt(2 * C_inner * np.abs(tt - t0_inner))
        ))
    Chi[m] = (-cw * np.sign(rr[m] - Rp) * eta_array[m] + eta_tilde_inner) / (tt[m] - t0_inner)

    """
    if (rr - Rp) <= -x_match*l:
        if t1 < (tf_inner + t0_inner):
            if eta1 <= eta_inner[-1] and eta1 >= eta_inner[0]:
    """
    interp_inner = RectBivariateSpline(eta_inner, t0_inner + time_inner, solution_inner, kx=1, ky=1)

    m = np.logical_and(inner_mask, 
        np.logical_and(before_N_wave_mask,
        np.logical_and(eta_inner[0] < eta_array,
            eta_array < eta_inner[-1]
        )))
    Chi[m] = -interp_inner(eta_array[m], tt[m], grid=False)


    return Chi

# get the density and velocity perturbations at the grid point from chi
def _get_dens_vel(rr, Chi, gamma, Rp, cw, csp, hr, q, p, use_old_vel):

    g1  = _g(rr, Rp, hr, q, p)
    dnl = Chi * 2 / (g1 * (gamma + 1))     # Eq. (11) Bollati et al. 2021
    
    if use_old_vel == True:
        Lfu = _Lambda_fu(rr, Rp, csp, hr, gamma, q, p)
        Lfv = _Lambda_fv(rr, Rp, csp, hr, gamma, q, p)
        unl = np.sign(rr - Rp) * Lfu * Chi           # Eq. (23) Bollati et al. 2021
        vnl = np.sign(rr - Rp) * Lfv * Chi * (-cw) # Eq. (24) Bollati et al. 2021 (the sign of v is reversed if we change cw)
    else:
        #psi = (np.power(dnl + 1, (gamma-1)/2) - 1) * (gamma+1) / (gamma-1)
        psi = ((gamma+1) / (gamma-1)) * np.sign(dnl + 1) * ((np.abs(dnl + 1)) ** ((gamma - 1) / 2) - 1)
    
        # get constants
        dOmega_r = np.abs(csp * Rp**-1 * hr**-1 * ((rr / Rp)**(-3 / 2) - 1)) * rr
        c0 = csp * (rr / Rp)**(-q)
        
        # get perturbations
        unl = np.sign(rr - Rp) * (2 * c0) / (gamma+1) * psi
        vnl = (-cw) * c0 * unl / dOmega_r

    return dnl, unl, vnl

# plot the t coordinate as a function of radius
def _plot_r_t(params):
    r = np.linspace(params.r_planet, params.r_outer, 1000)
    times = []
    for rad in r:
        times.append(_t(rad, params.r_planet, params.hr_planet, params.q, params.p))

    plt.plot(r, times)
    plt.xlabel("r")
    plt.ylabel('t')
    plt.show()
