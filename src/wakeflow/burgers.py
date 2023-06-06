# burgers.py
# Written by Thomas Hilder and Francesco Bollati

"""
Contains the solve_burgers function used in the wake non-linear propagation calculations.
"""

import numpy                as np
import matplotlib.pyplot    as plt
import matplotlib.colors    as cl
import matplotlib.cm        as cm

# NOTE: contents are intended for internal use and should not be directly accessed by users

# function to solve burgers equation as part of the non-linear propagation of the planet wake
def _solve_burgers(
    eta, 
    profile, 
    gamma, 
    beta_p, 
    C, 
    CFL, 
    eta_tilde, 
    t0, 
    linear_solution, 
    linear_t, 
    show_teta, 
    tf_fac,
    t_edge
): 
    """Propagate the wake in (t,eta,chi) space by solving Eq. 10 from Bollati et al. 2021 using Godunov scheme.
    """

    #Get eta step
    deta = eta[1] - eta[0]
    #print(len(eta))

    # time required to develop N-wave for betap = 1
    tf_th = 300 

    # time required to display N-wave for generic betap, Eq. (39) Rafikov 2002
    tf  = tf_th / beta_p 
    tf *= tf_fac
    
    """
    # Eq. (18) Bollati et al. 2021
    eta_min = -eta_tilde - np.sqrt(2 * C * tf) - 3 
    eta_max = -eta_tilde + np.sqrt(2 * C * tf) + 3

    # extend the profile domain due to profile broadening during evolution 
    extr = eta[-1]

    while extr < eta_max:

        eta     = np.append(eta, eta[-1] + deta)
        profile = np.append(profile, 0)
        extr    = eta[-1]

        #if show_teta:
        #    linear_solution = np.append(profile, 0, axis=0)

    extr = eta[0]

    while extr > eta_min:

        eta     = np.insert(eta, 0, eta[0] - deta)
        profile = np.insert(profile, 0, 0)
        extr    = eta[0]

        #if show_teta:
        #    linear_solution = np.insert(linear_solution, 0, 0, axis=0)
    """
    # number of centers
    Neta = len(eta)

    a = eta[0] - deta/2

    # cells edges
    x = np.zeros(Neta + 1) 
    for i in range(0, Neta + 1):
        x[i] = a + i * deta

    solution = [np.zeros((Neta), dtype=float)]
    time = [0]

    # linear solution as initial condition
    solution[0] = profile

    # define flux vector
    F = np.zeros(Neta + 1) #there is a flux at each cell edge

    # calculate flux according to Godunov scheme, vectorised
    def _NumericalFluxVector(uL, uR):
        """Calculate flux according to the Godunov scheme.
        """

        FL = 0.5 * uL**2
        FR = 0.5 * uR**2

        # compute the shock speed
        s = 0.5 * (uL + uR)

        # See "GodunovNumericalFlux".
        # This returns the same value without branching
        # through the if statements.
        return  (uL >= uR) * (
            (s > 0.0) * FL + (s <= 0.0) * FR
            ) + (uL < uR) * (
                (uL > 0.0) * (FL) +
                (uL <= 0.0) * (uR < 0.0) * FR
            )

    # initialise
    lapsed_time = 0
    counter     = 0
    """
    # time integration
    while lapsed_time < tf:

        # increment
        counter += 1

        # time step calculation
        dt = min(deta * CFL / (max(abs(solution[-1])) + 1e-8), 0.02)

        # update time
        time.append(time[-1] + dt)
        lapsed_time += dt

        # compute the interior fluxes
        F[1:Neta] = _NumericalFluxVector(solution[-1][0:-1], solution[-1][1:])

        # compute the left boundary flux
        uL   = solution[-1][-1]
        uR   = solution[-1][0]
        F[0] = _NumericalFluxVector(uL, uR)

        # compute the right boundary flux
        uR      = solution[-1][0]
        uL      = solution[-1][Neta - 1]
        F[Neta] = _NumericalFluxVector(uL, uR)

        solution.append(solution[-1][0:Neta] - dt / deta * (F[1:Neta+1] - F[0:Neta]))

        # Eq. (18) Bollati et al. 2021
        eta_minus = eta_tilde - np.sqrt(2 * C * tf)
        eta_plus = eta_tilde + np.sqrt(2 * C * tf)

        if lapsed_time >= t_edge:
            break
    """
    # time integration until boundary of the disc
    #print(t_edge)
    while lapsed_time <= t_edge:

        # increment
        counter += 1

        # time step calculation
        dt = min(deta * CFL / (max(abs(solution[-1])) + 1e-8), 0.02)

        # update time
        time.append(time[-1] + dt)
        lapsed_time += dt

        # compute the interior fluxes
        F[1:Neta] = _NumericalFluxVector(solution[-1][0:-1], solution[-1][1:])

        # compute the left boundary flux
        uL   = solution[-1][-1]
        uR   = solution[-1][0]
        F[0] = _NumericalFluxVector(uL, uR)

        # compute the right boundary flux
        uR      = solution[-1][0]
        uL      = solution[-1][Neta - 1]
        F[Neta] = _NumericalFluxVector(uL, uR)

        solution.append(solution[-1][0:Neta] - dt / deta * (F[1:Neta+1] - F[0:Neta]))

        # Implementing condition for exiting the loop in the asymptotic limit, when N wave and numerical solutions coincide up to 5%
        # Eq. (18) Bollati et al. 2021.
        # Left and right extrema of N wave profile
        eta_minus = eta_tilde - np.sqrt(2 * C * (time[-1] - t0))
        eta_plus  = eta_tilde + np.sqrt(2 * C * (time[-1] - t0))

        # Maximum and minimum values of N wave Eq. (17) Bollati et al. 2021.
        N_wave_max = (eta_plus  - eta_tilde) / (time[-1]-t0)
        N_wave_min = (eta_minus - eta_tilde) / (time[-1]-t0)
        # Conditions on minimum and maximum amplitude up to 5%
        condition_max_amplitude = np.abs((N_wave_max - np.max(solution[-1])) / np.max(solution[-1])) <= 0.05
        condition_min_amplitude = np.abs((N_wave_min - np.min(solution[-1])) / np.min(solution[-1])) <= 0.05 #to be improved, now works only for outer disc. Or not?
        condition_amplitude     = np.logical_and(condition_max_amplitude,condition_min_amplitude)

        # Conditions on left and right extrema up to 5%
        condition_max_width = np.abs(eta_plus  - eta[np.argmax(solution[-1])]) <= 0.05
        condition_min_width = np.abs(eta_minus - eta[np.argmin(solution[-1])]) <= 0.05
        condition_width     = np.logical_and(condition_max_width, condition_min_width)


        if lapsed_time >= tf:
            if np.logical_and(condition_amplitude, condition_width):
                break

    solution = np.array(solution).transpose()
    time     = np.array(time)
   
    if False:
        #plots for debugging: plot \eta profiles to check the evolution of the solution
        #colorbar
        n = len(time)
        norm = cl.Normalize(vmin=time.min(), vmax=time.max())
        cmap = cm.ScalarMappable(norm=norm, cmap=cm.viridis)
        cmap.set_array([])
        
        plt.figure(figsize=(15,5))
        
        idxp = int(n/5) #index to plot 5 profiles equidistant in \t
        plt.plot(eta, solution[:,0*idxp], label = "$t=t_0+%.1lf$ "%time[0*idxp],
                color = cmap.to_rgba(time[0*idxp] + 1), ls='-.')#+str(0*dt))
        plt.plot(eta, solution[:,1*idxp], label = "$t=t_0+%.1lf$ "%time[1*idxp],
                color = cmap.to_rgba(time[1*idxp] + 1), ls='-.')#+str(round(dt*Nt/10,2)))
        plt.plot(eta, solution[:,2*idxp], label = "$t=t_0+%.1lf$ "%time[2*idxp],
                color = cmap.to_rgba(time[2*idxp] + 1), ls='-.')#+str(round(dt*Nt/6,2)))
        plt.plot(eta, solution[:,3*idxp], label = "$t=t_0+%.1lf$ "%time[3*idxp],
                color = cmap.to_rgba(time[3*idxp] + 1), ls='-.')#+str(round(dt*Nt/3,2)))
        plt.plot(eta, solution[:,4*idxp], label = "$t=t_0+%.1lf$ "%time[4*idxp],
                color = cmap.to_rgba(time[4*idxp] + 1), ls='-.')#+str(round(T,2)))
        plt.plot(eta, solution[:,-1], label = "$t=t_0+%.1lf$ "%time[-1],
                color = cmap.to_rgba(time[-1]), ls='-.')#+str(round(T,2)))

        #check last profile with analytic N wave to see if asymptotic limit already reached
        if eta_tilde <= 0:
            eta_plus = eta_tilde + np.sqrt(2*C*(time[-1]-t0))
            eta_minus = eta_tilde - np.sqrt(2*C*(time[-1]-t0))
            #plt.plot(eta, np.where(np.logical_and(eta>eta_minus,eta<eta_plus),(eta-eta_tilde)/(time[-1]-t0),0), label='N')
        else:
            eta_plus = eta_tilde + np.sqrt(2*C*(time[-1]-t0))
            eta_minus = eta_tilde - np.sqrt(2*C*(time[-1]-t0))
            #plt.plot(eta, np.where(np.logical_and(eta>eta_minus,eta<eta_plus),-(eta-eta_tilde)/(time[-1]-t0),0), label='N')
        
        plt.legend()
        plt.xlim(-40,40)
        plt.xlabel(r"$\eta$")
        plt.ylabel(r"$\chi(t,\eta)$")
        plt.grid(True)
        cb = plt.colorbar(cmap, label = r't', shrink = 1)
        #cb.ax.tick_params(labelsize=20)
        cb.ax.set_ylabel(r't')#, fontsize=20)
        cb.ax.yaxis.set_label_position('left')
        cb.ax.set_aspect('auto')
        #Trying to get nice double colorbars, failing miserably
        """
        cb2 = cb.ax.twinx() 
        cb2.set_ylim([time[0], time[-1]])
        
        rticks=np.zeros(np.shape(time))
        for i in range(len(time)):
            if i%idxp == 0 or i == len(time)-1:
                def invert_t(x, Rp, hr, q, p, m_p, m_th):
                    return _t_vector(x, Rp, hr, q, p, m_p, m_th) - time[i]
                rticks[i] = fsolve(invert_t, args=(Rp, hr, q, p, m_p, m_th), x0 = 120)
        cb2.set_yticks([time[0*idxp],time[1*idxp],time[2*idxp],time[3*idxp],time[4*idxp],time[-1]])
        cb2.set_yticklabels([round(rticks[0*idxp]),round(rticks[1*idxp]),round(rticks[2*idxp]),
                            round(rticks[3*idxp]),round(rticks[4*idxp]),round(rticks[-1])])
        #cb2.set_yticks([tic for tic in pdfp])
        cb2.set_ylabel(r'R/R$_{\rm p}$', labelpad=8)#, fontsize=20)
        #cb2.tick_params(labelsize=20)
        """
        if eta_tilde <= 0:
            plt.title(r'$\chi$ "evolution" $r > r_p$')
        else:
            plt.title(r'$\chi$ "evolution" $r < r_p$')
        plt.show()
        
        #Like before, check last profile with analytic N wave to see if asymptotic limit already reached
        plt.plot(eta, solution[:,-1], label = "$t=t_0+%.1lf$ "%time[-1])#+str(round(T,2)))
        #plt.plot(eta, solution[:,-1], label = "$t=t_0+$ ")#+str(round(T,2)))
        if eta_tilde <= 0:
            plt.plot(eta, np.where(np.logical_and(eta>eta_minus,eta<eta_plus),(eta-eta_tilde)/(time[-1]-t0),0), label='N')
        else:
            plt.plot(eta, np.where(np.logical_and(eta>eta_minus,eta<eta_plus),(eta-eta_tilde)/(time[-1]-t0),0), label='N')
        plt.legend()
        plt.xlabel(r"$\eta$")
        plt.ylabel(r"$\chi(t,\eta)$")
        plt.grid(True)
        #plt.xticks([-28, -27, -26, -25, -24, 0, 20, 21, 22])
        plt.show()
        
    #plt.plot(time, rticks)

    if show_teta: # combining linear and non-linear solution and plotting
        
        # scale linear solution
        linear_solution = linear_solution * (gamma + 1) * beta_p / 2**(3 / 4)
    
        # add linear solution in (t,eta) to non-linear solution array
        total_solution = np.concatenate((linear_solution, solution), axis=1)
        total_time     = np.concatenate((linear_t,        time + t0))
    
        fig, ax = plt.subplots(1)
        cont = ax.contourf(total_time, eta, total_solution, levels=np.arange(-4, 4, 0.05), cmap='RdBu')
        for c in cont.collections:
            c.set_rasterized(True)
        plt.colorbar(cont, label=r'$\chi$')
        ax.set_xlim(0,10)
        ax.set_xlabel(r'$t$')
        ax.set_ylabel(r'$\eta$')
    #    #plt.savefig("teta_badjoin.pdf")
        plt.show()

    # The following will put the linear solution into returned solution for Chi
    #if False:
    #    Nt = np.shape(total_solution)[1]
    #
    #    solution_inner = np.zeros(total_solution.shape) # solution for r < Rp (and disc rotating counterclockwise)
    #    for i in range(Neta):
    #        for j in range(Nt):
    #            solution_inner[i,j] = total_solution[int(Neta-1-i),j]
    #
    #    eta_inner = - eta[::-1]
    #    return total_time, eta, total_solution, eta_inner, solution_inner

    return time, eta, solution
    
