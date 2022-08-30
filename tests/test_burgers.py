from wakeflow.burgers import _solve_burgers
import numpy              as np
import matplotlib.pyplot  as plt

def test_burgers_area_conservation():

    eta       = np.linspace(-20, 20, 100)
    profile   = -np.sin(np.pi*eta / 40)
    gamma     = 1
    beta_p    = 1
    C         = 50
    CFL       = 0.5
    eta_tilde = 0
    t0        = 0
    tf_fac    = 1

    time, eta, sol = _solve_burgers(
        eta,
        profile,
        gamma,
        beta_p,
        C,
        CFL,
        eta_tilde,
        t0,
        0,
        0,
        False,
        tf_fac
    )

    # take every 100th solution
    sol_restr  = sol [:,::100]
    time_restr = time[::100]

    #plt.imshow(sol_restr)
    #plt.show()

    #plt.plot(eta, sol_restr[:,0])
    #plt.plot(eta, sol_restr[:,10])
    #plt.plot(eta, sol_restr[:,-1])
    #plt.show()

    # calculate area under curve with trapezoidal rule
    areas = np.trapz(sol_restr, eta, axis=0)

    # check that all areas are very close to zero
    for area in areas:
        assert area < 1e-14