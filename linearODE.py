import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.integrate as inte

#from odeintw import odeintw
from matplotlib.colors import LogNorm

# DE system to be solved:
"""
dv1/dt = v2
dv2/dt = - [k_y^2 (t^2 + 1) + 4/9]*v1 + M*[ -v*(v^2 + 4)/(v^2 + 1)^(3/2) ]
"""

# IC Values:
"""
v  = 0 at t = -tau_max
v' = 0 at t = -tau_max
"""

# Solution parameters

M_p = 1
m_th = 1

# Fourier space grid parameters

N_x = 2**8     # must be even
N_y = 2**8     # must be even

tau_max = np.sqrt(N_y * np.pi / 8)
k_y_nyquist = 8


k_y = np.fft.fftfreq(N_y) * (2*k_y_nyquist)

tau = np.fft.fftfreq(N_x) * (2*tau_max)
t = np.sort(tau)

# DE functions

def solve_DE_system(k_y):

    def vfunc(v, t, k_y, M):
        v1, v2 = v
        return [ v2, - (k_y**2*(t**2 + 1) + 4/9)*v1 + M*( t*(t**2 + 4)/((t**2 + 1)**(3/2)) ) ]

    def vjac(v, t, k, M):
        jac = np.array([
            [0, 1],
            [- (k**2*(t**2 + 1) + 4/9), 0]
        ])
        return jac

    v0 = [0+0j, 0+0j]
    M = -1 * np.sign(k_y) * (M_p / m_th) * (2 * np.pi *1j / 3)

    v = None #odeintw(vfunc, v0, t, args=(k_y, M), Dfun=vjac)
    
    return v[:, 0]

# pitch angle filters

def y_filter(k_y):
    if np.abs(k_y) < k_y_nyquist/2:
        return 1
    elif np.abs(k_y) < k_y_nyquist:
        return 2*(1 - np.abs(k_y)/k_y_nyquist)
    else:
        return 0

def t_vector_filter(v, t_array):
    
    def t_filter(t):
        if t < tau_max/2:
            return 1
        elif t < tau_max:
            return 2*(1 - t/tau_max)
        else:
            return 0
    
    for i, tau in enumerate(t_array):
        v[i] *= t_filter(tau)
    
    return v

def sort_for_fft(v, N):
    return np.append(v[int(N/2):], v[:int(N/2)])

"""
v = np.ones(N_x)
print(t_vector_filter(v, t))
"""


v = []

for i, k in enumerate(k_y):

    # solve DE system
    vsol = solve_DE_system(k)

    # apply pitch angle filters
    vsol *= y_filter(k)
    vsol = t_vector_filter(vsol, t)

    # reorder for ifft
    sort_for_fft(vsol, N_x)
    
    # add solution and k_x
    v.append(vsol)

v = np.array(v)

print(v)

np.fft.ifft2(v)

plt.imshow(np.abs(v), norm=LogNorm())
plt.show()





"""
# DE plots

color1 = (0.5, 0.4, 0.3)
color2 = (0.2, 0.2, 1.0)
plt.plot(t, v[:, 0].real, color=color1, label='y1.real', linewidth=1.5)
plt.plot(t, v[:, 0].imag, '--', color=color1, label='y1.imag', linewidth=2)
plt.plot(t, v[:, 1].real, color=color2, label='y2.real', linewidth=1.5)
plt.plot(t, v[:, 1].imag, '--', color=color2, label='y2.imag', linewidth=2)
plt.xlabel(r'$\tau$')
plt.grid(True)
plt.legend(loc='best')
plt.show()

y = np.fft.fft(v)
plt.plot(y)
plt.show()
"""