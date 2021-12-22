
# === IMPORTS === #
from matplotlib.collections import PathCollection
import numpy                as np
import matplotlib.pyplot    as plt
import matplotlib.colors    as cs
import scipy.interpolate    as interp
import astropy.constants    as c
import transformations      as t
import copy                 as cp

# === PARAMS === #

# chi extraction points
CHI_RADII = [200, 220, 240]
CHI_SLICE_WIDTH_ANG = np.pi / 2

# phantom pixelmaps
PIXELMAP_LOC = "phantom/imlup_10mil"

# disc (R_INNER and R_OUTER must match pixelmap)
M_STAR      = 1.12
R_OUTER     = 300
R_INNER     = 30
R_REF       = 100
R_CRIT      = 150
P_INDEX     = 0.48
Q_INDEX     = 0.31
H_R         = 0.129
ROTATION_CW = False

# planet
M_PLANET    = 1.0
R_PLANET    = 139.207

# diagnostics
PLOTS = True
COMPONENTS_PLOT = True
CHI_2D_PLOT = True
CHI_SLICE_PLOT = True

# ==================================================================================================== #

# === constants === #

# all constants given in CGS units
M_SOL   = c.M_sun.cgs.value
M_JUP   = c.M_jup.cgs.value
G_CONST = c.G.cgs.value
AU      = c.au.cgs.value

# adiabatic index
GAMMA   = 1.6666667

# analytic linear box size (in units of Hill radius)
BOX_SIZE_HILL = 2

# === calculate needed quantities === #

# quantities at planet location
HR_PLANET = H_R * (R_PLANET / R_REF) ** (0.5 - Q_INDEX)
VK_PLANET = np.sqrt(G_CONST * M_STAR * M_SOL / (R_PLANET * AU))
CS_PLANET = VK_PLANET * HR_PLANET

# thermal mass in units of M_jup
M_THERMAL = (2 / 3) * HR_PLANET ** 3 * (M_STAR * M_SOL / M_JUP)
print(M_THERMAL)

# planet mass in units of thermal mass
BETA_P = M_PLANET / M_THERMAL

# length scale linear regime
L_SCALE = (2 / 3) * HR_PLANET * R_PLANET

# box size in physical units (au)
BOX_SIZE_AU = L_SCALE * BOX_SIZE_HILL

# clockwise constant
if ROTATION_CW:
    CW = 1
else:
    CW = -1

# calculate edge of linear regime
R_BOX = R_PLANET + 2 * L_SCALE
T_BOX = t.t(R_BOX, R_PLANET, HR_PLANET, Q_INDEX, P_INDEX)
CHI_RADII = [R_PLANET - 2 * L_SCALE] + [R_BOX] + CHI_RADII

# plotting parameters
POLAR_EXTENT = [R_INNER, R_OUTER, -np.pi, np.pi]
ASPECT = 50

# === print results info === #

print(f"Embedded planet has mass of {round(BETA_P, 3)} M_th, and linear regime ends at t = {round(T_BOX, 3)}.")

# === construct wake shape === #

print("Constructing analytic wake shape...")

# function for finding phi of the wake given a radius, vectorised and mod (-pi,pi)
def phi_wake_func(R):

    term1 =  ((R / R_PLANET)**(Q_INDEX - 0.5)) / (Q_INDEX - 0.5)
    term2 = -((R / R_PLANET)**(Q_INDEX + 1)) / (Q_INDEX + 1)
    term3 = -3 / ((2*Q_INDEX - 1) * (Q_INDEX + 1))

    phi = -CW * np.sign(R - R_PLANET) * (HR_PLANET**-1) * (term1 + term2 + term3)
    
    try:
        for i, p in enumerate(phi):
            if p < -np.pi:
                phi[i] += 2*np.pi
            if p >  np.pi:
                phi[i] -= 2*np.pi
    except:
        if phi < -np.pi:
            phi += 2*np.pi
        if phi >  np.pi:
            phi -= 2*np.pi

    return phi

# function for eta transformation, vectorised and mod (-pi,pi)
def mod_pi_Eta(r, phi):

    coeff = 1.5 / HR_PLANET
    phi_w = phi_wake_func(r)
    deltaphi = phi - phi_w

    return coeff * deltaphi

# constructing r,phi points along wake
r_wake = np.linspace(R_INNER, R_OUTER, 200)
phi_wake = phi_wake_func(r_wake)

# constructing x,y points along wake
x_wake = r_wake * np.cos(phi_wake)
y_wake = r_wake * np.sin(phi_wake)

# === read in phantom pixelmaps === #

print("Reading in Phantom pixelmaps and reconstructing grid...")

# read in midplane arrays, in cgs units
vr   = np.loadtxt(f"{PIXELMAP_LOC}/vr.pix")
vphi = np.loadtxt(f"{PIXELMAP_LOC}/vphi.pix")
rho  = np.loadtxt(f"{PIXELMAP_LOC}/rho.pix")

# get phantom grid
PH_GRID_LENGTH = vr.shape[0]
x_ph = np.linspace(-R_OUTER, R_OUTER, PH_GRID_LENGTH)
y_ph = np.linspace(-R_OUTER, R_OUTER, PH_GRID_LENGTH)
X_ph, Y_ph = np.meshgrid(x_ph, y_ph)

print("Interpolating Phantom results onto a regularly spaced polar grid...")

# interpolate over phantom grid
vr_interp   = interp.RectBivariateSpline(y_ph, x_ph, vr)
vphi_interp = interp.RectBivariateSpline(y_ph, x_ph, vphi)
rho_interp  = interp.RectBivariateSpline(y_ph, x_ph, rho)

# regularly spaced cylindrical grid
r_reg   = np.linspace(R_INNER, R_OUTER, PH_GRID_LENGTH)
phi_reg = np.linspace(-np.pi,    np.pi, PH_GRID_LENGTH)
R_reg, PHI_reg = np.meshgrid(r_reg, phi_reg)

# convert regularly space cylindrical grid to cartesian irregular grid
X_reg = R_reg * np.cos(PHI_reg)
Y_reg = R_reg * np.sin(PHI_reg)

# evalutate interpolation on polar grid
vr_polar   = vr_interp.ev  (Y_reg, X_reg)
vphi_polar = vphi_interp.ev(Y_reg, X_reg)
rho_polar  = rho_interp.ev (Y_reg, X_reg)

# convert velocities to cgs (density already cgs)
vr_polar   *= 1e5
vphi_polar *= 1e5
vr         *= 1e5
vphi       *= 1e5

# test plots
if PLOTS:
    plt.title(r"$v_r$")
    plt.imshow(vr_polar, cmap="RdBu", vmin=-0.6*min(abs(vr_polar.min()),abs(vr_polar.max())), vmax=0.6*min(abs(vr_polar.min()),abs(vr_polar.max())), origin="lower")
    plt.colorbar()
    plt.show()
    plt.title(r"$v_\phi$")
    plt.imshow(vphi_polar, cmap="inferno", vmin=0, vmax=0.5*abs(vphi_polar.max()), origin="lower")
    plt.colorbar()
    plt.show()
    plt.title(r"$\rho$")
    plt.imshow(rho_polar, cmap="inferno", norm=cs.LogNorm(vmin=np.abs(rho_polar).min(), vmax=0.01*rho_polar.max()), origin="lower")
    plt.colorbar()
    plt.show()

# === background density and Keplerian rotation subtraction === #

print("Subtracting azimuthally averaged density and rotation...")

AZ_WIDTH = 1

# calculate the mean with a lambda
w = AZ_WIDTH / 2
f = lambda r : rho_polar[(R_reg >= r-w) & (R_reg < r+w) & (PHI_reg < 7*np.pi/4) & (PHI_reg > np.pi/4)].mean()
r  = np.linspace(R_INNER, R_OUTER, num=int(R_OUTER))
mean = np.vectorize(f)(r)
mean_interp_func = interp.interp1d(r, mean, kind='linear')

# density at reference radius
RHO_REF = f(R_REF)
RHO_REF_TAPER = RHO_REF / np.exp(-(R_REF / R_CRIT)**(2-P_INDEX))

rho_power_law = RHO_REF * (r / R_REF)**(-P_INDEX)
rho_taper     = RHO_REF_TAPER * (r / R_REF)**(-P_INDEX) * np.exp(-(r / R_CRIT)**(2-P_INDEX))

# test plot
if PLOTS:
    fig,ax=plt.subplots()
    ax.plot(r, mean, label="Density Profile")
    ax.plot(r, rho_power_law, label="Power Law")
    ax.plot(r, rho_taper, label="Exponential Taper")
    ax.set_title("azimuthally averaged density")
    ax.set_yscale('log')
    ax.legend(loc="best")
    plt.show()

# subtract azimuthal average
averaged = mean_interp_func(R_reg)
rho_pert = (rho_polar - averaged) / averaged

# calculate Keplerian velocities for plotting
v_keplerian = np.sqrt(G_CONST * M_STAR * M_SOL / (r * AU))

# calculate the mean with a lambda
f = lambda r : vphi_polar[(R_reg >= r-w) & (R_reg < r+w) & (PHI_reg < 7*np.pi/4) & (PHI_reg > np.pi/4)].mean()
r  = np.linspace(R_INNER, R_OUTER, num=int(R_OUTER))
mean = np.vectorize(f)(r)
mean_interp_func = interp.interp1d(r, mean, kind='linear')

# test plot
if PLOTS:
    fig,ax=plt.subplots()
    plt.plot(r, v_keplerian, label="Keplerian")
    ax.plot(r, mean, label="Rotation curve")
    ax.set_title("azimuthally averaged rotation")
    #ax.set_yscale("log")
    ax.legend(loc="best")
    plt.show()

# subtract azimuthal average
averaged = mean_interp_func(R_reg)
vphi_pert = (vphi_polar - averaged)

# radial velocities don't need subtraction
vr_pert = cp.copy(vr_polar)

if COMPONENTS_PLOT:
    R_EXC = 50
    COMPONENTS_TITLES = [r"$\rho$",r"$v_r$",r"$v_\phi$",]
    fig, ax = plt.subplots(1, 3, figsize=[15,8], sharex='all', sharey='all')

    cw0 = ax[0].imshow(rho_pert,  cmap="RdBu", vmin=-0.1*min(abs(rho_pert.min()),abs(rho_pert.max())), vmax=0.1*min(abs(rho_pert.min()),abs(rho_pert.max())), 
        extent=POLAR_EXTENT, aspect=ASPECT, origin="lower")

    cw1 = ax[1].imshow(vr_pert,   cmap="RdBu", vmin=-0.5*min(abs(vr.min()),abs(vr.max())), vmax=0.5*min(abs(vr.min()),abs(vr.max())), 
        extent=POLAR_EXTENT, aspect=ASPECT, origin="lower")

    cw2 = ax[2].imshow(vphi_pert, cmap="RdBu", vmin=-0.5*min(abs(vr.min()),abs(vr.max())), vmax=0.5*min(abs(vr.min()),abs(vr.max())), 
        extent=POLAR_EXTENT, aspect=ASPECT, origin="lower")

    ax[0].set_ylabel(r"$\phi$ [rad]")
    cbars = [cw0, cw1, cw2]
    CBAR_LABELS = ["[g/cm^2]", "[cm/s]", "[cm/s]"]
    for i in range(3):
        ax[i].set_title(COMPONENTS_TITLES[i])
        f#ig.colorbar(cbars[i], cax=ax[i], extend="both", label=CBAR_LABELS[i], orientation="horizontal")
        ax[i].set_xlabel(r"$r$ [au]")
    plt.savefig("components.pdf", dpi=200)
    plt.show()

if False:
    plt.title(r"$v_r$")
    plt.imshow(vr_pert, cmap="RdBu", vmin=-0.6*min(abs(vr.min()),abs(vr.max())), vmax=0.6*min(abs(vr.min()),abs(vr.max())), origin="lower")
    plt.colorbar()
    plt.show()
    plt.title(r"$\Delta v_\phi$")
    R_EXC = 50
    plt.imshow(vphi_pert, cmap="RdBu", vmin=-0.1*min(abs(vphi_pert[R_reg>R_EXC].min()),abs(vphi_pert[R_reg>R_EXC].max())), 
        vmax=0.1*min(abs(vphi_pert[R_reg>R_EXC].min()),abs(vphi_pert[R_reg>R_EXC].max())), origin="lower")
    plt.colorbar()
    plt.show()
    plt.title(r"$\Delta \rho$")
    plt.imshow(rho_pert, cmap="RdBu", vmin=-0.1*min(abs(rho_pert.min()),abs(rho_pert.max())), vmax=0.1*min(abs(rho_pert.min()),abs(rho_pert.max())), origin="lower")
    plt.colorbar()
    plt.show()

# === extract chi values from each component === #

print("Extracting chi from each perturbation component...")

# density
chi_rho  = rho_pert  * ((GAMMA + 1) / 2) * t.g(R_reg, R_PLANET, HR_PLANET, Q_INDEX, P_INDEX)

# velocities
chi_vr   = vr_pert   / (np.sign(R_reg - R_PLANET) * t.Lambda_fu(R_reg, R_PLANET, CS_PLANET, HR_PLANET, GAMMA, Q_INDEX, P_INDEX))
chi_vphi = vphi_pert / (np.sign(R_reg - R_PLANET) * t.Lambda_fv(R_reg, R_PLANET, CS_PLANET, HR_PLANET, GAMMA, Q_INDEX, P_INDEX) * -CW)

if PLOTS:
    plt.close("all")
    plt.figure(figsize=[10,10])
    col = plt.imshow(chi_rho, cmap="RdBu", vmin=-1, vmax=1, extent=POLAR_EXTENT, aspect=ASPECT, origin="lower")
    plt.scatter(r_wake, phi_wake, c="k", marker='.', s=4)
    plt.colorbar(col, extend="both")
    plt.title("chi_dens")
    plt.show()
    plt.figure(figsize=[10,10])
    col = plt.imshow(chi_vr, cmap="RdBu", vmin=-1, vmax=1, extent=POLAR_EXTENT, aspect=ASPECT, origin="lower")
    plt.scatter(r_wake, phi_wake, c="k", marker='.', s=4)
    plt.colorbar(col, extend="both")
    plt.title("chi_vr")
    plt.show()
    plt.figure(figsize=[10,10])
    col = plt.imshow(chi_vphi, cmap="RdBu", vmin=-1, vmax=1, extent=POLAR_EXTENT, aspect=ASPECT, origin="lower")
    plt.scatter(r_wake, phi_wake, c="k", marker='.', s=4)
    plt.colorbar(col, extend="both")
    plt.title("chi_vphi")
    plt.show()

# === find slices at the radii we are probing === #

print("Extracting slices from given radii...")

# find indices where radii we selected occur in our grid
indices = []
for radius in CHI_RADII:
    for i, rad_sample in enumerate(r_reg):
        if radius <= rad_sample:
            indices.append(i)
            break

# find phi values along wake for selected radii
PHI_CHI = phi_wake_func(np.array(CHI_RADII))

# find indices for where to cut phi for each chi slice
phi_min_indices = []
phi_max_indices = []
phi_min_values  = []
phi_max_values  = []
for phi in PHI_CHI:
    min_phi = phi - (CHI_SLICE_WIDTH_ANG / 2)
    max_phi = phi + (CHI_SLICE_WIDTH_ANG / 2)
    for i, phi_sample in enumerate(phi_reg):
        if min_phi <= phi_sample:
            phi_min_indices.append(i)
            break
    for i, phi_sample in enumerate(phi_reg):
        if max_phi <= phi_sample:
            phi_max_indices.append(i)
            break
    phi_min_values.append(min_phi)
    phi_max_values.append(max_phi)

# create phi slices
phi_slices = []
for k in range(len(phi_min_indices)):
    i, j = phi_min_indices[k], phi_max_indices[k]
    phi_slices.append(phi_reg[i:j])

# create chi slices
chi_rho_slices  = []
chi_vr_slices   = []
chi_vphi_slices = []
for k, index in enumerate(indices):
    i, j = phi_min_indices[k], phi_max_indices[k]
    chi_rho_slices.append(chi_rho[i:j,index])
    chi_vr_slices.append(chi_vr[i:j,index])
    chi_vphi_slices.append(chi_vphi[i:j,index])

# === find t and phi values corresponding to the wake at the radii we are probing === #

print("Transforming slices to (t,eta) space...")

T_CHI = [t.t(radius, R_PLANET, HR_PLANET, Q_INDEX, P_INDEX) for radius in CHI_RADII]

TITLES = [r"$\chi (\rho)$", r"$\chi (v_r)$", r"$\chi (v_\phi)$"]

if CHI_2D_PLOT:

    fig, ax = plt.subplots(1, 3, figsize=[15,6], sharex='all', sharey='all')
    cw = ax[0].imshow(chi_rho,  cmap="RdBu", vmin=-1, vmax=1, extent=POLAR_EXTENT, aspect=ASPECT, origin="lower")
    ax[1].imshow(chi_vr,            cmap="RdBu", vmin=-1, vmax=1, extent=POLAR_EXTENT, aspect=ASPECT, origin="lower")
    ax[2].imshow(chi_vphi,          cmap="RdBu", vmin=-1, vmax=1, extent=POLAR_EXTENT, aspect=ASPECT, origin="lower")
    ax[0].set_ylabel(r"$\phi$ [rad]")
    for i in range(3):
        ax[i].scatter(r_wake, phi_wake, c="k", marker='.', s=2, alpha=0.2)
        ax[i].vlines(CHI_RADII, phi_min_values, phi_max_values, colors=['grey'], linestyles=['--'], alpha=0.5)
        ax[i].set_title(TITLES[i])
        ax[i].set_xlabel(r"$r$ [au]")
    plt.savefig("chis.pdf", dpi=200)
    plt.show()

# === find eta values for each slice === #

eta_slices = []
for i, radius in enumerate(CHI_RADII):
    phi_values = phi_slices[i]
    eta_vals = mod_pi_Eta(radius, phi_values)
    eta_slices.append(eta_vals)

# === extracting analytic form for chi initial condition === #

print("Extracting chi from analytic solution...")

# read perturbations from files and get arrays
analytic_perts = np.load("linear_perturbations.npy")
analytic_mesh  = np.load("linear_perturbations_mesh.npy")
v_r_an   = analytic_perts[0]
v_phi_an = analytic_perts[1]
rho_an   = analytic_perts[2]
X_an     = analytic_mesh[0]
Y_an     = analytic_mesh[1]

# linear perturbations read in grid
x = X_an[0,:]
y = Y_an[:,0]

# set height of cut
Y_MOD = 4

# cut square box grid in linear regime
x_cut = x[np.argmin(x <  BOX_SIZE_HILL)]
y_cut = y[np.argmin(y < -Y_MOD*BOX_SIZE_HILL) : np.argmin(y < Y_MOD*BOX_SIZE_HILL) + 1]

# convert x and y values to r and phi with phi in (-pi,pi)
phi_cut = np.arctan2(y_cut * L_SCALE, R_BOX)

# find cut indicies 
x_cut_i  = np.argmin(x <  BOX_SIZE_HILL)
y_cut_i1 = np.argmin(y < -Y_MOD*BOX_SIZE_HILL)
y_cut_i2 = np.argmin(y <  Y_MOD*BOX_SIZE_HILL) + 1

# cut perturbation arrays
cut_v_r   = v_r_an  [y_cut_i1:y_cut_i2, x_cut_i]
cut_v_phi = v_phi_an[y_cut_i1:y_cut_i2, x_cut_i]
cut_rho   = rho_an  [y_cut_i1:y_cut_i2, x_cut_i]

# scale perturbations to physical units
cut_v_r   *= CS_PLANET
cut_v_phi *= CS_PLANET

# scale perturbations by mass
cut_rho   *= BETA_P
cut_v_r   *= BETA_P
cut_v_phi *= BETA_P

# extract chi from density
chi_rho_an = cut_rho * ((GAMMA + 1) / 2) * t.g(R_BOX, R_PLANET, HR_PLANET, Q_INDEX, P_INDEX)

# extract chi from velocities
chi_vr_an   = cut_v_r   / (np.sign(R_BOX - R_PLANET) * t.Lambda_fu(R_BOX, R_PLANET, CS_PLANET, HR_PLANET, GAMMA, Q_INDEX, P_INDEX))
chi_vphi_an = cut_v_phi / (np.sign(R_BOX - R_PLANET) * t.Lambda_fv(R_BOX, R_PLANET, CS_PLANET, HR_PLANET, GAMMA, Q_INDEX, P_INDEX) * -CW)

# convert points to eta space
#eta_cut = mod_pi_Eta(R_BOX, phi_cut)
eta_cut = y_cut + 0.5 * x_cut**2 * np.sign(x_cut)

# === chi slice plots === #

if CHI_SLICE_PLOT:

    plt.figure(figsize=[10,10])
    plt.plot(eta_cut, chi_rho_an, lw=0.8, label=r"$\rho$", zorder=3)
    plt.plot(eta_cut, chi_vr_an, lw=0.8, label=r"$v_r$", zorder=2)
    plt.plot(eta_cut, chi_vphi_an, lw=0.8, label=r"$v_\phi$", zorder=1)
    plt.title(f"t = {round(T_BOX, 2)}, r = {round(R_BOX, 2)}")
    plt.xlabel(r"$\eta$")
    plt.ylabel(r"$\chi$")
    plt.legend(loc="upper left")
    plt.savefig(f"analytic_chi_t_{int(round(T_BOX,0))}.pdf")
    plt.show()

    N_SLICES = len(T_CHI)
    for i in range(N_SLICES):
        plt.figure(figsize=[10,10])
        plt.plot(eta_slices[i], chi_rho_slices[i],  lw=0.8, label=r"$\rho$", zorder=3)
        plt.plot(eta_slices[i], chi_vr_slices[i],   lw=0.8, label=r"$v_r$", zorder=2)
        plt.plot(eta_slices[i], chi_vphi_slices[i], lw=0.8, label=r"$v_\phi$", zorder=1)
        plt.title(f"t = {round(T_CHI[i], 2)}, r = {round(CHI_RADII[i], 2)}")
        plt.xlabel(r"$\eta$")
        plt.ylabel(r"$\chi$")
        plt.legend(loc="upper left")
        plt.savefig(f"chi_t_{int(round(T_CHI[i],0))}.pdf")
        plt.show()

    plt.figure(figsize=[10,10])
    plt.plot(eta_cut, chi_rho_an, lw=0.8, label=r"analytic")
    plt.plot(eta_slices[1], chi_rho_slices[1],  lw=0.8, label=r"sim")
    plt.title(f"t = {round(T_BOX, 2)}, r = {round(R_BOX, 2)}")
    plt.xlabel(r"$\eta$")
    plt.ylabel(r"$\chi$")
    plt.legend(loc="upper left")
    plt.show()