import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
import os, sys, subprocess
from pymcfost.parameters import Params
from pymcfost.disc_structure import Disc

def make_mcfost_parameter_file(parameters):
    
    # grab run parameters
    p = parameters

    # setting up an mcfost parameters object
    mp = Params("mcfost/ref3.0_3D.para")

    # --- mcfost parameters that are not set by config.yaml

    mp.phot.nphot_T = 1.28e+07          # set number of photons for temp. calculation
    mp.simu.compute_SED = False         # don't compute SED

    mp.map.nx = 1001                    # output resolution
    mp.map.ny = 1001

    mp.simu.image_symmetry = False      # we don't have any symmetries
    mp.simu.central_symmetry = False
    mp.simu.axial_symmetry = False

    mp.mol.molecule[0].n_trans = 1      # unsure
    mp.grid.n_rad_in = 1                # unsure

    mp.mol.v_turb = 0.05                # turbulence

    mp.zones[0].edge = 0.0              # not using smoothed inner and outer edge
    mp.zones[0].Rc = 0.0                # critical radius unused so set to 0

    mp.map.RT_n_az = 1                  # we are only using one azimuth
    mp.map.RT_ntheta = 1                # we are only using one inclination

    # --- mcfost parameters that are set by config.yaml

    mp.grid.n_rad = p.n_r               # grid geometry
    mp.grid.n_az = p.n_phi
    mp.grid.nz = p.n_z
    mp.zones[0].Rin = p.r_inner
    mp.zones[0].Rout = p.r_outer

    mp.zones[0].h0 = p.h_ref            # scale height at reference radius
    mp.zones[0].Rref = p.r_ref          # reference radius
    mp.zones[0].flaring_exp = p.beta    # flaring index: h \propto r^beta

    mp.map.distance = p.distance        # distance

    mp.stars[0].Teff = p.temp           # star temperature

    #mp.map.PA = p.PA - 90               # PA angle
    mp.map.PA = p.PA

    #mp.map.RT_imin = -1*p.inclination   # inclination
    #mp.map.RT_imax = -1*p.inclination
    mp.map.RT_imin = p.inclination   # inclination
    mp.map.RT_imax = p.inclination

    #mp.map.RT_az_min = -1*p.PAp         # azimuth used to rotate to correct planet position
    #mp.map.RT_az_max = -1*p.PAp
    mp.map.RT_az_min = p.PAp         # azimuth used to rotate to correct planet position
    mp.map.RT_az_max = p.PAp

    mp.mol.molecule[0].v_max = p.v_max  # max velocity
    mp.mol.molecule[0].nv = p.n_v       # number of velocity channels

    # write mcfost parameter file
    mp.writeto(f"{p.system}/{p.name}/mcfost/mcfost_{p.name}.para")

def make_mcfost_grid_data(parameters):

    print("Generating MCFOST grid data...")

    # grab run parameters
    p = parameters

    # generate grid data by running mcfost (in mcfost directory)
    working_dir = os.getcwd()
    os.chdir(f"{p.system}/{p.name}/mcfost/")
    subprocess.call(["rm", "-rf", "data_disk", "data_disk_old"])
    subprocess.call(["mcfost", f"mcfost_{p.name}.para", "-disk_struct"], stdout=subprocess.DEVNULL)
    os.chdir(working_dir)

    print("Done")

def read_mcfost_grid_data(parameters):

    # grab run parameters
    p = parameters
    
    # reading disk data
    with HiddenPrints():
        mcfost_disk = Disc(f"./{p.system}/{p.name}/mcfost/")

    # getting radii and heights
    r = mcfost_disk.r()[:, p.n_z:, :]
    z = mcfost_disk.z()[:, p.n_z:, :]

    return r, z

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout