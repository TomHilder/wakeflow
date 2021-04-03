import yaml
import sys
import os
import shutil as sh
import numpy as np


def run_setup():
    """Runs initial setup. Returns Parameters object with parameters from config file.
    """

    # check that code has been loaded correctly
    if len(sys.argv) != 2:
        print("Error: Correct usage is '$python wakeflow.py config_file.yaml'")
        print('Exiting')
        sys.exit(1)

    # create Parameters object using configuration file
    config_file = sys.argv[1]
    params = Parameters(config_file)

    # do sanity checks and exit if not passed
    if params.do_sanity_checks() != True:
        print('Exiting')
        sys.exit(1)

    # create results directory
    path = params.name + '/'
    os.makedirs(path, exist_ok=True)

    # copy configuration file used to results directory
    sh.copy(config_file, path)

    return params


class Constants:
    def __init__(self):

        from astropy import constants as c

        # all constants given in CGS units
        self.m_solar = c.M_sun.cgs.value
        self.m_jupiter = c.M_jup.cgs.value
        self.G_const = c.G.cgs.value
        self.au = c.au.cgs.value


class Parameters(Constants):
    def __init__(self, config_file):

        print(f"Reading parameters from {config_file} ")

        # inherit constants
        super().__init__()

        # read in config file
        config = yaml.load(open(config_file), Loader=yaml.FullLoader)

        # run_info parameters
        self.name = str(config["run_info"]["name"])
        self.system = str(config["run_info"]["system"])
        self.date = str(config["run_info"]["date"])

        # disk parameters
        self.m_star = float(config["disk"]["m_star"])
        self.m_planet = float(config["disk"]["m_planet"])
        self.r_outer = float(config["disk"]["r_outer"])
        self.r_inner = float(config["disk"]["r_inner"])
        self.r_planet = float(config["disk"]["r_planet"])
        self.r_ref = float(config["disk"]["r_ref"])
        self.q = float(config["disk"]["q"])
        self.p = float(config["disk"]["p"])
        self.hr = float(config["disk"]["hr"])
        self.cw_rotation = bool(config["disk"]["cw_rotation"])

        # angles parameters
        self.inclination = float(config["angles"]["inclination"])
        self.PA = float(config["angles"]["PA"])
        self.PAp = float(config["angles"]["PAp"])

        # grid parameters
        self.grid_type = str(config["grid"]["type"])
        self.n_x = int(config["grid"]["n_x"])
        self.n_y = int(config["grid"]["n_y"])
        self.n_r = int(config["grid"]["n_r"])
        self.n_phi = int(config["grid"]["n_phi"])
        self.r_log = bool(config["grid"]["r_log"])
        self.make_3D = bool(config["grid"]["make_3D"])

        # plot parameters
        self.make_plots = bool(config["plotting"]["make_plots"])
        self.show_plots = bool(config["plotting"]["show_plots"])

        # mcfost parameters
        self.make_cube = bool(config["mcfost"]["make_cube"])
        self.run_mcfost = bool(config["mcfost"]["run_mcfost"])
        self.pymcfost_plots = bool(config["mcfost"]["pymcfost_plots"])

        # physical parameters
        self.gamma = float(config["physical"]["adiabatic_index"])
        self.malpha = float(config["physical"]["damping_malpha"])

        # numerical parameters
        self.CFL = float(config["numerical"]["CFL"])
        self.scale_box = float(config["numerical"]["scale_box"])

        # get flaring at r_planet
        self.hr_planet = self.hr * (self.r_planet / self.r_ref) ** (0.5 - self.q)

    def do_sanity_checks(self):

        # check that planet mass does not exceed thermal mass
        m_thermal = 0.6666667 * self.hr_planet ** 3 * (self.m_star * self.m_solar / self.m_jupiter)
        #print('Thermal mass = ', m_thermal)
        if self.m_planet > m_thermal:
            if not self.warning("Planet mass exceeds thermal mass. This may break the solution."):
                return False

        # check grid type
        if self.grid_type != "cartesian" and self.grid_type != "cylindrical":
            print("Error: Please choose a valid grid type (cartesian or cylindrical)")
            return False
        
        # check settings OK if mcfost is to be run
        if self.run_mcfost:
            if not self.make_cube:
                print("Error: You must select make_cube = True to run mcfost")
                return False
            if not self.r_log:
                print("Error: You must select r_log = True to run mcfost")
                return False
            if self.grid_type != "cylindrical":
                print("Error: You must use a cylindrical grid to run mcfost")
                return False

        # check linear box scale factor
        if self.scale_box != 1:
            if not self.warning("Changing linear box scale factor can cause strange results."):
                return False

        # check CFL
        if self.CFL > 0.5:
            if not self.warning("CFL chosen > 0.5, this will likely break the numerical PDE solver."):
                return False

        print("Parameters Ok. Continuing... ")
        return True
      
    def warning(self, warning_msg):
        statement = "Warning: " + warning_msg + " Continue? [y/n]: "
        cont = input(statement)
        if cont != "y":
            return False
        else:
            return True
