import yaml, sys, os
import shutil as sh
import numpy as np
from mcfost_interface import make_mcfost_parameter_file, make_mcfost_grid_data, read_mcfost_grid_data

def run_setup():

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

    # check if directory for system exists, if not create it
    system_path = f"{params.system}/"
    os.makedirs(system_path, exist_ok=True)

    # create results directory
    results_path = f"{params.system}/{params.name}/"
    os.makedirs(results_path, exist_ok=True)

    # copy configuration file used to results directory
    sh.copy(config_file, results_path)

    # run mcfost grid setup if needed
    if params.run_mcfost or params.grid_type == "mcfost":

        # make directory for mcfost outputs
        mcfost_path = f"{params.system}/{params.name}/mcfost_output/"
        os.makedirs(mcfost_path, exist_ok=True)

        # generate mcfost parameter file
        make_mcfost_parameter_file(params)

        # generate mcfost grid data to run analytics on
        make_mcfost_grid_data(params)

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
        self.rho_ref = float(config["disk"]["dens_ref"])
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
        self.n_z = int(config["grid"]["n_z"])
        self.r_log = bool(config["grid"]["r_log"])

        # plot parameters
        self.make_plots = bool(config["plotting"]["make_plots"])
        self.show_plots = bool(config["plotting"]["show_plots"])
        self.synthetic_velocity_channels = list(config["plotting"]["synthetic_velocity_channels"])

        # mcfost parameters
        self.run_mcfost = bool(config["mcfost"]["run_mcfost"])
        self.temp = float(config["mcfost"]["temp_star"])
        self.distance = float(config["mcfost"]["distance"])
        self.v_max = float(config["mcfost"]["v_max"])
        self.n_v = int(config["mcfost"]["n_v"])
        self.pymcfost_plots = bool(config["mcfost"]["pymcfost_plots"])
        self.velocity_channels = list(config["mcfost"]["velocity_channels"])

        # physical parameters
        self.gamma = float(config["physical"]["adiabatic_index"])
        self.malpha = float(config["physical"]["damping_malpha"])

        # numerical parameters
        self.CFL = float(config["numerical"]["CFL"])
        self.scale_box = float(config["numerical"]["scale_box"])
        self.scale_box_ang = float(config["numerical"]["scale_box_ang"])

        # get flaring at r_planet
        self.hr_planet = self.hr * (self.r_planet / self.r_ref) ** (0.5 - self.q)

        # get length scale l
        self.l = (2/3) * self.hr_planet * self.r_planet

        # get sound speed at planet radius
        v_kep_planet = np.sqrt(self.G_const*self.m_planet*self.m_solar / (self.r_planet*self.au))
        self.c_s_planet = v_kep_planet*self.hr_planet

        # clockwise rotation factor
        if self.cw_rotation == True:
            self.a_cw = -1
        else:
            self.a_cw = 1

        # get height scale at reference radius
        self.h_ref = self.hr * self.r_ref

        # get flaring exponent beta (exponent for h, NOT h/r) for mcfost
        self.beta = 1 + 0.5 - self.q

        # get thermal mass
        self.m_thermal = 0.6666667 * self.hr_planet ** 3 * (self.m_star * self.m_solar / self.m_jupiter)


    def do_sanity_checks(self):

        # check that planet mass does not exceed thermal mass
        if self.m_planet > self.m_thermal:
            if not self.warning("Planet mass exceeds thermal mass. This may break the solution."):
                return False

        # check grid type
        if self.grid_type != "cartesian" and self.grid_type != "cylindrical" and self.grid_type != "mcfost":
            print("Error: Please choose a valid grid type (cartesian or cylindrical or mcfost)")
            return False
        
        # check settings OK if mcfost is to be run  -- NEEDS TO BE UPDATED
        if self.run_mcfost:
            if self.grid_type != "mcfost":
                print("Error: You must use mcfost grid to run mcfost")
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
        statement = f"Warning: {warning_msg} Continue? [y/n]: "
        cont = input(statement)
        if cont != "y":
            return False
        else:
            return True
