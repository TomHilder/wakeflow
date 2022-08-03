import yaml, sys, os
import shutil   as sh
import numpy    as np

def load_config_file(config_file, default_config_dict=None):

    # read in config file as dictionary
    config_dict = yaml.load(open(config_file), Loader=yaml.FullLoader)

    if default_config_dict is not None:
        for key in config_dict.keys():
            if key not in default_config_dict.keys():
                raise Exception(f"{key} is not a valid parameter.")

    return config_dict

def write_config_file(config_dict, directory, filename):

    with open(f'{directory}{filename}', 'w') as yaml_file:
        yaml.dump(config_dict, yaml_file, default_flow_style=False)

def run_setup(param_dict, default_param_dict=None, overwrite=False):

    if default_param_dict is not None:
        for key in param_dict.keys():
            if key not in default_param_dict.keys():
                raise Exception(f"{key} is not a valid parameter.")

    params = Parameters(param_dict)

    # do sanity checks and exit if not passed
    params.do_sanity_checks()

    # check if directory for system exists, if not create it
    system_path = f"{params.system}/"
    os.makedirs(system_path, exist_ok=True)

    # create results directory path
    results_path = f"{params.system}/{params.name}/"

    # check if results directory already exists, ask to overwrite if yes
    results_exist = os.path.isdir(results_path)
    if results_exist is True:
        if overwrite:
            print("Overwriting previous results")
            sh.rmtree(results_path)
            os.makedirs(results_path)
        else:
            raise Exception("Results already exist for run name. Either choose a different name, or run with overwrite=True")

    # create individual directories for each planet mass result
    if params.m_planet_array is not None:
        for mass in params.m_planet_array:
            individual_result_path = f"{params.system}/{params.name}/{mass}Mj"
            os.makedirs(individual_result_path, exist_ok=True)
    else:
        individual_result_path = f"{params.system}/{params.name}/{params.m_planet}Mj"
        os.makedirs(individual_result_path, exist_ok=True)

    # write parameters used to a file in the results directory
    write_config_file(param_dict, results_path, f"{params.name}_config.yaml")

    # run mcfost grid setup if needed
    if params.run_mcfost or params.grid_type == "mcfost":

        from .mcfost_interface import make_mcfost_parameter_file, make_mcfost_grid_data

        # make directory for mcfost outputs
        mcfost_path = f"{params.system}/{params.name}/mcfost/"
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
        self.m_solar    = c.M_sun.cgs.value
        self.m_jupiter  = c.M_jup.cgs.value
        self.G_const    = c.G.cgs.value
        self.au         = c.au.cgs.value


class Parameters(Constants):
    def __init__(self, config):

        #print(f"Reading parameters from {config_file} \n ")

        # inherit constants
        super().__init__()

        # read in config file
        #config = yaml.load(open(config_file), Loader=yaml.FullLoader)

        # run_info parameters
        self.name   = str(config["name"])
        self.system = str(config["system"])

        # disk parameters
        self.m_star = float(config["m_star"])

        try:
            self.m_planet       = float(config["m_planet"])
            self.m_planet_array = None
        except:
            self.m_planet       = None
            self.m_planet_array = list(config["m_planet"])

        self.r_outer  = float(config["r_outer"])
        self.r_inner  = float(config["r_inner"])
        self.r_planet = float(config["r_planet"])
        self.r_ref    = float(config["r_ref"])
        self.r_c      = float(config["r_c"])
        self.q        = float(config["q"])
        self.p        = float(config["p"])
        self.hr       = float(config["hr"])
        self.rho_ref  = float(config["dens_ref"])

        # override given rotation, only calculated anticlockwise
        self.cw_rotation      = False
        self.user_cw_rotation = bool(config["cw_rotation"])

        # angles parameters
        self.inclination = float(config["inclination"])
        self.PA          = float(config["PA"])
        self.PAp         = float(config["PAp"])

        # grid parameters
        self.grid_type = str (config["type"])
        self.n_x       = int (config["n_x"])
        self.n_y       = int (config["n_y"])
        self.n_r       = int (config["n_r"])
        self.n_phi     = int (config["n_phi"])
        self.n_z       = int (config["n_z"])
        self.r_log     = bool(config["r_log"])

        # plot parameters
        self.make_midplane_plots   = bool(config["make_midplane_plots"])
        self.show_midplane_plots   = bool(config["show_midplane_plots"])
        self.show_teta_debug_plots = bool(config["show_teta_debug_plots"])

        # results parameters
        self.use_planet         = bool(config["include_planet"])
        self.include_linear     = bool(config["include_linear"])
        self.save_perturbations = bool(config["save_perturbations"])
        self.save_total         = bool(config["save_total"])
        self.write_FITS         = bool(config["write_FITS"])
        self.dimensionless      = bool(config["dimensionless"])

        # mcfost parameters
        self.run_mcfost = bool (config["run_mcfost"])
        self.temp       = float(config["temp_star"])
        self.distance   = float(config["distance"])
        self.v_max      = float(config["v_max"])
        self.n_v        = int  (config["n_v"])

        # pymcfost parameters
        #self.pymcfost_plots = bool(config["pymcfost"]["pymcfost_plots"])
        #self.velocity_channels = list(config["pymcfost"]["velocity_channels"])
        #self.beam = float(config["pymcfost"]["beam"])
        #self.obs_dir = str(config["pymcfost"]["obs_dir"])
        #self.sim_dir = str(config["pymcfost"]["sim_dir"])
        #self.v_system = float(config["pymcfost"]["v_system"])

        # physical parameters
        self.gamma  = float(config["adiabatic_index"])
        self.malpha = float(config["damping_malpha"])

        # numerical parameters
        self.CFL           = float(config["CFL"])
        self.scale_box     = float(config["scale_box"])
        self.scale_box_ang = float(config["scale_box_ang"])
        self.box_warp      = bool (config["box_warp"])
        self.use_box_IC    = bool (config["use_box_IC"])
        self.tf_fac        = float(config["tf_fac"])

        # get flaring at r_planet
        self.hr_planet = self.hr * (self.r_planet / self.r_ref) ** (0.5 - self.q)

        # get length scale l
        self.l = (2/3) * self.hr_planet * self.r_planet

        # get sound speed at planet radius and reference radius
        v_kep_planet    = np.sqrt(self.G_const*self.m_star*self.m_solar / (self.r_planet*self.au))
        v_kep_ref       = np.sqrt(self.G_const*self.m_star*self.m_solar / (self.r_ref*self.au))
        self.c_s_planet = v_kep_planet*self.hr_planet
        self.c_s_0      = v_kep_ref*self.hr 

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

        # get dimensions
        if self.dimensionless:
            self.U_len  = self.r_ref * self.au
            self.U_mass = self.m_solar
            self.U_time = np.sqrt(self.U_len**3 / (self.G_const * self.U_mass)) 
            self.U_vel  = self.U_len / self.U_time / 1E5

    def do_sanity_checks(self):

        print("\n* Performing sanity checks on model parameters:")

        if self.m_planet != None:
            print(f"M_thermal = {self.m_thermal:.3f} M_Jup")
            print(f"M_planet  = {(self.m_planet / self.m_thermal):.3f} M_th")

            # check that planet mass does not exceed thermal mass
            if self.m_planet > self.m_thermal:
                print("WARNING: M_planet > M_thermal. This may break the solution.")
                
        else:
            print(f"M_thermal = {self.m_thermal:.3f} M_Jup")
            print(f"M_planets = {np.array2string(np.array(self.m_planet_array)/self.m_thermal, precision=3, floatmode='fixed')} M_th")
            
            
            #print(f"Planet masses are: ", end='')
            #print(np.array(self.m_planet_array)/self.m_thermal, end='')
            #print(" thermal masses")

            # check that planet mass does not exceed thermal mass
            exceed = False
            for mass in self.m_planet_array:
                if mass > self.m_thermal:
                    exceed = True
            if exceed ==  True:
                print("At least one planet mass exceeds thermal mass. This may break the solution. ")

        # check grid type
        if self.grid_type != "cartesian" and self.grid_type != "cylindrical" and self.grid_type != "mcfost":
            raise Exception("Invalid grid type. Choose either cartesian or cylindrical or mcfost)")
        
        # check settings OK if mcfost is to be run
        if self.run_mcfost:
            if self.grid_type != "mcfost":
                raise Exception("Cannot run mcfost without using mcfost grid")
            elif self.write_FITS != True:
                raise Exception("Cannot run mcfost without writing FITS file (ie. require write_FITS: True)")

        # check linear box scale factor
        if self.scale_box != 1:
            print("WARNING: Changing linear box scale factor can cause strange results.")

        # check CFL
        if self.CFL > 0.5:
            print("WARNING: CFL chosen > 0.5, this will likely break the numerical PDE solver.")

        # check box warp
        if not self.box_warp:
            print("WARNING: Choosing box_warp=False may invalidate your results.")

        # check IC read out
        if self.use_box_IC:
            print("WARNING: Choosing use_box_IC=True will almost certainly invalidate your results.")

        # debug plots t eta
        if self.show_teta_debug_plots:
            print("WARNING: Choosing show_teta_debug_plots=True may cause the run to fail.")

        print("Parameters Ok - continuing")
        return True
