# model_setup.py
# Written by Thomas Hilder

"""
Contains the Parameters class responsible for storing the parameters used in a Wakeflow model. Additionally, contains
functions used read in said parameters, check their validity, and create/replace the results directory.
"""

import yaml, os
import shutil   as sh
import numpy    as np

# NOTE: contents are intended for internal use and should not be directly accessed by users

# class for storing physical constants, will be inherited by _Parameters class
class _Constants:
    """
    Constants class for storing physical constants needed by Wakeflow. Inherited by _Parameters.
    """

    def __init__(self):

        from astropy import constants as c

        # all constants given in CGS units
        self.m_solar    = c.M_sun.cgs.value
        self.m_jupiter  = c.M_jup.cgs.value
        self.G_const    = c.G.cgs.value
        self.au         = c.au.cgs.value

# class for storing parameters to be used in a Wakeflow model
class _Parameters(_Constants):
    """
    Class for storing Wakeflow model parameters to be easily handed around by various parts of Wakeflow.
    """

    # read in parameters from dictionary to class attributes, calculate some needed quantities
    def __init__(self, config: dict) -> None:
        """Inherit constants, read in parameters to dictionary, calculate quantities like thermal mass and scale height.
        """

        # inherit constants
        super().__init__()

        # run_info parameters
        self.name   = str(config["name"])
        self.system = str(config["system"])

        # disk parameters
        self.m_star = float(config["m_star"])

        # user can provide either single planet mass or list of planet masses
        try:
            self.m_planet       = float(config["m_planet"])
            self.m_planet_array = None
        except:
            self.m_planet       = None
            self.m_planet_array = list(config["m_planet"])

        self.r_outer    = float(config["r_outer"])
        self.r_inner    = float(config["r_inner"])
        self.r_planet   = float(config["r_planet"])
        self.phi_planet = float(config["phi_planet"])
        try:
            self.r_ref = float(config["r_ref"])
        except:
            self.r_ref = self.r_planet
        self.r_c      = float(config["r_c"])
        self.z_max    = float(config["z_max"])
        self.q        = float(config["q"])
        self.dens_p   = float(config["p"])
        self.p        = self.dens_p + self.q - 1.5
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
        self.grid_type = str (config["grid_type"])
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

        # physical parameters
        self.gamma  = float(config["adiabatic_index"])
        self.malpha = float(config["damping_malpha"])

        # numerical parameters
        self.CFL                  = float(config["CFL"])
        self.smooth_box           = float(config["smooth_box"])
        self.scale_box_left       = float(config["scale_box_left"])
        self.scale_box_right      = float(config["scale_box_right"])
        self.scale_box_ang_top    = float(config["scale_box_ang_top"])
        self.scale_box_ang_bottom = float(config["scale_box_ang_bottom"])
        self.box_warp             = bool (config["box_warp"])
        self.use_box_IC           = bool (config["use_box_IC"])
        self.tf_fac               = float(config["tf_fac"])
        self.r_cut_inner_fac      = float(config["r_cut_inner_fac"])
        
        # Choice of physics
        self.use_old_vel          = bool(config["use_old_vel"])
        self.lin_type             = str (config["lin_type"])
        self.nl_wake              = bool(config["nl_wake"])
        self.vr_evolution         = bool(config["vr_evolution"])
        
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
        #if self.cw_rotation == True:
        #    self.a_cw = -1
        #else:
        self.a_cw = 1
        
        # check if phi planet is not zero
        if self.phi_planet != 0.:
            self.rot_interp = True
        else:
            self.rot_interp = False

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

    # checks that the combination of parameters the user has provided is sensical
    def _do_sanity_checks(self) -> bool:
        """Called to check that the parameters specified by the user will not break the results/run. Some parameter combinations are incompatible.
        """

        print("\n* Performing checks on model parameters:")

        if self.m_planet != None:
            print(f"M_thermal = {self.m_thermal:.3f} M_Jup")
            print(f"M_planet  = {(self.m_planet / self.m_thermal):.3f} M_th")

            # check that planet mass does not exceed thermal mass
            if self.m_planet > self.m_thermal:
                print("WARNING: M_planet > M_thermal. This may break the solution.")
                
        else:
            print(f"M_thermal = {self.m_thermal:.3f} M_Jup")
            print(f"M_planets = {np.array2string(np.array(self.m_planet_array)/self.m_thermal, precision=3, floatmode='fixed')} M_th")

            # check that planet mass does not exceed thermal mass
            exceed = False
            for mass in self.m_planet_array:
                if mass > self.m_thermal:
                    exceed = True
            if exceed ==  True:
                print("WARNING: At least one planet mass exceeds thermal mass. This may break the solution. ")

        # check grid type
        if self.grid_type != "cartesian" and self.grid_type != "cylindrical" and self.grid_type != "mcfost":
            raise Exception("Invalid grid type. Choose either cartesian or cylindrical or mcfost)")

        # check settings OK if mcfost is to be run
        if self.run_mcfost:
            if self.grid_type != "mcfost":
                raise Exception("Cannot run mcfost without using mcfost grid")
            elif self.write_FITS != True:
                raise Exception("Cannot run mcfost without writing FITS file (ie. require write_FITS: True)")

        # check for planet rotation
        if self.rot_interp is True and self.grid_type != "cartesian":
            raise Exception("Currently you must choose grid_type='cartesian' to use non-zero phi_planet.")

        # check if box smoothing is enabled
        if self.smooth_box:
            print("WARNING: Using smooth_box=True can cause strange results.")
            if self.grid_type == "cylindrical":
                raise Exception("You must choose grid_type='cartesian' or 'mcfost' to use smooth_box=True")

        # check linear box scale factor
        if self.scale_box_left != 1 or self.scale_box_right != 1:
            print("WARNING: Changing linear box scale factor can cause strange results.")

        # check linear box scale factor
        if self.scale_box_left != self.scale_box_right:
            print("WARNING: Using a different linear box scale factor for left and right edge can cause strange results.")

        # check linear box scale factor
        if self.scale_box_ang_top != self.scale_box_ang_bottom:
            print("WARNING: Using a different linear box scale factor for top and bottom edge can cause strange results.")

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

        # check velocity formulas
        if self.use_old_vel:
            print("WARNING: Choosing use_old_vel=True may cause a different velocity output.")
        print("Parameters Ok - continuing")
        return True
            
        # check linear perturbations input from user    
        if self.lin_type != "shearing_sheet" and self.lin_type != "global" and self.lin_type != "simulation":
            raise Exception("Invalid linear perturbation type. Choose either global or simulation or shearing_sheet)")
            
        # Warning on shearing sheet approximation
        if self.lin_type == "shearing_sheet":
            print("WARNING: The shearing sheet approximation may be invalid for a given choice of parameters. This may lead to incorrect results")

 # read in .yaml file, check keys correspond to parameters and return dictionary of parameters
def _load_config_file(config_file: str, default_config_dict: dict = None) -> dict:
    """Reads .yaml parameter file into a dictionary so that Wakeflow may parse it.
    """

    # read in config file as dictionary
    with open(config_file, "r") as file:
        config_dict = yaml.load(file, Loader=yaml.FullLoader)

    # check that keys in config_dict correspond to actual Wakeflow parameters
    if default_config_dict is not None:
        for key in config_dict.keys():
            if key not in default_config_dict.keys():
                raise Exception(f"{key} is not a valid parameter.")
            
    

    return config_dict

# write dict to .yaml file
def _write_config_file(config_dict: dict, directory: str, filename: str) -> None:
    """Writes parameters dictionary to .yaml file.
    """

    with open(f'{directory}{filename}', 'w') as yaml_file:
        yaml.dump(config_dict, yaml_file, default_flow_style=False)

def _run_setup(param_dict: dict, overwrite: bool = False) -> _Parameters:
    """Perform setup by generating parameters object, checking the parameters are okay, creating results directory,
    writing .yaml file with parameters used to results dir, and calling MCFOST to generate grid data if needed.
    """

    params = _Parameters(param_dict)

    # do sanity checks and exit if not passed
    params._do_sanity_checks()

    # check if anything is being saved, if yes create directory
    if params.save_perturbations or params.save_total or params.make_midplane_plots:

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
        _write_config_file(param_dict, results_path, f"{params.name}_config.yaml")

    else:

        pass

    # run mcfost grid setup if needed
    if params.run_mcfost or params.grid_type == "mcfost":

        from .mcfost_interface import _make_mcfost_parameter_file, _make_mcfost_grid_data

        # make directory for mcfost outputs
        mcfost_path = f"{params.system}/{params.name}/mcfost/"
        os.makedirs(mcfost_path, exist_ok=True)

        # generate mcfost parameter file
        _make_mcfost_parameter_file(params)

        # generate mcfost grid data to run analytics on
        _make_mcfost_grid_data(params)

    return params


