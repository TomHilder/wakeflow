import yaml
import sys
import numpy as np

class Constants:
    def __init__(self):

        # all constants given in CGS units
        self.m_solar = 1.988e33
        self.m_jupiter = 1.898e30


class Parameters(Constants):
    def __init__(self, config_file):

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
        self.generate_3D = bool(config["grid"]["generate_3D"])

        # mcfost parameters
        self.generate_cube = bool(config["mcfost"]["generate_cube"])
        self.run_mcfost = bool(config["mcfost"]["run_mcfost"])

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
            if not self._warning("Planet mass exceeds thermal mass. This may break the solution."):
                return False

        # check grid type
        if self.grid_type != "cartesian" and self.grid_type != "cylindrical":
            print("Error: Please choose a valid grid type (cartesian or cylindrical)")
            return False
        
        # check settings OK if mcfost is to be run
        if self.run_mcfost:
            if not self.generate_cube:
                print("Error: You must select generate_cube = True to run mcfost")
                return False
            if not self.r_log:
                print("Error: You must select r_log = True to run mcfost")
                return False
            if self.grid_type != "cylindrical":
                print("Error: You must use a cylindrical grid to run mcfost")
                return False

        # check linear box scale factor
        if self.scale_box != 1:
            if not self._warning("Changing linear box scale factor can cause strange results."):
                return False

        # check CFL
        if self.CFL > 0.5:
            if not self._warning("CFL chosen > 0.5, this will likely break the numerical PDE solver."):
                return False

        print("Continuing")
        return True
      
    def _warning(self, warning_msg):
        statement = "Warning: " + warning_msg + " Continue? [y/n]: "
        cont = input(statement)
        if cont != "y":
            return False
        else:
            return True

            
        



p = Parameters('config.yaml')
if p.do_sanity_checks() != True:
    print('Exiting')
    sys.exit()
