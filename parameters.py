import yaml
import sys
import numpy as np

class Constants:
    def __init__(self):

        # all constants given in CGS units
        self.m_solar = 1.988e33
        self.m_jupiter = 1.898e30

# load constants
c = Constants()

class Parameters:
    def __init__(self, config_file):

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

        # sanity checks
        m_thermal = 0.6666667 * self.hr_planet ** 3 * (self.m_star * c.m_solar / c.m_jupiter)
        print('mthermal = ', m_thermal)
        if self.m_planet > m_thermal:
            cont = input("Planet mass exceeds thermal mass. This may break the solution. Continue? [y/n]: ")
            if cont != "y":
                print('Exiting')
                sys.exit()
            else:
                print("Continuing")

Parameters('config.yaml')
