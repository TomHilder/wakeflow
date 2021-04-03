from setup import Constants
import numpy as np

c = Constants()

class Grid:
    def __init__(self, parameters):

        # grab parameters object
        self.p = parameters

        # initialise grid properties
        self.info = {
            "Type": None,
            "Log_r": None,
            "Size": [0, 0, 0],
            "Contains": "Empty"
            }

        # initialise arrays
        self.v_r = None 
        self.v_phi = None 
        self.rho_pert = None

    def make_keplerian_disk(self):

        print("Making Keplerian disk ")

        # get radii
        if self.info["Type"] == "cartesian":
            self.get_r_phi_coords()
            _r = self.R_xy
        else:
            _r = self.R

        # Keplerian velocities
        self.v_r = np.zeros((self.info["Size"][0], self.info["Size"][1]))
        self.v_phi = np.sqrt(c.G_const * self.p.m_star * self.p.m_solar / (_r * self.p.au))

        # pressure correction
        _corr = (1 - (2*self.p.q + self.p.p) * self.p.hr**2 * (_r / self.p.r_ref)**(1 - 2 * self.p.q))**(1 / 2)
        self.v_phi *= _corr

        # update grid info
        self.info["Contains"] = "Keplerian velocity field"

    def make_grid(self):

        print(f"Constructing {self.p.grid_type} Grid ")

        # make cartesian grid
        if self.p.grid_type == "cartesian":
            self.make_cartesian_grid()

        # make cylindrical grid
        elif self.p.grid_type == "cylindrical":
            self.make_cylindrical_grid()


    def make_cartesian_grid(self):
        
        # make grid from specifications in parameter file
        self.x = np.linspace(-self.p.r_outer, self.p.r_outer, self.p.n_x)
        self.y = np.linspace(-self.p.r_outer, self.p.r_outer, self.p.n_y)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        # update grid info
        self.info["Type"] = "cartesian"
        self.info["Size"][0] = self.x.shape[0]
        self.info["Size"][1] = self.y.shape[0]

    def make_cylindrical_grid(self):
         
        # make grid from specifications in parameter file
        if self.p.r_log:
            self.r = np.geomspace(self.p.r_inner, self.p.r_outer, self.p.n_r)
        else:
            self.r = np.linspace(self.p.r_inner, self.p.r_outer, self.p.n_r)    
        self.phi = np.linspace(0, 2*np.pi, self.p.n_phi)
        self.R, self.PHI = np.meshgrid(self.r, self.phi)

        # update grid info
        self.info["Type"] = "cylindrical"
        self.info["Size"][0] = self.r.shape[0]
        self.info["Size"][1] = self.phi.shape[0]
        self.info["Log_r"] = self.p.r_log

    def get_r_phi_coords(self):

        # often need the equivalent (r,phi) coordinates when using cartesian grid
        self.r_xy = np.sqrt(self.x**2 + self.y**2)
        self.phi_xy = np.arctan(self.y / self.x)
        self.R_xy, self.PHI_xy = np.meshgrid(self.r_xy, self.phi_xy)

    def merge_linear(self):
        pass

    def merge_nonlinear(self):
        pass

    def extrapolate_to_3D_grid(self):
        pass