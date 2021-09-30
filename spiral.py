### Import Libraries 

import  yaml
import  numpy               as      np
import  matplotlib.pyplot   as      plt
import  pymcfost            as      pmc

from    matplotlib.colors   import  LogNorm
from    matplotlib          import  ticker, cm
from    casa_cube           import  casa_cube as casa



### Model parameters

# wakeflow parameter file
wakeflow_params = "HD_163296_secondary_kinks/smol_planet/config_hd163_v2.yaml"

# observations fits file (CO cube)
observations = "HD_163296_secondary_kinks/Observation/lines.fits"



### Plotting parameters

# angles in degrees for wakeflow model -- must match the observational ones!
inclination     = -225
PA              = 45
planet_az       = 45

# which channels to plot



### Read in wakeflow model

# import wakeflow model configuration
model_params    = yaml.load(open(wakeflow_params), Loader=yaml.FullLoader)

# check grid is Cartesian
if str(model_params["grid"]["type"]) != "cartesian":
    raise Exception("Must use wakeflow model with cartesian grid")
else:
    print("Reading wakeflow parameter file")
    
# check perturbations were saved
if not bool(model_params["results"]["save_perturbations"]):
  raise Exception("Must have saved perturbations")  
else:
    pass

# importing parameters
model_name      = str(model_params["run_info"]["name"])
model_system    = str(model_params["run_info"]["system"])

# check if multiple planet masses, ask which to use if multiple
try:
    model_m_planet = float(model_params["disk"]["m_planet"])
except:
    model_m_planet = list(model_params["disk"]["m_planet"])
    if len(model_m_planet) == 1:
        model_m_planet = model_m_planet[0]
    elif len(model_m_planet) > 1:
        print("Multiple planet masses detected: ", model_m_planet)
        num = input(f"Select model to use: (0-{len(model_m_planet)-1})")
        model_m_planet = model_m_planet[num]
    else:
        raise Exception("Planet mass missing")

# grab model directory
model_dir = f"{model_system}/{model_name}/{model_m_planet}Mj"

# read in arrays, and grab midplane -- resultant shape is (X, Y)
z_index = 0
X       = np.load(f"{model_dir}/X.npy")             [:,z_index,:]
Y       = np.load(f"{model_dir}/Y.npy")             [:,z_index,:]
v_r     = np.load(f"{model_dir}/delta_v_r.npy")     [:,z_index,:]
v_phi   = np.load(f"{model_dir}/delta_v_phi.npy")   [:,z_index,:]



### Angles for rotating velocity fields

# get angles needed in radians
PA          *= np.pi / 180.
planet_az   *= np.pi / 180.
inclination *= np.pi / 180.

# auxiliary angle
xi = np.arctan(np.tan(planet_az) / np.cos(inclination))
if planet_az == np.pi/2 or planet_az == 3*np.pi/2:
    xi = planet_az
elif planet_az > np.pi/2 and planet_az < 3*np.pi/2:
    xi += np.pi
elif planet_az > 3*np.pi/2:
    xi += 2*np.pi




###  get cartesian components of the velocity fields instead of radial and azimuthal

# get meshgrids for polar coordinates
R = np.sqrt(X**2 + Y**2)
PHI = np.arctan2(Y, X)

# perform transformations
v_x = -v_phi * np.sin(PHI) + v_r * np.cos(PHI)
v_y =  v_phi * np.cos(PHI) + v_r * np.sin(PHI)
v_z = np.zeros(v_x.shape)





"""
def rotate_meshgrid():

        for i in range(s.Nr):
            for j in range(s.Nphi):
                # rotation around disc normal axis, needed to match the required PAp
                s.X[j,i], s.Y[j,i], temp = np.dot( Rz(s.xi), np.array([s.X[j,i],s.Y[j,i],0]) )

                # rotation araound (sky plane) x-axis to get the required inclination!
                s.X[j,i], s.Y[j,i], temp = np.dot( Rx(s.i), np.array([s.X[j,i],s.Y[j,i],0]) )

                # rotation around disc normal axis to get the given PA!
                s.X[j,i], s.Y[j,i], temp = np.dot( Rz(s.PA), np.array([s.X[j,i],s.Y[j,i], 0]) )

        for i in range(200):
            s.disc_edge[i,0],s.disc_edge[i,1], temp = np.dot( Rx(s.i), np.array([s.disc_edge[i,0],s.disc_edge[i,1],0]) )
            s.disc_edge[i,0],s.disc_edge[i,1], temp = np.dot( Rz(s.PA), np.array([s.disc_edge[i,0],s.disc_edge[i,1],0]) )


def rotate_velocity_field(v_field):

    for i in range(s.Nr):
        for j in range(s.Nphi):

            # rotation around disc normal axis, needed to match the required PAp
            v_field[j,i,:] = np.dot( Rz(s.xi), v_field[j,i,:] )

            # rotation araound (sky plane) x-axis to get the required inclination!
            v_field[j,i,:] = np.dot( Rx(s.i), v_field[j,i,:] )

            # rotation around disc normal axis to get the given PA!
            v_field[j,i,:] = np.dot( Rz(s.PA), v_field[j,i,:] )


    return v_field


def get_channel_maps(vz_field):

    Nvch = len(s.vchs)
    channel_maps = 100 * np.ones(vz_field.shape) # Arbitrary large value
    for i in range(len(channel_maps[:,0])):
        for j in range(len(channel_maps[0,:])):
            for ch in range(Nvch):
                if (vz_field[i,j] > s.vchs[ch] - s.delta_vch and vz_field[i,j] < s.vchs[ch] + s.delta_vch):
                    channel_maps[i,j] = vz_field[i,j]

    return channel_maps
"""


### Functions for rotation matrices
def Rx(a):
   R = [[1, 0, 0],
        [0, np.cos(a),-np.sin(a)],
        [0, np.sin(a), np.cos(a)]
        ]
   return R

def Ry(b):  
   R = [[ np.cos(b), 0, np.sin(b)],
        [0, 1, 0],
        [-np.sin(b), 0, np.cos(b)]
        ]
   return R

def Rz(g):
   R = [[ np.cos(g),-np.sin(g), 0],
        [ np.sin(g), np.cos(g), 0],
        [0, 0, 1]
        ]
   return R