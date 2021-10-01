# Import Libraries

import yaml
import numpy as np
import matplotlib.pyplot as plt
import pymcfost as pmc

from matplotlib.colors import LogNorm
from matplotlib import ticker, cm
from casa_cube import casa_cube as casa


def main():
    ### Model parameters

    # wakeflow parameter file
    wakeflow_params = "HD_163296_secondary_kinks/smol_planet/config_hd163_v2.yaml"

    # observations fits file (CO cube)
    observations = "HD_163296_secondary_kinks/Observation/lines.fits"

    # use perturbations or total velocity? "delta" or "total"
    v_type = "delta"

    ### Plotting parameters

    # angles in degrees for wakeflow model -- must match the observational ones!
    inclination = -225
    PA = 45
    planet_az = 45

    # which channels to plot


    ### Read in wakeflow model

    # import wakeflow model configuration
    model_params = yaml.load(open(wakeflow_params), Loader=yaml.FullLoader)

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
    model_name = str(model_params["run_info"]["name"])
    model_system = str(model_params["run_info"]["system"])

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
    X = np.load(f"{model_dir}/X.npy")[:, z_index, :]
    Y = np.load(f"{model_dir}/Y.npy")[:, z_index, :]
    v_r = np.load(f"{model_dir}/{v_type}_v_r.npy")[:, z_index, :]
    v_phi = np.load(f"{model_dir}/{v_type}_v_phi.npy")[:, z_index, :]

    plt.contourf(X, Y, 1e3*v_r, levels=np.linspace(-200,200,199), cmap="RdBu")
    plt.title("v_r")
    plt.colorbar()
    plt.show()
    
    plt.contourf(X, Y, 1e3*v_phi, levels=np.linspace(-200,200,199), cmap="RdBu")
    plt.title("v_phi")
    plt.colorbar()
    plt.show()

    ### Angles for rotating velocity fields

    # get angles needed in radians
    PA *= np.pi / 180.
    planet_az *= np.pi / 180.
    inclination *= np.pi / 180.

    ### Get cartesian components of the velocity fields instead of radial and azimuthal

    # get meshgrids for polar coordinates
    R = np.sqrt(X**2 + Y**2)
    PHI = np.arctan2(Y, X)

    # perform transformations
    v_x = -v_phi * np.sin(PHI) + v_r * np.cos(PHI)
    v_y = v_phi * np.cos(PHI) + v_r * np.sin(PHI)
    v_z = np.zeros(v_x.shape)

    # define velocity field
    v_field = np.array([v_x, v_y, v_z])

    plt.contourf(X, Y, 1e3*v_field[0,:,:], levels=np.linspace(-200,200,199), cmap="RdBu")
    plt.title("v_x")
    plt.colorbar()
    plt.show()
    
    plt.contourf(X, Y, 1e3*v_field[1,:,:], levels=np.linspace(-200,200,199), cmap="RdBu")
    plt.title("v_y")
    plt.colorbar()
    plt.show()

    ### Rotate velocity field and grid
    
    print("Rotating velocity fields")

    # grid shape from x grid
    N_x = X.shape[0]
    N_y = X.shape[1]
    
    # check y grid has same shape
    assert Y.shape[0] == N_x
    assert Y.shape[1] == N_y

    # loop over all points
    for i in range(N_x):
        for j in range(N_y):

            # rotate around the normal axis of the disk, corresponding the planet_az angle
            rot_pl_z = Rot(-planet_az, "z")
            X[i,j], Y[i,j], _   = np.dot(rot_pl_z, [X[i,j], Y[i,j], 0])
            v_field[:,i,j]      = np.dot(rot_pl_z, v_field[:,i,j])

            # rotate around the x-axis of the sky plane to match the inclination
            rot_in_x = Rot(inclination, "x")
            X[i,j], Y[i,j], _   = np.dot(rot_in_x, [X[i,j], Y[i,j], 0])
            v_field[:,i,j]      = np.dot(rot_in_x, v_field[:,i,j])

            # rotate around the normal axis of the sky plane to match the PA
            rot_pa_z = Rot(PA, "z")
            X[i,j], Y[i,j], _   = np.dot(rot_pa_z, [X[i,j], Y[i,j], 0])
            v_field[:,i,j]      = np.dot(rot_pa_z, v_field[:,i,j])
    
    ### plot z-axis velocities
    plt.contourf(X, Y, 1e3*v_field[2,:,:], levels=np.linspace(-200,200,199), cmap="RdBu")
    #plt.contourf(X, Y, 1e3*v_field[2,:,:], levels=np.linspace(-5000,5000,199), cmap="RdBu")
    plt.colorbar()
    plt.show()


def Rot(ang, ax='x'):
    """Function for rotation matrices"""
    
    if ax == "x":
        return [
            [1, 0, 0],
            [0, np.cos(ang), -np.sin(ang)],
            [0, np.sin(ang), np.cos(ang)]
        ]
    elif ax == "y":
        return [
            [np.cos(ang), 0, np.sin(ang)],
            [0, 1, 0],
            [-np.sin(ang), 0, np.cos(ang)]
        ]
    elif ax == "z":
        return [
            [np.cos(ang), -np.sin(ang), 0],
            [np.sin(ang), np.cos(ang), 0],
            [0, 0, 1]
        ]
    else:
        raise ValueError("ax must be x, y or z")
    
    
    
if __name__ == '__main__':
    main()