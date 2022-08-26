from wakeflow               import WakeflowModel
from wakeflow.grid          import _Grid
from wakeflow.model_setup   import _Parameters
from copy                   import copy
import numpy                as np
import pytest

def cart_grid(dimless=False):

    # create dictionary of parameters
    model = WakeflowModel()
    model.configure(grid_type="cartesian", dimensionless=dimless)
    param_dict = model.model_params

    # intialise grid with parameters object
    grid = _Grid(_Parameters(param_dict))

    # check initial properties
    assert grid.info == {
        'Type': None, 
        'Log_r': None, 
        'Size': [0, 0, 0], 
        'Contains': 'Empty'
    }

    # setup grid geometry
    grid._make_grid()

    return grid

def cyli_grid(dimless=False):

    # create dictionary of parameters
    model = WakeflowModel()
    model.configure(grid_type="cylindrical", dimensionless=dimless)
    param_dict = model.model_params

    # intialise grid with parameters object
    grid = _Grid(_Parameters(param_dict))

    # check initial properties
    assert grid.info == {
        'Type': None, 
        'Log_r': None, 
        'Size': [0, 0, 0], 
        'Contains': 'Empty'
    }

    # setup grid geometry
    grid._make_grid()

    return grid

def test_cartesian():

    grid = cart_grid()

    # check grid geometry
    assert grid.x.min() == grid.y.min() == -500
    assert grid.x.max() == grid.y.max() ==  500

    # check grid properties
    assert grid.info == {
        'Type': 'cartesian', 
        'Log_r': None, 
        'Size': [400, 50, 400], 
        'Contains': 'Empty'
    }

    grid._get_r_phi_coords()

    assert grid.R_xy.max() == np.sqrt(2 * 500**2)
    assert grid.PHI_xy.max() == -grid.PHI_xy.min()
    assert (np.pi - grid.PHI_xy.max()) / np.pi < 1e-3

def test_cylindrical():

    grid = cyli_grid()

    # check grid geometry
    assert grid.r.min() == 100
    assert grid.r.max() == 500
    assert -grid.phi.min() == grid.phi.max() == np.pi

    # check grid properties
    assert grid.info == {
        'Type': 'cylindrical', 
        'Log_r': False, 
        'Size': [160, 50, 200], 
        'Contains': 'Empty'
    }

def test_dimless():

    ca_g = cart_grid(dimless=True)
    cy_g = cyli_grid(dimless=True)

    # fill fields with zeros
    ca_g._make_empty_disk()
    cy_g._make_empty_disk()

    # fill cartesian fields with ones
    ca_g.v_r   = np.ones(ca_g.v_r  .shape)
    ca_g.v_phi = np.ones(ca_g.v_phi.shape)
    ca_g.rho   = np.ones(ca_g.rho  .shape)

    # fill cylindrical fields with ones
    cy_g.v_r   = np.ones(cy_g.v_r  .shape)
    cy_g.v_phi = np.ones(cy_g.v_phi.shape)
    cy_g.rho   = np.ones(cy_g.rho  .shape)

    # get velocity unit
    assert ca_g.p.U_vel == cy_g.p.U_vel
    U_vel = ca_g.p.U_vel

    # create arrays full of 1 / unit with right shapes
    ca_u_arr = np.full(ca_g.v_r.shape, 1/U_vel)
    cy_u_arr = np.full(cy_g.v_r.shape, 1/U_vel)

    # remove dimensions on both
    ca_g._remove_dimensions()
    cy_g._remove_dimensions()

    # check velocities
    assert np.array_equal(ca_g.v_r,   ca_u_arr)
    assert np.array_equal(ca_g.v_phi, ca_u_arr)
    assert np.array_equal(cy_g.v_r,   cy_u_arr)
    assert np.array_equal(cy_g.v_phi, cy_u_arr)

    # check densities are all still ones
    assert np.array_equal(ca_g.rho, np.ones(ca_g.rho.shape))
    assert np.array_equal(cy_g.rho, np.ones(cy_g.rho.shape))

    # check lengths cartesian
    assert ca_g.x.min() == ca_g.y.min() == -2
    assert ca_g.x.max() == ca_g.y.max() ==  2

    # check lengths cylindrical
    assert cy_g.r.min() == 0.4
    assert cy_g.r.max() == 2
    assert -cy_g.phi.min() == cy_g.phi.max() == np.pi

def test_merge():

    ca_g   = cart_grid()
    ca_g_2 = copy(ca_g)
    cy_g   = cyli_grid()

    with pytest.raises(Exception):
        ca_g._merge_grids(cy_g) # give grid with wrong parameters
        ca_g._merge_grids(None) # give wrong object

    # fill one grid with background and one with zeros
    ca_g  ._make_keplerian_disk()
    ca_g_2._make_empty_disk()

    # add the background to the zeros
    ca_g_2._merge_grids(ca_g)

    assert np.array_equal(ca_g_2.v_r,   ca_g.v_r)
    assert np.array_equal(ca_g_2.v_phi, ca_g.v_phi)
    assert np.array_equal(ca_g_2.rho,   ca_g.rho)
