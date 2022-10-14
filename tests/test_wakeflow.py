from wakeflow import WakeflowModel
import pkg_resources

# test that you can build the model object
def test_initialisation():
    model = WakeflowModel()
    assert type(model) == WakeflowModel

# test that the default configuration is correct
def test_default_user_config():

    # configure default model
    model = WakeflowModel()
    model.configure()

    # check default params
    assert model.model_params["name"]                == "default_results_dir"
    assert model.model_params["system"]              == "default_parent_dir"
    assert model.model_params["m_star"]              == 1.0
    assert model.model_params["m_planet"]            == 0.1
    assert model.model_params["r_outer"]             == 500
    assert model.model_params["r_inner"]             == 100
    assert model.model_params["r_planet"]            == 250
    assert model.model_params["r_ref"]               == None
    assert model.model_params["r_c"]                 == 0
    assert model.model_params["q"]                   == 0.25
    assert model.model_params["p"]                   == 1.0
    assert model.model_params["hr"]                  == 0.10
    assert model.model_params["dens_ref"]            == 1.0
    assert model.model_params["cw_rotation"]         == False
    assert model.model_params["grid_type"]           == "cartesian"
    assert model.model_params["n_x"]                 == 400
    assert model.model_params["n_y"]                 == 400
    assert model.model_params["n_r"]                 == 200
    assert model.model_params["n_phi"]               == 160
    assert model.model_params["n_z"]                 == 50
    assert model.model_params["r_log"]               == False
    assert model.model_params["make_midplane_plots"] == True
    assert model.model_params["show_midplane_plots"] == True
    assert model.model_params["dimensionless"]       == False
    assert model.model_params["include_planet"]      == True
    assert model.model_params["include_linear"]      == True
    assert model.model_params["save_perturbations"]  == True
    assert model.model_params["save_total"]          == True
    assert model.model_params["write_FITS"]          == False
    assert model.model_params["run_mcfost"]          == False
    assert model.model_params["inclination"]         == -225
    assert model.model_params["PA"]                  == 45
    assert model.model_params["PAp"]                 == 45
    assert model.model_params["temp_star"]           == 9250
    assert model.model_params["distance"]            == 101.5
    assert model.model_params["v_max"]               == 3.2
    assert model.model_params["n_v"]                 == 40

# check that the default developer and numerical values are correct
def test_default_dev_config():

    # configure default model
    model = WakeflowModel()
    model.configure()

    # check default params
    assert model.model_params["adiabatic_index"]       == 1.6666667
    assert model.model_params["damping_malpha"]        == 0.0
    assert model.model_params["CFL"]                   == 0.5
    assert model.model_params["scale_box_l"]           == 1.0
    assert model.model_params["scale_box_r"]           == 1.0
    assert model.model_params["scale_box_ang_b"]       == 1.0
    assert model.model_params["scale_box_ang_t"]       == 1.0
    assert model.model_params["tf_fac"]                == 1.0
    assert model.model_params["show_teta_debug_plots"] == False
    assert model.model_params["box_warp"]              == True
    assert model.model_params["use_box_IC"]            == False

# test that user configuration overrides default
def test_manual_config():

    # configure model
    model = WakeflowModel()
    model.configure(m_planet=0.1, r_planet=100)

    # check it worked
    assert model.model_params["m_planet"] == 0.1
    assert model.model_params["r_planet"] == 100

# test that wakeflow can be configured from .yaml file using the default file
def test_file_config():

    # find default config file
    file_loc = pkg_resources.resource_filename('wakeflow', 'data/default_config.yaml')

    # configure model with .yaml file containing defaults
    model = WakeflowModel()
    model.configure_from_file(file_loc)

    # check default params
    assert model.model_params["name"]                == "default_results_dir"
    assert model.model_params["system"]              == "default_parent_dir"
    assert model.model_params["m_star"]              == 1.0
    assert model.model_params["m_planet"]            == 0.1
    assert model.model_params["r_outer"]             == 500
    assert model.model_params["r_inner"]             == 100
    assert model.model_params["r_planet"]            == 250
    assert model.model_params["r_ref"]               == None
    assert model.model_params["r_c"]                 == 0
    assert model.model_params["q"]                   == 0.25
    assert model.model_params["p"]                   == 1.0
    assert model.model_params["hr"]                  == 0.10
    assert model.model_params["dens_ref"]            == 1.0
    assert model.model_params["cw_rotation"]         == False
    assert model.model_params["grid_type"]           == "cartesian"
    assert model.model_params["n_x"]                 == 400
    assert model.model_params["n_y"]                 == 400
    assert model.model_params["n_r"]                 == 200
    assert model.model_params["n_phi"]               == 160
    assert model.model_params["n_z"]                 == 50
    assert model.model_params["r_log"]               == False
    assert model.model_params["make_midplane_plots"] == True
    assert model.model_params["show_midplane_plots"] == True
    assert model.model_params["dimensionless"]       == False
    assert model.model_params["include_planet"]      == True
    assert model.model_params["include_linear"]      == True
    assert model.model_params["save_perturbations"]  == True
    assert model.model_params["save_total"]          == True
    assert model.model_params["write_FITS"]          == False
    assert model.model_params["run_mcfost"]          == False
    assert model.model_params["inclination"]         == -225
    assert model.model_params["PA"]                  == 45
    assert model.model_params["PAp"]                 == 45
    assert model.model_params["temp_star"]           == 9250
    assert model.model_params["distance"]            == 101.5
    assert model.model_params["v_max"]               == 3.2
    assert model.model_params["n_v"]                 == 40

# check that running the default model works
def test_run_default():
    model = WakeflowModel()
    # turn saving off and don't show plots
    model.configure(
        n_x=100,
        n_y=100,
        show_midplane_plots=False, 
        make_midplane_plots=False, 
        save_perturbations=False, 
        save_total=False
    )
    model.run()

def test_run_cylindrical():
    model = WakeflowModel()
    # turn saving off and don't show plots
    model.configure(
        grid_type="cylindrical",
        n_r=100,
        n_phi=100,
        show_midplane_plots=False, 
        make_midplane_plots=False, 
        save_perturbations=False, 
        save_total=False
    )
    model.run()

def test_run_cylindrica_dev_option():
    model = WakeflowModel()
    # turn saving off and don't show plots
    model.configure(
        grid_type="cylindrical",
        n_r=100,
        n_phi=100,
        r_log=True,
        show_midplane_plots=False, 
        make_midplane_plots=False, 
        save_perturbations=False, 
        save_total=False
    )
    model.model_params["box_warp"] = False
    model.run()
