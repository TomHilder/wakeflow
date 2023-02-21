from astropy                import constants        as c
from wakeflow.model_setup   import _Constants, _Parameters, _load_config_file, _write_config_file, _run_setup
from wakeflow               import WakeflowModel
import pytest, pkg_resources, os

def test_constants():
    const = _Constants()
    assert const.m_solar   == c.M_sun.cgs.value
    assert const.m_jupiter == c.M_jup.cgs.value
    assert const.G_const   == c.G.cgs.value
    assert const.au        == c.au.cgs.value

def test_params():

    # make parameter object
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params
    params = _Parameters(param_dict)

    # check it worked
    assert type(params) == _Parameters

def test_multiple_planet_mass():

     # make parameter object
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params
    param_dict["m_planet"] = [0.1, 0.2]
    params = _Parameters(param_dict)

    # check it worked
    assert type(params) == _Parameters

def test_sanity_good():

    # make parameter object
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params
    params = _Parameters(param_dict)

    # check that sanity checks returns True for default parameters
    assert params._do_sanity_checks()

def test_sanity_invalid_grid():

    # make parameter dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # change grid to invalid option
    param_dict["grid_type"] = "bogus"
    params = _Parameters(param_dict)

    # check that sanity checks fail
    with pytest.raises(Exception):
        params._do_sanity_checks()

def test_sanity_invalid_mcfost_grid():

    # make parameter dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # change to invalid option of cartesian grid with mcfost
    param_dict["grid_type"] = "cartesian"
    param_dict["run_mcfost"] = True
    params = _Parameters(param_dict)

    # check that sanity checks fail
    with pytest.raises(Exception):
        params._do_sanity_checks()

def test_sanity_invalid_mcfost_no_fits():

    # make parameter dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # try to run mcfost without writing fits file
    param_dict["grid_type"]  = "mcfost"
    param_dict["write_fits"] = False
    param_dict["run_mcfost"] = True
    params = _Parameters(param_dict)
    
    # check that sanity checks fail
    with pytest.raises(Exception):
        params._do_sanity_checks()

def test_m_thermal_warning(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["m_planet"] = 100
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: M_planet > M_thermal" in captured.out

def test_m_thermal_warning_multiple_masses(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["m_planet"] = [0.01, 100]
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: At least one planet mass exceeds thermal mass. This may break the solution. " in captured.out

def test_box_scale_l_warning(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["scale_box_l"] = 1.01
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: Changing linear box scale factor can cause strange results." in captured.out
    
def test_box_scale_r_warning(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["scale_box_r"] = 1.01
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: Changing linear box scale factor can cause strange results." in captured.out
    
def test_box_scale_lr_warning(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["scale_box_l"] = 1.01
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: Using a different linear box scale factor for left and right edge can cause strange results." in captured.out
    
def test_box_scale_tb_warning(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["scale_box_ang_t"] = 1.01
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: Using a different linear box scale factor for top and bottom edge can cause strange results." in captured.out


def test_CFL_warning(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["CFL"] = 0.9
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: CFL chosen > 0.5, this will likely break the numerical PDE solver." in captured.out

def test_box_warp_warning(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["box_warp"] = False
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: Choosing box_warp=False may invalidate your results." in captured.out

def test_box_IC_warning(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["use_box_IC"] = True
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: Choosing use_box_IC=True will almost certainly invalidate your results." in captured.out

def test_show_teta_warning(capfd):

    # make params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # make planet huge
    param_dict["show_teta_debug_plots"] = True
    params = _Parameters(param_dict)

    # check that sanity checks warns about planet mass
    params._do_sanity_checks()
    captured = capfd.readouterr()
    assert "WARNING: Choosing show_teta_debug_plots=True may cause the run to fail." in captured.out

def test_load_config_file():

    # make default params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # read in default params file
    file_loc = pkg_resources.resource_filename('wakeflow', 'data/default_config.yaml')
    config_dict = _load_config_file(file_loc)

    assert param_dict == config_dict

def test_load_config_file_bad():

    # this is cheating a bit but instead of changing the parameter file we just pretend that the default parameter dictionary does not contain m_star

    # make default params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # delete an etry
    del param_dict["m_star"]

    # read in default params file
    file_loc = pkg_resources.resource_filename('wakeflow', 'data/default_config.yaml')
    with pytest.raises(Exception):
        config_dict = _load_config_file(file_loc, param_dict)

def test_write_config_yaml(tmp_path):

    # make default params dict
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # write file
    _write_config_file(param_dict, tmp_path, "/test_config.yaml")

    # read file again
    config_dict_read = _load_config_file(f"{tmp_path}/test_config.yaml")

    assert config_dict_read == param_dict

def test_run_setup(tmp_path):

    # make parameter object
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # change directory to temp
    os.chdir(tmp_path)

    # run setup
    params = _run_setup(param_dict)

def test_run_setup_multiple_masses(tmp_path):

    # make parameter object
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params
    param_dict["m_planet"] = [0.1, 0.2]

    # change directory to temp
    os.chdir(tmp_path)

    # run setup
    params = _run_setup(param_dict)

def test_run_setup_overwrite(tmp_path):

    # make parameter object
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # change directory to temp
    os.chdir(tmp_path)

    # run setup
    params = _run_setup(param_dict)

    # run setup with overwrite
    params2 = _run_setup(param_dict, overwrite=True)

def test_run_setup_bad_overwrite(tmp_path):

    # make parameter object
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # change directory to temp
    os.chdir(tmp_path)

    # run setup
    params = _run_setup(param_dict)

    # run setup with overwrite
    with pytest.raises(Exception):
        params2 = _run_setup(param_dict)
