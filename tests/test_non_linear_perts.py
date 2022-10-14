from wakeflow                   import WakeflowModel
from wakeflow.grid              import _Grid
from wakeflow.linear_perts      import _LinearPerts
from wakeflow.non_linear_perts  import _NonLinearPerts
from wakeflow.model_setup       import _Parameters
from pytest                     import approx
import numpy                        as np
import pytest

def default_Params():

    # create dictionary of parameters
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # intialise grid with parameters object
    return _Parameters(param_dict)

def default_LinearPerts():

    # default params
    def_params = default_Params()

    # get linear perts and do cut
    lp = _LinearPerts(def_params)
    lp._cut_box_annulus_segment()

    return lp

def default_NonLinearPerts():

    # parameters and linear perts
    def_params = default_Params()

    # setup grid
    nlp_grid = _Grid(def_params)
    nlp_grid._make_grid()
    nlp_grid._make_empty_disk()

    # nonlinear perts object
    return _NonLinearPerts(def_params, nlp_grid)

def test_extract_ICs_default():

    # default nonlinear perts object
    nlp = default_NonLinearPerts()

    # default linear perts object
    lp = default_LinearPerts()

    # extract initial conditions from linear
    nlp._extract_ICs(lp)

    # check ICs
    assert nlp.eta_tilde_inner == approx(1.707850734468749)
    assert nlp.eta_tilde_outer == approx(3.1189508931861276)
    assert nlp.C_inner         == approx(0.1406077053627964)
    assert nlp.C_outer         == approx(0.08439583277062648)
    assert nlp.t0_inner        == approx(2.3034788704426026)
    assert nlp.t0_outer        == approx(1.608427450463994)

def test_alternate_extract_ICs_default():
    # this should not work

    # default nonlinear perts object
    nlp = default_NonLinearPerts()

    # default linear perts object
    lp = default_LinearPerts()

    # extract initial conditions from linear using alternate method
    with pytest.raises(UnboundLocalError):
        nlp._extract_ICs_ann(lp)

def test_alternate_extract_ICs_works():

    # create dictionary of parameters
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # extend annulus
    param_dict["scale_box_ang_b"] = 4.0
    param_dict["scale_box_ang_t"] = 4.0

    # intialise grid with parameters object
    params = _Parameters(param_dict)

    # get linear perts and do cut
    lp = _LinearPerts(params)
    lp._cut_box_annulus_segment()

    # setup grid
    nlp_grid = _Grid(params)
    nlp_grid._make_grid()
    nlp_grid._make_empty_disk()

    # nonlinear perts object
    nlp = _NonLinearPerts(params, nlp_grid)

    # extract initial conditions from linear using alternate method
    nlp._extract_ICs_ann(lp)

    # check ICs
    assert nlp.eta_tilde_inner == approx(1.6859265632683058)
    assert nlp.eta_tilde_outer == approx(2.627559133795374)
    assert nlp.C_inner         == approx(0.19149042489930776)
    assert nlp.C_outer         == approx(0.02144357556464586)
    assert nlp.t0_inner        == approx(2.3034788704426026)
    assert nlp.t0_outer        == approx(1.608427450463994)
