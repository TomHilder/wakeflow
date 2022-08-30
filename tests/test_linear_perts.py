from wakeflow               import WakeflowModel
from wakeflow.linear_perts  import _LinearPerts
from wakeflow.model_setup   import _Parameters
from pytest                 import approx
import numpy                    as np

def default_Params():
    # create dictionary of parameters
    model = WakeflowModel()
    model.configure()
    param_dict = model.model_params

    # intialise grid with parameters object
    return _Parameters(param_dict)

def test_read_linear_perts():
    lp = _LinearPerts(default_Params())

def test_cut_square():
    # we test this even though it is deprecated because it is still occasionally useful

    # get linear perts and do cut
    lp = _LinearPerts(default_Params())
    lp._cut_box_square()

    # check edges of cut
    assert lp.x_cut.min() == approx(-1.987 )
    assert lp.x_cut.max() == approx( 2.001 )
    assert lp.y_cut.min() == approx(-11.928)
    assert lp.y_cut.max() == approx( 12.026)

def test_cut_annulus():

    # get linear perts and do cut
    lp = _LinearPerts(default_Params())
    lp._cut_box_annulus_segment()

    assert np.abs(lp.R_ann.min() - 216.666) < 1e-3
    assert np.abs(lp.R_ann.max() - 283.333) < 1e-3
    assert np.abs(lp.PHI_ann.min() + 0.154) < 1e-3
    assert np.abs(lp.PHI_ann.max() - 0.154) < 1e-3