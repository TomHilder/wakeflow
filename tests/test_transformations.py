from wakeflow.transformations import _mod2pi, _phi_wake, _Eta, _t, _g, _get_dens_vel, _Eta_vector, _t_vector
from pytest                   import approx
import numpy as np

def test_phi_wake():
    assert _phi_wake(1.5, 1, 0.1, 0.2, -1) == approx(-1.4048135831564323)
    assert _phi_wake(0.2, 1, 0.1, 0.2, -1) == approx( 13.563185995888238)
    assert _phi_wake(1.5, 1, 0.1, 0.2,  1) == approx( 1.4048135831564323)
    assert _phi_wake(0.2, 1, 0.1, 0.2,  1) == approx(-13.563185995888238)

def test_Eta():
    assert _Eta(1.5, 0.2, 1, 0.1, 0.2, -1) == approx( 24.07220374734649)
    assert _Eta(0.2, 0.2, 1, 0.1, 0.2, -1) == approx(-11.952230722935989)
    assert _Eta(1.5, 0.2, 1, 0.1, 0.2,  1) == approx(-18.072203747346485)
    assert _Eta(0.2, 0.2, 1, 0.1, 0.2,  1) == approx( 17.95223072293599)
    assert _Eta(1.5,   5, 1, 0.1, 0.2,  1) == approx(-40.31998335504028)

def test_Eta_vector():
    rad = np.array([0.2, 1.5])
    assert _Eta_vector(rad, 0.2, 1, 0.1, 0.2, -1)[0] == approx(-11.952230722935989)
    assert _Eta_vector(rad, 0.2, 1, 0.1, 0.2, -1)[1] == approx( 24.07220374734649)
    assert _Eta_vector(rad, 0.2, 1, 0.1, 0.2,  1)[0] == approx( 17.95223072293599)
    assert _Eta_vector(rad, 0.2, 1, 0.1, 0.2,  1)[1] == approx(-18.072203747346485)

def test_mod2pi():
    assert _mod2pi(      np.pi) == approx(np.pi)
    assert _mod2pi(     -np.pi) == approx(np.pi)
    assert _mod2pi(  3/4*np.pi) == approx((3/4) * np.pi)
    assert _mod2pi( -3/4*np.pi) == approx((5/4) * np.pi)
    assert _mod2pi( 45.2*np.pi) == approx(3.769911184307759)
    assert _mod2pi(-45.2*np.pi) == approx(2.5132741228718274)

def test_t():
    assert _t(1.5, 1, 0.1, 0.2,   1) == approx(34.671246598896204)
    assert _t(0.2, 1, 0.1, 0.2,   1) == approx(761.9720999691256)
    assert _t(1.5, 1, 0.1, 0.2, 0.5) == approx(32.23043221899377)
    assert _t(0.2, 1, 0.1, 0.2, 0.5) == approx(1016.399602991283)

def test_t_vector():
    rad = np.array([0.2, 1.5])
    assert _t_vector(rad, 1, 0.1, 0.2,   1)[0] == approx(761.9720999691256)
    assert _t_vector(rad, 1, 0.1, 0.2,   1)[1] == approx(34.671246598896204)
    assert _t_vector(rad, 1, 0.1, 0.2, 0.5)[0] == approx(1016.399602991283)
    assert _t_vector(rad, 1, 0.1, 0.2, 0.5)[1] == approx(32.23043221899377)

def test_g(): # r, Rp, hr, q, p
    assert _g(1.5, 1, 0.1, 0.2,   1) == approx(0.49329351377521874)
    assert _g(0.2, 1, 0.1, 0.2,   1) == approx(0.19101495094609305)
    assert _g(1.5, 1, 0.1, 0.2, 0.5) == approx(0.5459190128004762)
    assert _g(0.2, 1, 0.1, 0.2, 0.5) == approx(0.12773939655074654)

def test_get_dens_vel():
    assert _get_dens_vel(1.5, 1.0, 5/3, 1, -1, 0.1, 0.1, 0.2, 1, False)   == approx((1.5203929892776895, 0.09983362447771839, 0.013468458590558391))
    assert _get_dens_vel(0.2, 1.0, 5/3, 1, -1, 0.1, 0.1, 0.2,   1, False) == approx((3.9263942235163563, -0.29038210913525486, -0.019677575285259364))
    assert _get_dens_vel(1.5, 1.0, 5/3, 1, -1, 0.1, 0.1, 0.2, 0.5, False) == approx((1.3738301513856814, 0.09239015780274291, 0.01246426763579161))
    assert _get_dens_vel(0.2, 1.0, 5/3, 1, -1, 0.1, 0.1, 0.2, 0.5, False) == approx((5.871328816729225, -0.37299792971384055, -0.025275988472730864))
