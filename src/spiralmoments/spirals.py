# spirals.py
# Written by Thomas Hilder

class Spiral():
    
    # coordinates of spiral
    _rad = None
    _phi = None
    _z   = None
    
    @property
    def rad_phi(self):
        return self._rad, self._phi
    
    @property
    def height(self):
        return self._z