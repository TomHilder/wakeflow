from setup import run_setup
from grid import Grid
from linear_perts import LinearPerts
import numpy as np

params = run_setup()

g1 = Grid(params)

g1.make_grid()
g1.make_keplerian_disk()
g1.show_disk2D(40)

print(g1.info)

lin = LinearPerts(params)
lin.cut_box_annulus_segment()