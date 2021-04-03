from setup import run_setup
from grid import Grid

p = run_setup()

g = Grid(p)

g.make_grid()
g.make_keplerian_disk()

#print(g.info)