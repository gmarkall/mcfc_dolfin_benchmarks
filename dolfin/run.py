from advdiff import DolfinSimulation
from initial_condition import val
from parameters import *

# Load mesh
mesh = "../mesh/cdisk.xml"

sim = DolfinSimulation(D, A, t, dt, endtime, mesh, val)

sim.run(1, False)
