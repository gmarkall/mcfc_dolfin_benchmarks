from advdiff import DolfinSimulation
from initial_condition import val
from parameters import *


# Number of element orders to benchmark for
num_orders = 1

# Load mesh
mesh = "../mesh/cdisk.xml"

sim = DolfinSimulation(D, A, t, dt, endtime, mesh, val)

for i in range(num_orders):
    order = i+1
    sim.run(order)

