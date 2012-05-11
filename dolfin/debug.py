from dolfin import *
from math import sqrt, exp, pi

# Parameters

D = 0.1       # Diffusivity
A = 0.1       # Normalisation
t = 0.1       # Starting time
dt = 0.00015625  # Timestep
endtime = 0.2 # Simulation finish time

# Initial condition

def val(X, t):
    x = X[0]-0.5
    y = X[1]-0.5
    r = sqrt(x*x + y*y)
    d = 0.05
    if r<0.25:
        return A*(exp((-r**2)/(4*d*t))/(4*pi*d*t))
    else:
        return 0.0

class InitialCondition(Expression):
    def __init__(self, fn):
        self._fn = fn

    def eval(self, values, x):
        values[0] = self._fn(x, 0.1)

def simulation(D, A, t, dt, endtime, mesh, initial):

    # Added due to mesh not conforming to UFC numbering
    mesh.order()

    # Create FunctionSpaces
    T = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)


    # Initialise source function and previous solution function
    u0 = Function(T)
    u1 = Function(T)
    concentration = InitialCondition(initial)
    u1.assign(concentration)
    u, v = TrialFunction(T), TestFunction(T)

    Mass = v*u*dx
    d = -dt*D*dot(grad(v),grad(u))*dx
    diff_matrix = Mass - 0.5*d
    diff_rhs = action(Mass + 0.5*d, u0)

    out_file = File("temperature.pvd")
    out_file << (u1, t)

    fcp = { "optimize": True }

    # Time-stepping
    while t < endtime:
        u0.assign(u1)
        print "Max:", u0.vector().max()
        solve(diff_matrix==diff_rhs, u1, form_compiler_parameters=fcp)
        t += dt

    out_file << (u1, t)

mesh = UnitSquare(32, 32)

simulation(D, A, t, dt, endtime, mesh, val)
