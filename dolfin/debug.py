from dolfin import *
from math import sqrt, exp, pi

# Parameters

D = 0.1       # Diffusivity
A = 0.1       # Normalisation
t = 0.1       # Starting time
dt = 0.00015625  # Timestep
endtime = 0.2 # Simulation finish time
mesh = UnitSquare(32, 32)

class InitialCondition(Expression):

    def eval(self, values, X):
        x = X[0]-0.5
        y = X[1]-0.5
        r = sqrt(x*x + y*y)
        t = 0.1
        d = 0.05
        if r<0.25:
            values[0] = A*(exp((-r**2)/(4*d*t))/(4*pi*d*t))
        else:
            values[0] = 0.0

# Create FunctionSpaces
T = FunctionSpace(mesh, "CG", 1)
V = VectorFunctionSpace(mesh, "CG", 1)

# Initialise source function and previous solution function
u0 = Function(T)
u1 = Function(T)
u, v = TrialFunction(T), TestFunction(T)

concentration = InitialCondition()
u1.assign(concentration)

Mass = v*u*dx
d = -dt*D*dot(grad(v),grad(u))*dx
diff_matrix = Mass - 0.5*d
diff_rhs = action(Mass + 0.5*d, u0)

out_file = File("temperature.pvd")
out_file << (u1, t)

# Time-stepping
while t < endtime:
    u0.assign(u1)
    print "Max:", u0.vector().max()
    solve(diff_matrix==diff_rhs, u1, form_compiler_parameters={ "optimize": True })
    t += dt

out_file << (u1, t)

