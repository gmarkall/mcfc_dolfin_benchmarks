from dolfin import *
from initial_condition import val
from parameters import *

#parameters["num_threads"] = 6
#parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

class InitialCondition(Expression):
    def __init__(self, fn):
        self._fn = fn

    def eval(self, values, x):
        values[0] = self._fn(x, 0.1)

def simulation(D, A, t, dt, endtime, mesh, initial):

    # Added due to mesh not conforming to UFC numbering, allegedly
    mesh.order()

    # Create FunctionSpaces
    T = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)

    # Create velocity Function
    velocity = Constant( (5.0, 0.0) )

    concentration = InitialCondition(initial)

    # Initialise source function and previous solution function
    u0 = Function(T)
    u1 = Function(T)
    u1.assign(concentration)

    # Test and trial functions
    u, v = TrialFunction(T), TestFunction(T)

    # Advection
    Mass = v*u*dx
    adv_rhs = (v*u0+dt*dot(grad(v),velocity)*u0)*dx

    #Diffusion
    d = -dt*D*dot(grad(v),grad(u))*dx
    diff_matrix = Mass - 0.5*d
    diff_rhs = action(Mass + 0.5*d, u0)

    out_file = File("temperature.pvd")
    out_file << (u1, t)

    t_a_a = Timer("advection assembly")
    t_a_s = Timer("advection solve")
    t_d_a = Timer("diffusion assembly")
    t_d_s = Timer("diffusion solve")

    t_loop = Timer("Timestepping loop")

    # Time-stepping
    while t < endtime:

        # Copy soln from prev.
        u0.assign(u1)
        print "Max:", u0.vector().max()

        # Advection
        t_a_a.start()
        M = assemble(Mass)
        b = assemble(adv_rhs)
        t_a_a.stop()
        t_a_s.start()
        solve(M, u1.vector(), b, "cg", "jacobi")
        t_a_s.stop()

        # Copy solution from advection
        u0.assign(u1)

        # Diffusion
        t_d_a.start()
        M = assemble(diff_matrix)
        b = assemble(diff_rhs)
        t_d_a.stop()
        t_d_s.start()
        solve(M, u1.vector(), b, "cg", "jacobi")
        t_d_s.stop()
        #if save_output:
        #    out_file << (u1, t)
     
        # Next timestep
        t += dt

    t_loop.stop()

    list_timings()
    
    out_file << (u1, t)

    return t_loop.value()


# Load mesh
mesh = Mesh("../mesh/cdisk.xml")

simulation(D, A, t, dt, endtime, mesh, val)
