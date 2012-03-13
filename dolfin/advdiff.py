from dolfin import *

def initial_condition(mesh, t, val):
    nv = mesh.num_vertices()
    coordinates = mesh.coordinates()
    concentration = Vector(nv)
    for i in range(nv):
        X = coordinates[i]
        concentration[i] = val(X, t)
    return concentration

class DolfinSimulation(object):

    def __init__(self, D, A, t, dt, endtime, mesh, initial):
        self.D = D
        self.A = A
        self.t = t
        self.dt = dt
        self.endtime = endtime
        self.mesh = Mesh(mesh)
        self.initial = initial

    def run(self, order, save_output=True):
        # Grab parameters
        D = self.D
        A = self.A
        t = self.t
        dt = self.dt
        endtime = self.endtime
        mesh = self.mesh
        initial = self.initial

        # Create FunctionSpaces
        T = FunctionSpace(mesh, "CG", order)
        V = VectorFunctionSpace(mesh, "CG", order)

        # Create velocity Function
        velocity = Constant( (5.0, 0.0) )

        concentration = initial_condition(mesh, t, initial)

        # Initialise source function and previous solution function
        u0 = Function(T)
        u1 = Function(T, concentration)

        # Test and trial functions
        u, v = TrialFunction(T), TestFunction(T)

        # Advection
        Mass = v*u*dx
        adv_rhs = (v*u0+dt*dot(grad(v),velocity)*u0)*dx

        #Diffusion
        d = -dt*D*dot(grad(v),grad(u))*dx
        diff_matrix = Mass - 0.5*d
        diff_rhs = action(Mass + 0.5*d, u0)

        if save_output:
            out_file = File("temperature.pvd")
            out_file << (u1, t)

        t_a = Timer("advection")
        t_d = Timer("diffusion")

        t_loop = Timer("Timestepping loop")

        # Time-stepping
        while t < endtime:

            # Copy soln from prev.
            u0.assign(u1)
            
            # Advection
            t_a.start()
            M = assemble(Mass)
            b = assemble(adv_rhs)
            solve(M, u1.vector(), b)
            t_a.stop()

            # Copy solution from advection
            u0.assign(u1)

            # Diffusion
            t_d.start()
            M = assemble(diff_matrix)
            b = assemble(diff_rhs)
            solve(M, u1.vector(), b)
            t_d.stop()
            if save_output:
                out_file << (u1, t)
         
            # Next timestep
            t += dt

        t_loop.stop()
        if save_output:
            out_file << (u1, t)
