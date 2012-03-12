from dolfin import *

# Load mesh and subdomains
mesh = Mesh("cdisk.xml")

nv = mesh.num_vertices()
print "Num vertices ", nv

coordinates = mesh.coordinates()
concentration = Vector(nv)

for i in range(mesh.num_vertices()):
    x = coordinates[i][0]+0.5
    y = coordinates[i][1]
    r = sqrt(x*x + y*y)
    if r<0.25:
        concentration[i] = 1.0
    else:
        concentration[i] = 0.0

# Create FunctionSpaces
T = FunctionSpace(mesh, "CG", 1)
V = VectorFunctionSpace(mesh, "CG", 1)

# Create velocity Function
velocity = Constant( (1.0, 0.0) )

# Initialise source function and previous solution function
u0 = Function(T)
u1 = Function(T, concentration)

# Test and trial functions
u, v = TrialFunction(T), TestFunction(T)

Mass = v*u*dx
a= assemble(Mass)

# params
endtime = 1.0
dt = 0.01
t=0.0

rhs = (v*u0+dt*dot(grad(v),velocity)*u0)*dx

out_file = File("temperature.pvd")

plot(u1)
interactive()

# Set intial condition

# Time-stepping
while t < endtime:

    # Copy soln from prev.
    u0.assign(u1)
    
    # Assemble vector and apply boundary conditions
    b = assemble(rhs)

    # Solve the linear system (re-use the already factorized matrix A)
    solve(a, u1.vector(), b)

    # Plot solution
    plot(u1)

    # Save the solution to file
    out_file << (u1, t)

    # Move to next interval and adjust boundary condition
    t += dt

# Hold plot
interactive()
