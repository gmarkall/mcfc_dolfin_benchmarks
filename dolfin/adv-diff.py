from dolfin import *

# Load mesh and subdomains
mesh = Mesh("../mesh/cdisk.xml")

nv = mesh.num_vertices()
print "Num vertices ", nv

coordinates = mesh.coordinates()
concentration = Vector(nv)

# Parameters
D = 0.05     # Diffusivity
A = 0.1      # Normalisation
t = 0.1      # Starting time
dt = 0.00125 # Timestep

endtime = 0.7

for i in range(mesh.num_vertices()):
    x = coordinates[i][0]+0.5
    y = coordinates[i][1]
    r = sqrt(x*x + y*y)
    if r<0.25:
        concentration[i] = A*(exp((-r**2)/(4*D*t))/(4*pi*D*t))
    else:
        concentration[i] = 0.0

# Create FunctionSpaces
T = FunctionSpace(mesh, "CG", 1)
V = VectorFunctionSpace(mesh, "CG", 1)

# Create velocity Function
velocity = Constant( (5.0, 0.0) )

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


plot(u1)
interactive()

# Set intial condition

# Time-stepping
while t < endtime:

    # Copy soln from prev.
    u0.assign(u1)
    
    # Advection
    M = assemble(Mass)
    b = assemble(adv_rhs)
    solve(M, u1.vector(), b)

    # Diffusion
    M = assemble(diff_matrix)
    b = assemble(diff_rhs)
    solve(M, u1.vector(), b)

    # Plot solution
    plot(u1)

    # Save the solution to file
    #out_file << (u1, t)

    # Move to next interval and adjust boundary condition
    t += dt

plot(u1)
# Hold plot
interactive()
