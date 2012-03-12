from mpi4py import *
from dolfin import *
from math import exp, pi, sqrt

mesh = Mesh("cdisk.xml")

nv = mesh.num_vertices()
print "Num vertices", nv

coordinates = mesh.coordinates()
concentration = Vector(nv)

diffusivity = 0.1
t=0.1
dt = 0.1
end = 2.0

for i in range(nv):
    x = coordinates[i][0]
    y = coordinates[i][1]
    r = sqrt(x*x + y*y)
    D = 0.1
    A = 0.1
    concentration[i] = A*(exp((-r**2)/(4*D*t))/(4*pi*D*t))

T = FunctionSpace(mesh, "CG", 1)

u0 = Function(T)
u1 = Function(T, concentration)

q, p = TestFunction(T), TrialFunction(T)

M=q*p*dx
d=-dt*diffusivity*dot(grad(q),grad(p))*dx

A=M-0.5*d
rhs=action(M+0.5*d,u1)

out_file = File("temperature.pvd")

plot(u1)
interactive()


while t<end:

    u0.assign(u1)
    a = assemble(A)
    L = assemble(rhs)
    solve(a, u1.vector(), L)
    plot(u1)

    out_file << (u1, t)
    t += dt

interactive()

