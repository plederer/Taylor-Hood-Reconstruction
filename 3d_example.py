from ngsolve import *
from netgen.csg import unit_cube
from ngsolve.internal import visoptions
from ngsolve.internal import viewoptions

from reclib import *

mesh = Mesh(unit_cube.GenerateMesh (maxh=0.2))
nu = 1e-6

def POW(x,y):
        erg = x
        for i in range(1,y):
            erg *= x
        return erg

#right hand side for: zeta = x^2(1-x)^2y^2(1-y)^2z^2(1-z)^2, u = curl(zeta); p = x^5+y^5+z^4 - 1/2

f1 = -nu * (4 * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * POW(x,2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) + 4 * POW(x,2) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) - 4 * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) - 4 * POW(x,2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) - 4 * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) - 4 * POW(x,2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) + 12 * POW(x,2) * POW((x-1),2) * (y - 1) * POW(z,2) * POW((z-1),2) + 12 * POW(x,2) * POW((x-1),2) * y * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * POW((y-1),2) * z * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * z * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * POW((y-1),2) * POW(z,2) * (z - 1) - 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW(z,2) * (z - 1) + 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW((z-1),2) + 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) + 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * POW((z-1),2) + 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) - 12 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * (z - 1) - 12 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z + 16 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * z * (z - 1) + 16 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * z * (z - 1) + 16 * x * (x - 1) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 16 * x * (x - 1) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) - 16 * x * (x - 1) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) - 16 * x * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) - 16 * POW(x,2) * POW((x-1),2) * y * (y - 1) * z * POW((z-1),2) - 16 * POW(x,2) * POW((x-1),2) * y * (y - 1) * POW(z,2) * (z - 1)) + 5 * POW(x,4)


f2 =-nu * (4 * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) - 12 * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 12 * x * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * x * POW((x-1),2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * x * POW((x-1),2) * POW(y,2) * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * (x - 1) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * (x - 1) * POW(y,2) * POW(z,2) * POW((z-1),2) - 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW((z-1),2) - 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) - 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW((z-1),2) - 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) + 4 * POW(x,2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) + 4 * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) + 4 * POW(x,2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) + 4 * POW(x,2) * POW((x-1),2) * POW((y-1),2) * z * POW((z-1),2) + 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * z * POW((z-1),2) + 4 * POW(x,2) * POW((x-1),2) * POW((y-1),2) * POW(z,2) * (z - 1) + 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW(z,2) * (z - 1) + 12 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * (z - 1) + 12 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z - 16 * x * POW((x-1),2) * y * (y - 1) * POW(z,2) * POW((z-1),2) - 16 * POW(x,2) * (x - 1) * y * (y - 1) * POW(z,2) * POW((z-1),2) - 16 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * (z - 1) - 16 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * z * (z - 1) + 16 * x * (x - 1) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) + 16 * x * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) + 16 * POW(x,2) * POW((x-1),2) * y * (y - 1) * z * POW((z-1),2) + 16 * POW(x,2) * POW((x-1),2) * y * (y - 1) * POW(z,2) * (z - 1)) + 5 * POW(y,4)


f3 = -nu * (-4 * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) + 12 * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 12 * x * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * x * POW((x-1),2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * x * POW((x-1),2) * POW(y,2) * POW(z,2) * POW((z-1),2) + 4 * POW(x,2) * (x - 1) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * POW(x,2) * (x - 1) * POW(y,2) * POW(z,2) * POW((z-1),2) + 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW((z-1),2) + 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) + 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW((z-1),2) + 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) - 12 * POW(x,2) * POW((x-1),2) * (y - 1) * POW(z,2) * POW((z-1),2) - 12 * POW(x,2) * POW((x-1),2) * y * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) - 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) + 16 * x * POW((x-1),2) * y * (y - 1) * POW(z,2) * POW((z-1),2) + 16 * POW(x,2) * (x - 1) * y * (y - 1) * POW(z,2) * POW((z-1),2) + 16 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * (z - 1) + 16 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * z * (z - 1) - 16 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * z * (z - 1) - 16 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * z * (z - 1) - 16 * x * (x - 1) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 16 * x * (x - 1) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2)) + 5 * POW(z,4)

force = CoefficientFunction((f1,f2,f3))

u_1 = 2 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 2 * POW(x,2) * POW((x-1),2) * y*y* (y - 1) * POW(z,2) * POW((z-1),2) - 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) - 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1)
u_2 = -2 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 2 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) + 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1)
u_3 = 2 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 2 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 2 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2)

u_exact = CoefficientFunction( (u_1,u_2,u_3))

u_ex_deriv = CoefficientFunction( (4 * x * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * POW(x,2) * (x - 1) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * x * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) + 4 * POW(x,2) * (x - 1) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) - 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) - 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) - 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) - 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1),2 * POW(x,2) * POW((x-1),2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 8 * POW(x,2) * POW((x-1),2) * y * (y - 1) * POW(z,2) * POW((z-1),2) + 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * z * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * z * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * (z - 1) - 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * (z - 1),4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * z * POW((z-1),2) + 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * (z - 1) + 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * z * POW((z-1),2) + 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * (z - 1) - 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW((z-1),2) - 8 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * (z - 1) - 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2)))

v_ex_deriv = CoefficientFunction( (-2 * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 8 * x * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 2 * POW(x,2) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) + 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) + 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) + 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1),-4 * x * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * x * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * (x - 1) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * (x - 1) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) + 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * z * POW((z-1),2) + 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * z * POW((z-1),2) + 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * (z - 1) + 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * (z - 1),-4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) - 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) - 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) - 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) + 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW((z-1),2) + 8 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * (z - 1) + 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2)))

w_ex_deriv = CoefficientFunction( (2 * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 8 * x * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 2 * POW(x,2) * POW(y,2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * x * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * (x - 1) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 4 * x * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) - 4 * POW(x,2) * (x - 1) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2),4 * x * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * x * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) + 4 * POW(x,2) * (x - 1) * y * POW((y-1),2) * POW(z,2) * POW((z-1),2) + 4 * POW(x,2) * (x - 1) * POW(y,2) * (y - 1) * POW(z,2) * POW((z-1),2) - 2 * POW(x,2) * POW((x-1),2) * POW((y-1),2) * POW(z,2) * POW((z-1),2) - 8 * POW(x,2) * POW((x-1),2) * y * (y - 1) * POW(z,2) * POW((z-1),2) - 2 * POW(x,2) * POW((x-1),2) * POW(y,2) * POW(z,2) * POW((z-1),2),4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) + 4 * x * POW((x-1),2) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) + 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * z * POW((z-1),2) + 4 * POW(x,2) * (x - 1) * POW(y,2) * POW((y-1),2) * POW(z,2) * (z - 1) - 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * z * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * y * POW((y-1),2) * POW(z,2) * (z - 1) - 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * z * POW((z-1),2) - 4 * POW(x,2) * POW((x-1),2) * POW(y,2) * (y - 1) * POW(z,2) * (z - 1)))


order=2

#Taylor Hood
Vx = H1(mesh, order=order, dirichlet=[1,2,3,4,5,6])
Vy = H1(mesh, order=order, dirichlet=[1,2,3,4,5,6])
Vz = H1(mesh, order=order, dirichlet=[1,2,3,4,5,6])
Q = H1(mesh, order=order-1)

N2 = FESpace("number", mesh)

X = FESpace([Vx,Vy,Vz, Q,N2])

u,v,w, p, lam2 = X.TrialFunction()
wu,wv,ww,wp, mu2 = X.TestFunction()

gradu = u.Deriv()
gradv = v.Deriv()
gradw = w.Deriv()
gradwu = wu.Deriv()
gradwv = wv.Deriv()
gradww = ww.Deriv()

a = BilinearForm(X, symmetric=True)
a += SymbolicBFI( nu*(gradu*gradwu + gradv*gradwv + gradw*gradww) 
                  - (gradu[0] + gradv[1]+gradw[2]) * wp 
                  - (gradwu[0] + gradwv[1]+gradww[2]) * p + p*mu2 + wp*lam2)    

f = LinearForm(X)
f += SymbolicLFI( force[0]*wu + force[1]*wv + force[2]*ww, bonus_intorder = 10)

with TaskManager():
    a.Assemble()
    f.Assemble()

sol = GridFunction(X)

#Reconstruction:
Sig = FESpace("hdivho", mesh, order=order)
W = FESpace("l2hoh1dofs", mesh, order=order-1)
N = FESpace("number", mesh)
H1space = FESpace("h1ho", mesh, order=order-1)

Xr = FESpace([Sig, W,N])

sigma, w, lam = Xr.TrialFunction()
tau, vr, mu= Xr.TestFunction()

div_simga = sigma.Deriv()
div_tau = tau.Deriv()

ar = BilinearForm(Xr, symmetric = False)
ar += SymbolicBFI(tau*sigma + w*div_tau + vr*div_simga + w*mu + lam*vr)

fr = LinearForm(Xr)
fr += SymbolicLFI( force * tau, bonus_intorder = 10 )

fdelta = LinearForm(X)
with TaskManager():
    ar.Assemble()
    fr.Assemble()
    fdelta.Assemble()

with TaskManager():
    RecV = ReconstructionVertex3D()
    RecV.Setup(X,Xr,ar,H1space)
    RecV.CalcTrans(fr,fdelta, ar, H1space)
    f.vec.data -= fdelta.vec


with TaskManager():
    res = sol.vec.CreateVector()
    res.data = f.vec-a.mat*sol.vec 
    sol.vec.data += a.mat.Inverse(X.FreeDofs(),inverse="pardiso")*res
    
vel = CoefficientFunction( (sol.components[0], sol.components[1],sol.components[2]) )
pres = sol.components[3]
Draw(mesh=mesh,cf=vel,name="velocity")
Draw(mesh=mesh,cf=pres,name="pressure")

p_exact = CoefficientFunction(x*x*x*x*x+y*y*y*y*y+z*z*z*z*z - 1/2)

Draw(mesh=mesh,cf=p_exact,name="p_exakt")
Draw(mesh=mesh,cf=u_exact,name="u_exakt")

with TaskManager():
    err = (pres-p_exact)*(pres-p_exact)
    error = Integrate(err, mesh, VOL)
    print("L2-Error-pressure:", sqrt(error))

    err = (vel-u_exact)*(vel-u_exact)
    error = Integrate(err, mesh, VOL)
    print("L2-Error-vel:",sqrt(error))

    h1err = (sol.components[0].Deriv()-u_ex_deriv)*(sol.components[0].Deriv()-u_ex_deriv) + (sol.components[1].Deriv()-v_ex_deriv)*(sol.components[1].Deriv()-v_ex_deriv) + (sol.components[2].Deriv()-w_ex_deriv)*(sol.components[2].Deriv()-w_ex_deriv)
    h1error = Integrate(h1err, mesh, VOL)
    print("h1-Error-vel:",sqrt(abs(h1error)))


visoptions.scalfunction='velocity:0'
visoptions.clipsolution = 'scal'
   
viewoptions.clipping.nx= -1
viewoptions.clipping.ny= -1
viewoptions.clipping.nz= -1
viewoptions.clipping.enable = 1
visoptions.subdivisions = 4      
