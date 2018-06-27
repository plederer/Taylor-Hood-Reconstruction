from ngsolve import *
from netgen.geom2d import unit_square
from ngsolve.internal import visoptions
from ngsolve.internal import viewoptions

from reclib import *

mesh = Mesh(unit_square.GenerateMesh (maxh=0.2))

order=3

#Taylor Hood
Vx = H1(mesh, order=order, dirichlet=[1,2,3,4])
Vy = H1(mesh, order=order, dirichlet=[1,2,3,4])
Q = H1(mesh, order=order-1)

N2 = FESpace("number", mesh)

X = FESpace([Vx,Vy,Q,N2])

u,v,p, lam2 = X.TrialFunction()
wu,wv,wp, mu2 = X.TestFunction()

gradu = u.Deriv()
gradv = v.Deriv()
gradwu = wu.Deriv()
gradwv = wv.Deriv()

nu = 1

a = BilinearForm(X, symmetric=True)
a += SymbolicBFI( nu*(gradu*gradwu + gradv*gradwv) 
                  - (gradu[0] + gradv[1]) * wp 
                  - (gradwu[0] + gradwv[1]) * p + p*mu2 + wp*lam2)



#right hand side for: zeta = x^2(1-x)^2y^2(1-y)^2, u = curl(zeta); p = x^6+y^6 - 2/7
force = CoefficientFunction((-nu * (4 * (1 - x) * (1-x) * y * (1 - y) * (1-y) - 16 * x * (1 - x) * y * (1 - y) * (1-y) + 4 * x * x * y * (1 - y) * (1-y) - 4 * (1 - x) * (1-x) * y * y  * (1 - y) + 16 * x * (1 - x) * y * y  * (1 - y) - 4 * x * x * y * y * (1 - y) - 12 * x * x * (1 - x) * (1 - x) * (1 - y) + 12 * x * x  * (1 - x) * (1-x) * y) + 6*x * x * x * x * x,-nu * (12 * (1 - x) * y * y * (1 - y) * (1-y) - 12 * x * y * y * (1 - y) * (1-y) - 4 * x * (1 - x) * (1-x) * (1 - y) * (1-y) + 16 * x * (1 - x) * (1-x) * y * (1 - y) - 4 * x * (1 - x) * (1-x) * y * y  + 4 * x * x * (1 - x) * (1 - y) * (1-y) - 16 * x * x * (1 - x) * y * (1 - y) + 4 * x * x  * (1 - x) * y * y) + 6 *y* y *y*y*y))

f = LinearForm(X)
f += SymbolicLFI( force[0]*wu + force[1]*wv , bonus_intorder = 10)


sol = GridFunction(X)

sol.components[0].Set(0*x,BND)
sol.components[1].Set(0*y,BND)



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

ar = BilinearForm(Xr, symmetric = False) #symmetric = False is important!!!
ar += SymbolicBFI(tau*sigma + w*div_tau + vr*div_simga + w*mu + lam*vr)

sig = GridFunction(Xr)

with TaskManager():
    a.Assemble()
    f.Assemble()
    ar.Assemble()
    
    sol.vec.data = a.mat.Inverse(X.FreeDofs(),inverse="umfpack")*f.vec
    RecV = ReconstructionVertex()
    RecV.Setup(X,Xr,ar,H1space)
    RecV.Calc(sol,sig)

vel = CoefficientFunction( (sol.components[0], sol.components[1]) )
velR = vel - sig.components[0]

div_vel = sol.components[0].Deriv()[0] + sol.components[1].Deriv()[1]
div_velR = div_vel - sig.components[0].Deriv()

Draw(vel, mesh, "velocity")
Draw(velR, mesh, "R_velocity")
Draw(div_vel, mesh, "div_vel")
Draw(div_velR, mesh, "div_velR")

pres = sol.components[2]

visoptions.scalfunction='velocity:0'
visoptions.subdivisions=4
