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

nu = 1e-6

a = BilinearForm(X, symmetric=True)
a += SymbolicBFI( nu*(gradu*gradwu + gradv*gradwv) 
                  - (gradu[0] + gradv[1]) * wp 
                  - (gradwu[0] + gradwv[1]) * p + p*mu2 + wp*lam2)

a.Assemble()

#right hand side for: zeta = x^2(1-x)^2y^2(1-y)^2, u = curl(zeta); p = x^6+y^6 - 2/7
force = CoefficientFunction((-nu * (4 * (1 - x) * (1-x) * y * (1 - y) * (1-y) - 16 * x * (1 - x) * y * (1 - y) * (1-y) + 4 * x * x * y * (1 - y) * (1-y) - 4 * (1 - x) * (1-x) * y * y  * (1 - y) + 16 * x * (1 - x) * y * y  * (1 - y) - 4 * x * x * y * y * (1 - y) - 12 * x * x * (1 - x) * (1 - x) * (1 - y) + 12 * x * x  * (1 - x) * (1-x) * y) + 6*x * x * x * x * x,-nu * (12 * (1 - x) * y * y * (1 - y) * (1-y) - 12 * x * y * y * (1 - y) * (1-y) - 4 * x * (1 - x) * (1-x) * (1 - y) * (1-y) + 16 * x * (1 - x) * (1-x) * y * (1 - y) - 4 * x * (1 - x) * (1-x) * y * y  + 4 * x * x * (1 - x) * (1 - y) * (1-y) - 16 * x * x * (1 - x) * y * (1 - y) + 4 * x * x  * (1 - x) * y * y) + 6 *y* y *y*y*y))

f = LinearForm(X)
f += SymbolicLFI( force[0]*wu + force[1]*wv , bonus_intorder = 10)
f.Assemble()

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
ar.Assemble()

fr = LinearForm(Xr)
fr += SymbolicLFI( force * tau , bonus_intorder = 10)
fr.Assemble()

fdelta = LinearForm(X)
fdelta.Assemble()

#RecEl = ReconstructionElement()
#RecEl.Setup(X,Xr,ar,H1space)
#RecEl.CalcTrans(fr,fdelta)

RecV = ReconstructionVertex()
RecV.Setup(X,Xr,ar,H1space)
RecV.CalcTrans(fr,fdelta)

f.vec.data -= fdelta.vec

sol.vec.data = a.mat.Inverse(X.FreeDofs(),inverse="umfpack")*f.vec

vel = CoefficientFunction( (sol.components[0], sol.components[1]) )
pres = sol.components[2]

p_exact = CoefficientFunction(x*x*x*x*x*x+y*y*y*y*y*y - 2/7)
u_exact = CoefficientFunction( (x*x*(1-x)*(1-x) *(2*y*(1-y)*(1-y) - 2* y*y*(1-y)), -1*y*y*(1-y)*(1-y) *(2*x*(1-x)*(1-x) - 2*x*x*(1-x))))

u_ex_deriv = CoefficientFunction( (  (2*x*(1-x)*(1-x) - 2*x*x*(1-x)) * (2*y*(1-y)*(1-y) - 2* y*y*(1-y)) , x*x*(1-x)*(1-x) * (2*(1-y)*(1-y) - 8*y*(1-y)+2*y*y)  ))
v_ex_deriv = CoefficientFunction(  ( -y*y*(1-y)*(1-y) * (2*(1-x)*(1-x) - 8*x*(1-x)+2*x*x) , -1*(2*x*(1-x)*(1-x) - 2*x*x*(1-x)) * (2*y*(1-y)*(1-y) - 2* y*y*(1-y))  ))

Draw(mesh=mesh,cf=p_exact,name="p_exakt")
Draw(mesh=mesh,cf=u_exact,name="u_exakt")

Draw(mesh=mesh,cf=sol.components[1].Deriv(),name="v_deriv")
Draw(mesh=mesh,cf=v_ex_deriv,name="v_deriv_exakt")

Draw(mesh=mesh,cf=pres,name="pressure")
Draw(mesh=mesh,cf=vel,name="velocity")

err = (pres-p_exact)*(pres-p_exact)
error = Integrate(err, mesh, VOL)
print("L2-Error-pressure:", sqrt(error))

err = (vel-u_exact)*(vel-u_exact)
error = Integrate(err, mesh, VOL)
print("L2-Error-vel:",sqrt(error))

h1err = (sol.components[0].Deriv()-u_ex_deriv)*(sol.components[0].Deriv()-u_ex_deriv) + (sol.components[1].Deriv()-v_ex_deriv)*(sol.components[1].Deriv()-v_ex_deriv)
h1error = Integrate(h1err, mesh, VOL)
print("h1-Error-vel:",sqrt(abs(h1error)))

visoptions.scalfunction='velocity:0'
visoptions.subdivisions=4
