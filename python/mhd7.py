import numpy as np
from numpy import gradient as grad
import numpy.matlib
import matplotlib.pyplot as plt
import pandas as pd 
#pd.read_csv()

def parx(f):
    return grad(f,dx,dy,edge_order = 2)[1]

def pary(f):
    return grad(f,dx,dy,edge_order = 2)[0]    
###### argument ######
nx = 801
ny = 201
lx = 2.
ly = 4.
x = np.linspace(-lx,lx,nx)
y = np.linspace(0,ly,ny)
dx = lx*2/(nx-1)
dx2 = dx*dx
dy = ly/(ny-1)
dy2 = dy*dy
eta = 1e-6
mu = 1e-6
gamma = 5./3.
beta = 0.5
b0 = 1.
l = 0.16
p00 = beta*0.5*b0**2
T0 = 0.25
Dn = 1e-5
Dp = 1e-5
dt = 1e-3
tstep = 100000
tic = 0
tin = 100
###### initial ######
By0 = np.matlib.repmat(b0*np.tanh(x/l),ny,1)
Bz0 = np.matlib.repmat(0.5,ny,nx)
Bx0 = np.matlib.repmat(0.*x,ny,1)
B20 = By0**2 + Bz0**2 + Bx0**2
psi0 = np.matlib.repmat(-l*b0*np.log(np.cosh(x/l)),ny,1) + 1e-3*np.dot(np.sin(.5*np.pi*y).reshape(ny,1),np.exp(-(x/l)**2).reshape(1,nx))
p0 = p00 + b0**2*.5 - .5*B20**2
p0 = np.matlib.repmat(np.tanh(-2*x)+3.5,ny,1)
n0 = p0/p0
vx0 = np.matlib.repmat(0.*x,ny,1)
vy0 = np.matlib.repmat(0.*x,ny,1)
vz0 = np.matlib.repmat(0.*x,ny,1)  
###### main ######
Bx,By,Bz,B2,psi,p,n,vx,vy,vz = Bx0,By0,Bz0,B20,psi0,p0,n0,vx0,vy0,vz0

for fk in np.arange(0,tstep):
    dBydx = parx(By)
    dBydy = pary(By)
    ddBy = parx(parx(By)) + pary(pary(By))

    dBxdx = parx(Bx)
    dBxdy = pary(Bx)
    ddBx = parx(parx(Bx)) + pary(pary(Bx))

    dBzdx = parx(Bz)
    dBzdy = pary(Bz)
    ddBz = parx(parx(Bz)) + pary(pary(Bz))

    dB2dx = parx(B2)
    dB2dy = pary(B2)

    dpsidx = parx(psi)
    dpsidy = pary(psi)
    ddpsi = parx(parx(psi)) + pary(pary(psi))

    dpdx = parx(p)
    dpdy = pary(p)
    ddp = parx(parx(p)) + pary(pary(p))

    dndx = parx(n)
    dndy = pary(n)
    ddn = parx(parx(n)) + pary(pary(n))

    dvydx = parx(vy)
    dvydy = pary(vy)
    ddvy = parx(parx(vy)) + pary(pary(vy))

    dvxdx = parx(vx)
    dvxdy = pary(vx)
    ddvx = parx(parx(vx)) + pary(pary(vx))

    dvzdx = parx(vz)
    dvzdy = pary(vz)
    ddvz = parx(parx(vz)) + pary(pary(vz))
####equation###
    dndt = - vx*dndx - vy*dndy - n*(dvxdx + dvydy) + Dn*ddn
    dvxdt = - vx*dvxdx - vy*dvxdy - dpdx*.5*beta/n - dB2dx*.5/n + (Bx*dBxdx + By*dBxdy)/n + ddvx*mu/n
    dvydt = - vx*dvydx - vy*dvydy - dpdy*.5*beta/n - dB2dy*.5/n + (Bx*dBydx + By*dBydy)/n + ddvy*mu/n
    dvzdt = - vx*dvzdx - vy*dvzdy + (Bx*dBzdx + By*dBzdy)/n + ddvz*mu/n
    dpdt = - vx*dpdx - vy*dpdy - (dvxdx+dvydy)*gamma*p + ddp*Dp
    dpsidt = - vx*dpsidx - vy*dpsidy + ddpsi*eta
    dBzdt = -vx*dBzdx -vy*dBzdy + Bx*dvzdx + By*dvzdy - Bz*dvxdx -Bz*dvydy + ddBz*eta
###RK###
    n = n + dndt*dt 
    vx = vx + dvxdt*dt 
    vy = vy + dvydt*dt 
    vz = vz + dvzdt*dt  
    p = p + dpdt*dt
    psi = psi + dpsidt*dt   
    Bz = Bz + dBzdt*dt    
    Bx = pary(psi)
    By = -parx(psi)
    B2 = Bx**2 + By**2 + Bz**2


    tic += 1
    print(tic)
    if np.mod(tic,tin) == 0:
        df = pd.DataFrame(psi) 
        df.to_csv('psi'+str(tic/tin)+'.csv')
        df = pd.DataFrame(p)
        df.to_csv('pre'+str(tic/tin)+'.csv')
        df = pd.DataFrame(Bz)
        df.to_csv('Bz'+str(tic/tin)+'.csv')   
        df = pd.DataFrame(vx)
        df.to_csv('vx'+str(tic/tin)+'.csv')  
        df = pd.DataFrame(vy)
        df.to_csv('vy'+str(tic/tin)+'.csv')  
        df = pd.DataFrame(vz)
        df.to_csv('vz'+str(tic/tin)+'.csv')  
    if abs(p[100,400]) >= 1e20:
        exit()      