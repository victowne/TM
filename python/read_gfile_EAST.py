# !/usrs/bin/env python
# coding: utf-8
"""
made by  tangweikang & mahaojie
time : 2018 10 23
purpose : read gfile and construct bfiled surface coordinate
"""
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy import interpolate
from numpy import gradient as grad
from scipy.interpolate import interp1d
import numpy.matlib
def wsnd(data):
    datas = data.split('\n')
    dataa = []
    for fk in datas:
        temp = fk[0:16]
        dataa.append(temp)
        temp = fk[16:32]
        dataa.append(temp)
        temp = fk[32:48]
        dataa.append(temp)
        temp = fk[48:64]
        dataa.append(temp)
        temp = fk[64:80]
        dataa.append(temp)
    return dataa

def parr(f):
    return grad(f,dR,dZ,edge_order = 2)[1]

def parz(f):
    return grad(f,dR,dZ,edge_order = 2)[0]  

a = open('g022485.00580')
odata = a.read()
nr = int(odata[53:56])
nz = int(odata[57:60])
Rboxlen = odata[62:77]
Zboxlen = odata[78:93]
R0= odata[94:109]
Rmin = odata[110:125]
Z0 = odata[126:141]
#Raxis = odata[142:158]
#Zaxis = odata[158:174]
Psi_axis = float(odata[174:190])
Psi_bound = float(odata[190:206])
corp = np.linspace(Psi_axis,Psi_bound,129)
#B0 = odata[207:222]
#current = odata[224:239]
#Psi_axis1 = odata[239:255]
#xdum1 = odata[256:271]
#Raxis = odata[272:287]
#xdum2 = odata[288:303]
#Zaxis = odata[304:320]
#xdum3 = odata[320:336]
#Psi_bound = odata[336:352]
#xdum4 = odata[353:368]
#xdum5 = odata[369:384]
f = np.array(odata[385:2475].split(),dtype = 'float64')
p = np.array(odata[2475:4565].split(),dtype = 'float64')
ffprime = np.array(odata[4565:6655].split(),dtype = 'float64')
pprime = re.split('-| ',odata[6655:8745])
pprime = -np.array(pprime[1:len(pprime)],dtype = 'float64')
psi = wsnd(odata[8745:278329])
psi = np.array(psi[0:len(psi)-4],dtype = 'float64').reshape(nr,nz)
q = np.array(odata[278329:280419].split(),dtype = 'float64')
R = np.ones(nr)
Z = np.ones(nz)
for i in np.arange(0,129):
    R[i] = float(Rmin)+(float(Rboxlen)*i)/128
    Z[i] = float(Z0)-float(Zboxlen)*0.5+(float(Zboxlen)*i)/128
dR = R[1]-R[0]
dZ = Z[1]-Z[0]
Br =  parz(psi)
Bz = -parr(psi)
Jz = -parr(parr(psi)) - parz(parz(psi))
#---------------EAST
# nbound = odata[280422:280425]
# nlimiter = odata[280428:280430]

# bound = wsnd(odata[280431:283832])
# bound = np.array(bound[0:len(bound)],dtype = 'float64').reshape(105,2)

# lim = wsnd(odata[283833:285971])
# lim = np.array(lim[0:len(lim)-3],dtype = 'float64').reshape(66,2)
#---------------HL2A
bound = wsnd(odata[280431:283120])
bound = np.array(bound[0:len(bound)-4],dtype = 'float64').reshape(83,2)

lim = wsnd(odata[283121:284481])
lim = np.array(lim[0:len(lim)-1],dtype = 'float64').reshape(42,2)
#PSI
plt.contourf(R,Z,psi,100,alpha = 0.75,cmap="Blues_r")
plt.plot(bound[:,0],bound[:,1],label = 'boundary')
plt.plot(lim[:,0],lim[:,1],label = 'limiter')
plt.title('$\psi$  $countour$',fontsize=14)
plt.xlabel('R/m',fontsize=13)
plt.ylabel('Z/m',fontsize=13)
plt.legend()
plt.axis('equal')
plt.show()
#safety factor
plt.plot(corp,q)
plt.title('qprofile',fontsize=14)
plt.xlabel('$\Psi$',fontsize=13)
plt.ylabel('q',fontsize=13)
plt.show()
#f
plt.plot(corp,f)
plt.title('$Poloidal$ $Current$',fontsize=14)
plt.xlabel('$\Psi$',fontsize=13)
plt.ylabel('f/(T*m)',fontsize=13)
plt.show()
#pressure
plt.plot(corp,p)
plt.title('$Pressure$ $Profile$',fontsize=14)
plt.xlabel('$\Psi$',fontsize=13)
plt.ylabel('P/(Pa)',fontsize=13)
plt.show()
#pp'
plt.plot(corp,pprime)
plt.title('$p\'$',fontsize=14)
plt.xlabel('$\Psi$',fontsize=13)
plt.ylabel('dP/(Pa/T*m^2)',fontsize=13)
plt.show()
#ff'
plt.plot(corp,ffprime)
plt.title('$ff\'$',fontsize=14)
plt.xlabel('$\Psi$',fontsize=13)
plt.ylabel('ffpeimr/(T)',fontsize=13)
plt.show()
#Bvector
# plt.contour(R,Z,psi,30,cmap = 'Blues_r')
# plt.stream(R,Z,Br,Bz,units='x', pivot='tip', width=0.0022,scale=1 / 0.05)
plt.streamplot(R,Z,Br,Bz,density = 1.5,linewidth = 1.,arrowsize=0.5, arrowstyle='simple')
plt.title('$mag field$',fontsize=14)
plt.xlabel('R/m',fontsize=13)
plt.ylabel('Z/m',fontsize=13)
plt.axis('equal')
plt.show()
#Current
# plt.contour(R,Z,psi,30)
plt.contourf(R,Z,Jz,30,cmap = 'hot')
plt.streamplot(R,Z,Br,Bz,density = 1.5,linewidth = 1.,arrowsize=0.5, arrowstyle='simple')
plt.title('$J$',fontsize=14)
plt.xlabel('R/m',fontsize=13)
plt.ylabel('Z/m',fontsize=13)
plt.axis('equal')
plt.show()

#funcpsi = interpolate.RectBivariateSpline(R,Z,psi)
qin = interp1d(corp,q,kind = 'cubic')
pin = interp1d(corp,p,kind = 'cubic')
psi[psi > Psi_bound] = Psi_bound
plt.contourf(R,Z,qin(psi),20,cmap = 'Oranges_r')
plt.axis('equal')
plt.colorbar()
plt.show()
plt.contour(R,Z,pin(psi),15,cmap = 'hot')
plt.axis('equal')
plt.colorbar()
plt.show()