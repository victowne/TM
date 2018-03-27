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
from scipy.interpolate import griddata
import h5py
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

a = open('g048605.04686')
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


nbound = odata[280422:280425]
nlimiter = odata[280428:280430]

bound = wsnd(odata[280431:283832])
bound = np.array(bound[0:len(bound)],dtype = 'float64').reshape(105,2)

lim = wsnd(odata[283833:285971])
lim = np.array(lim[0:len(lim)-3],dtype = 'float64').reshape(66,2)


#funcpsi = interpolate.RectBivariateSpline(R,Z,psi)
qin = interp1d(corp,q,kind = 'cubic')
pin = interp1d(corp,p,kind = 'cubic')
# the = np.loadtxt('RZ_mid.dat')
# r,z = np.meshgrid(R,Z)
# psi[psi > Psi_bound] = np.nan
# thein = griddata(the[:,0:2],the[:,2],(r,z),method='nearest')
# Br =  parz(psi)
# Bz = -parr(psi)
# Jz = -parr(parr(psi)) - parz(parz(psi))
# Bphi = qin(psi)*(parz(psi)*parr(thein)-parr(psi)*parz(thein)) 
zb = np.loadtxt('RZ.dat')
zbt = np.ones(zb.shape)
nt = 1
zbt[0,:] = zb[0,:]
for i in np.arange(1,101):
    n = (i-1)*i*2+1
    for j in np.arange(1,i+1):
        zbt[nt,:] = zb[n,:]
        n += 1
        nt += 1

for i in np.arange(1,101):
    n = (i-1)*i*2+1+1*i
    for j in np.arange(1,i+1):
        zbt[nt,:] = zb[n,:]
        n += 1
        nt += 1
        
for i in np.arange(1,101):
    n = (i-1)*i*2+1+2*i
    for j in np.arange(1,i+1):
        zbt[nt,:] = zb[n,:]
        n += 1
        nt += 1
        
for i in np.arange(1,101):
    n = (i-1)*i*2+1+3*i
    for j in np.arange(1,i+1):
        zbt[nt,:] = zb[n,:]
        n += 1
        nt += 1
        
funcpsi = interpolate.RectBivariateSpline(R,Z,psi)
val = funcpsi(zbt[:,0],zbt[:,1],grid=False)
temp = np.ones((zbt.shape[0],zbt.shape[1]+1))
co = np.ones((1,3))
for i in np.arange(1,17):
    temp[:,0] = zbt[:,0]*np.cos((i-1)*0.125*np.pi)
    temp[:,1] = zbt[:,0]*np.sin((i-1)*0.125*np.pi)
    temp[:,2] = zbt[:,1]
    co = np.vstack((co,temp))
co = co[1:len(co),:]
for i in np.arange(1,5):
    val = np.hstack((val,val))
val.shape = (len(val),1)
f = h5py.File('test.h5')
del f['time_coordinates[0]/coordinates/values']
f['time_coordinates[0]/coordinates/values'] = co
del f['time_node_data[0]/node_data[0]/values']
f['time_node_data[0]/node_data[0]/values'] = val
del f['time_node_data[0]/node_data[1]/values']
f['time_node_data[0]/node_data[1]/values'] = val
del f['time_node_data[0]/node_data[2]/values']
f['time_node_data[0]/node_data[2]/values'] = val
del f['time_node_data[0]/node_data[3]/values']
f['time_node_data[0]/node_data[3]/values'] = val
del f['time_node_data[0]/node_data[4]/values']
f['time_node_data[0]/node_data[4]/values'] = val
del f['time_node_data[0]/node_data[5]/values']
f['time_node_data[0]/node_data[5]/values'] = val
del f['time_node_data[0]/node_data[6]/values']
f['time_node_data[0]/node_data[6]/values'] = val
del f['time_node_data[0]/node_data[7]/values']
f['time_node_data[0]/node_data[7]/values'] = val
del f['time_node_data[0]/node_data[8]/values']
f['time_node_data[0]/node_data[8]/values'] = val
del f['time_node_data[0]/node_data[9]/values']
f['time_node_data[0]/node_data[9]/values'] = val
del f['time_node_data[0]/node_data[10]/values']
f['time_node_data[0]/node_data[10]/values'] = val
del f['time_node_data[0]/node_data[11]/values']
f['time_node_data[0]/node_data[11]/values'] = val
del f['time_node_data[0]/node_data[12]/values']
f['time_node_data[0]/node_data[12]/values'] = val
del f['time_node_data[0]/node_data[13]/values']
f['time_node_data[0]/node_data[13]/values'] = val
f.close()



#custom1
val = funcpsi(zbt[:,0],zbt[:,1],grid=False)
val[0] = 1
n = 1
for k in np.arange(1,5):
    for i in np.arange(2,102):
        for j in np.arange(1,i):
            val[n] = i
            n += 1
for i in np.arange(1,5):
    val = np.hstack((val,val))
val.shape = (len(val),1)

#custom2
val = np.ones(323216)
n = 0
for k in np.arange(1,17):
    val[n] = corp[0]
    n += 1
    to = 0.125*np.pi*(k-1)
    for i in np.arange(2,102):
        for j in np.arange(1,i):
            po = np.pi*j/(2*i)
            val[n] = corp[i-1]*np.cos(2*po+1*to)
            n += 1
   
    for i in np.arange(2,102):
        for j in np.arange(1,i):
            po = np.pi*j/(2*i)+.5*np.pi
            val[n] = corp[i-1]*np.cos(2*po+1*to)
            n += 1
  
    for i in np.arange(2,102):
        for j in np.arange(1,i):
            po = np.pi*j/(2*i)+np.pi
            val[n] = corp[i-1]*np.cos(2*po+1*to)
            n += 1

    for i in np.arange(2,102):
        for j in np.arange(1,i):
            po = np.pi*j/(2*i)+1.5*np.pi
            val[n] = corp[i-1]*np.cos(2*po+1*to)
            n += 1
    
val.shape = (len(val),1)

#costom3
val[0] = 1
n = 1
for k in np.arange(1,5):
    for i in np.arange(2,102):
        if (i == 50 or i ==10): 
            for j in np.arange(1,i):
                val[n] = i
                n += 1
        else:
             for j in np.arange(1,i):
                val[n] = 0
                n += 1   
for i in np.arange(1,5):
    val = np.hstack((val,val))
val.shape = (len(val),1)





val = np.ones(323216)
n = 0
for k in np.arange(1,17):
    val[n] = 0*corp[0]
    n += 1
    to = 0.125*np.pi*(k-1)
    for i in np.arange(2,102):
        if i == 55:
            for j in np.arange(1,i):
                po = np.pi*j/(2*i)
                val[n] = corp[i-1]*np.cos(3*po+2*to)
                n += 1
        elif i == 71:
            for j in np.arange(1,i):
                po = np.pi*j/(2*i)
                val[n] = corp[i-1]*np.cos(2*po+1*to)
                n += 1
        else:
            for j in np.arange(1,i):
             
                val[n] = 0
                n += 1          
    for i in np.arange(2,102):
        if i == 55:
            for j in np.arange(1,i):
                po = np.pi*j/(2*i) + .5*np.pi
                val[n] = corp[i-1]*np.cos(3*po+2*to)
                n += 1
        elif i == 71:
            for j in np.arange(1,i):
                po = np.pi*j/(2*i) + .5*np.pi
                val[n] = corp[i-1]*np.cos(2*po+1*to)
                n += 1
        else:
            for j in np.arange(1,i):
             
                val[n] = 0
                n += 1        
    for i in np.arange(2,102):
        if i == 55:
            for j in np.arange(1,i):
                po = np.pi*j/(2*i) + 1*np.pi
                val[n] = corp[i-1]*np.cos(3*po+2*to)
                n += 1
        elif i == 71:
            for j in np.arange(1,i):
                po = np.pi*j/(2*i) + 1*np.pi
                val[n] = corp[i-1]*np.cos(2*po+1*to)
                n += 1
        else:
            for j in np.arange(1,i):
             
                val[n] = 0
                n += 1        
    for i in np.arange(2,102):
        if i == 55:
            for j in np.arange(1,i):
                po = np.pi*j/(2*i) + 1.5*np.pi
                val[n] = corp[i-1]*np.cos(3*po+2*to)
                n += 1
        elif i == 71:
            for j in np.arange(1,i):
                po = np.pi*j/(2*i) + 1.5*np.pi
                val[n] = corp[i-1]*np.cos(2*po+1*to)
                n += 1
        else:
            for j in np.arange(1,i):
             
                val[n] = 0
                n += 1        