import matplotlib as mpl
import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import style
import os

#print style.available      
plt.style.use('seaborn-paper')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.major.size'] = 0
mpl.rcParams['ytick.major.size'] = 0
os.chdir('C:\\Users\\vic_l\\Desktop\\p3e5e5n1e4o3e2fb0.20kp10kc3e6_10.07e4')
#Arguments
Imax = 499
L = 1
Lmax = 3
rs = 310 
drs = rs+L*(Imax+1)
tlast = 4000
tstep = 50
ttotal = (tlast-1000)*tstep
t = np.linspace(0,ttotal,tlast-999)
p = np.loadtxt('Profile0.dat')
c = abs(p[rs,1]/(p[rs,2]*p[rs,5]))

ww=[]
phase=[]
for i in range(1000,tlast+1):
    A = np.loadtxt('xy'+str(i)+'.dat')
    ww.append(np.sqrt(np.sqrt(np.square(A[drs-1,4])+np.square(A[drs-1,5]))*c)*4)
    a = A[:,4]
    b = A[:,5]
    a = a.reshape(Imax+1,Lmax+1,order='F')
    b = b.reshape(Imax+1,Lmax+1,order='F')
    a = a+1j*b
    phase.append(np.angle(a[rs-1,1]))

omega = np.ones(tlast-999)
for i in range(0,tlast-1000):
    if phase[i+1]-phase[i] < 4.5:
        omega[i] = (phase[i+1]-phase[i])/tstep
    else:
        omega[i] = (phase[i+1]-phase[i]-np.pi*2)/tstep
omega[tlast-1000] = omega[tlast-1001]

plt.subplot(211)
plt.plot(t,ww,'-',label='$w{_2}{_/}{_1}$')
#plt.semilogy(a[:,0],a[:,12],'-',label='$E{_3}{_/}{_2}$')
plt.grid(linestyle = '--')
plt.title(r'$Island$ $Width$&$Mode$ $Frequence$',fontsize=14)
#plt.xlabel(r'$t/\tau_a$',fontsize=13)
#plt.xlim((0,1))
#plt.xticks(np.linspace(-1, 1, 5))
plt.ylabel(r'$w/a$',fontsize=13)
#plt.ylim((0,1))
#plt.yticks([0, 0.5], ['$minimum$', 'normal'])
#plt.autoscale(tight=True)   
plt.legend(loc='upper right')
plt.subplot(212)
plt.plot(t,-omega,'r-',label='$\omega{_2}{_/}{_1}$')
#plt.semilogy(a[:,0],a[:,12],'-',label='$E{_3}{_/}{_2}$')
plt.grid(linestyle = '--')
#plt.title(r'$Mode$ $Frequence$',fontsize=14)
plt.xlabel(r'$t/\tau_a$',fontsize=13)
#plt.xlim((0,1))
#plt.xticks(np.linspace(-1, 1, 5))
plt.ylabel(r'$\omega$',fontsize=13)
#plt.ylim((0,1))
#plt.yticks([0, 0.5], ['$minimum$', 'normal'])
#plt.autoscale(tight=True)   
plt.legend(loc='upper right')

plt.savefig("hehe.png")
plt.show()