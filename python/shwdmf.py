import matplotlib as mpl
import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import style
import os
import pandas as pd

plt.switch_backend('agg')
#print style.available      
plt.style.use('seaborn-paper')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.major.size'] = 0
mpl.rcParams['ytick.major.size'] = 0
mpl.rcParams['mathtext.default'] = 'regular'
font = {'family': 'Times New Roman',
         'style': 'italic',
        'weight': 'normal',
         'color': 'black', 
          'size': 18,
        }

os.chdir('/home/wktang/testfb/p3e5e5n1e4o3e2fb0.25kp10kc1e6')
#Arguments
Imax = 499
L = 1
Lmax = 3
rs = 310 
drs = rs+L*(Imax+1)
tlast = 2000
tstep = 50
ttotal = (tlast-1000)*tstep
t = np.linspace(0,ttotal,tlast-999)
p = np.loadtxt('Profile0.dat')
c = abs(p[rs,1]/(p[rs,2]*p[rs,5]))
rmpmax = 9.9
ts = (rmpmax*1000+10000)/50

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

df = pd.DataFrame(ww) 
df.to_csv('width.csv') 
df = pd.DataFrame(omega) 
df.to_csv('omega.csv') 

rmp = np.ones(tlast-999)
for i in range(0,tlast-999):
    if i < 200:
        rmp[i] = 0
    elif 200 <= i <= ts:
        rmp[i] = 1e-3*(i*50-10000)
    else:
        rmp[i] = rmpmax

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(t,ww,'-',label='$w{_2}{_/}{_1}$')
plt.grid(linestyle = '--')
ax2 = ax1.twinx()
ax2.plot(t,rmp,'-',label = 'RMP',color = '#ffaa00ff')
#plt.semilogy(a[:,0],a[:,12],'-',label='$E{_3}{_/}{_2}$')
plt.title('Island Width&Mode Frequence',fontdict = font)
#plt.xlabel(r'$t/\tau_a$',fontsize=13)
#plt.xlim((0,1))
#plt.xticks(np.linspace(-1, 1, 5))
ax1.set_ylabel('w/a',fontdict = font)
#ax2.set_ylim([-0.8,10])
ax2.set_ylabel(r'rmp($10^{-4}$)',fontdict = font)
#plt.ylim((0,1))
#plt.yticks([0, 0.5], ['$minimum$', 'normal'])
#plt.autoscale(tight=True)   
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')
plt.subplot(212)
plt.plot(t,-omega,'r-',label=r'$\omega{_2}{_/}{_1}$')
#plt.semilogy(a[:,0],a[:,12],'-',label='$E{_3}{_/}{_2}$')
plt.grid(linestyle = '--')
#plt.title(r'$Mode$ $Frequence$',fontsize=14)
plt.xlabel(r'$t/\tau_a$',fontdict = font)
#plt.xlim((0,1))
#plt.xticks(np.linspace(-1, 1, 5))
plt.ylabel(r'$\omega$',fontdict = font)
#plt.ylim((0,1))
#plt.yticks([0, 0.5], ['$minimum$', 'normal'])
#plt.autoscale(tight=True)   
plt.legend(loc='upper left')

#plt.savefig("C:\\Users\\Administrator\\Desktop\\hehe.png")
plt.savefig("2.png")
# plt.show()
