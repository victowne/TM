import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import style

#print style.available      
plt.style.use('seaborn-paper')

#Arguments
Imax = 499
Lmax = 3
rs = 310 
tlast = 1380
tstep = 50
ttotal = (tlast-1000)*tstep
t = np.linspace(0,ttotal,tlast-999)

phase=[]
for i in range(1000,tlast+1):
    A = np.loadtxt('xy'+str(i)+'.dat')
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
        

plt.plot(t,-omega,'-',label='$\omega{_2}{_/}{_1}$')
#plt.semilogy(a[:,0],a[:,12],'-',label='$E{_3}{_/}{_2}$')


#plt.grid()
plt.title(r'$Mode$ $Frequence$',fontsize=14)
plt.xlabel(r'$t/\tau_a$',fontsize=13)
#plt.xlim((0,1))
#plt.xticks(np.linspace(-1, 1, 5))
plt.ylabel(r'$\omega$',fontsize=13)
#plt.ylim((0,1))
#plt.yticks([0, 0.5], ['$minimum$', 'normal'])
#plt.autoscale(tight=True)   
plt.legend()

ax = plt.gca()
#ax.spines['left'].set_color(16)
#ax.spines['right'].set_color('red')
#ax.spines['top'].set_color('red')
#ax.spines['bottom'].set_color(16)
#ax.xaxis.set_ticks_position('bottom')
#ax.yaxis.set_ticks_position('left')
#ax.set_fontsize(10)

#plt.savefig("C:\\Users\\vic_l\\Desktop\\test.png")
plt.show()