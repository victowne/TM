import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import style
import os

#print style.available      
plt.style.use('seaborn-paper')

os.chdir('C:\\Users\\vic_l\\Desktop')
a = np.loadtxt('enrdps.dat')
plt.semilogy(a[:,0],a[:,5],'-',label='$E{_2}{_/}{_1}$')
plt.semilogy(a[:,0],a[:,16],'-',label='$E{_3}{_/}{_2}$')


#plt.grid()
plt.title(r'$\mu=1e{^-}{^4}, \eta=5e{^-}{^5}$',fontsize=14)
plt.xlabel(r'$t/\tau_a$',fontsize=13)
#plt.xlim((0,1))
#plt.xticks(np.linspace(-1, 1, 5))
plt.ylabel(r'$E{_m}{_a}{_g}$',fontsize=13)
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