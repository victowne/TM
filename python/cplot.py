import matplotlib.pyplot as plt 
import numpy as np 
import matplotlib as mpl
import numpy.matlib
from matplotlib import style
  ######## Debug #############
#   import pandas as pd       #
#   df = pd.DataFrame(in_cnt) #
#   df.to_csv('a.txt')        #
  ############################
#def cart2pol(x, y):
#    rho = np.sqrt(x**2 + y**2)
#    phi = np.arctan2(y, x)
#    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)
#print style.available      
plt.style.use('seaborn-paper')
mpl.rcParams['xtick.major.size'] = 0
mpl.rcParams['ytick.major.size'] = 0

Lmax=3               
Imax=499
m=3

A = np.loadtxt('xy1600.dat')
p = np.loadtxt('Profile0.dat')
x = p[:,0]
eq = p[:,3]
y = np.linspace(0,2.*np.pi,630)
Y,X = np.meshgrid(y,x)
a = A[:,4]
b = A[:,5]
a = a.reshape(Imax+1,Lmax+1,order='F')
b = b.reshape(Imax+1,Lmax+1,order='F')
a = a+1j*b
k = m*np.arange(0,Lmax+1,1)
k = k.reshape(-1,1)      
b = 2*np.exp(1j*k*y)
b[0,:] = b[0,:]/2
pfunc = np.dot(a,b)
pfunc = np.real(pfunc)
pfunc0 = eq+.25*x*x/3
pfunc0 = pfunc0.reshape(-1,1)
pfunc0 = np.matlib.repmat(pfunc0,1,630)
pfunc = pfunc+pfunc0


X,Y = pol2cart(X,Y)
plt.contourf(Y, X, pfunc, 30, alpha=.99, cmap='plasma')
#plt.colorbar()
#plt.contour(Y, X, pfunc, 30, colors='gray', linewidth=.5)



#plt.grid()
plt.title(r'$\psi^*$',fontsize=14)
plt.xlabel(r'$\theta/rad$',fontsize=13)
#plt.xlim((0,1))
#plt.xticks(np.linspace(-1, 1, 5))
plt.ylabel(r'$r/a$',fontsize=13)
#plt.ylim((0,1))
#plt.yticks([0, 0.5], ['$minimum$', 'normal'])
#plt.autoscale(tight=True)   




#plt.savefig("C:\\Users\\vic_l\\Desktop\\test.png")
plt.show()