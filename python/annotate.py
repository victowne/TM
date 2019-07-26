import matplotlib.pyplot as plt 
import numpy as np 
import matplotlib as mpl
import numpy.matlib
from matplotlib import style
import h5py
#plt.switch_backend('agg')
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
mpl.rcParams['mathtext.default'] = 'regular'
font = {'family': 'Times New Roman',
         'style': 'italic',
        'weight': 'normal',
         'color': 'black', 
          'size': 18,
        }
font2 = {'family': 'Times New Roman',
         'style': 'normal',
        'weight': 'normal',
         'color': 'black', 
          'size': 18,
        }
font3 = {'family': 'Times New Roman',
         'style': 'italic',
        'weight': 'normal',
         'color': 'darkorange', 
          'size': 18,
        }

Lmax = 6               
Imax = 1000
m = 3
tstep = 50

ttotal = 15000
tstep = 20
tn = ttotal/tstep
t = np.linspace(0,ttotal,tn)
f = h5py.File('data.h5')
w = f['/width']
p = f['/phase']
omega = np.ones(tn)
for i in range(0,tn-1):
    if p[i+1]-p[i] < 2.:
        omega[i] = (p[i+1]-p[i])/tstep
    else:
        omega[i] = (p[i+1]-p[i]-np.pi*2)/tstep
omega[tn-1] = omega[tn-2]
plt.subplot(2,1,1)
plt.plot(t,omega[0:tn])
plt.plot([0,15000],[-9.6e-2,-9.6e-2],'--')
plt.plot([2000,2000],[0.02,-0.23],'k--',linewidth=1)
plt.plot([3900,3900],[0.02,-0.23],'k--',linewidth=1)
plt.plot([6030,6030],[0.02,-0.23],'k--',linewidth=1)
plt.annotate('RMP on',fontsize=10,xy=(2000,-0.1),xytext=(200,-0.17),arrowprops=dict(arrowstyle='<-',connectionstyle='arc3',lw='1'))
plt.text(480,-0.025,'I',fontdict=font2)
plt.text(2100,-0.025,'II',fontdict=font2)
plt.text(4000,-0.025,'III',fontdict=font2)
plt.text(6500,-0.025,'IV',fontdict=font2)
plt.text(13000,-0.13,r'$\omega_{rmp}$',fontdict=font3)
plt.xlabel(r'$t/\tau_{a}$',fontdict=font)
plt.ylabel(r'$\omega_{ntm}$',fontdict=font)
plt.xlim((0,15000))
plt.ylim((-0.23,0.02))
plt.tick_params(labelsize=11)

p = np.loadtxt('Profile0.dat')
x = p[:,0]
eq = p[:,3]
y = np.linspace(0,2.*np.pi,630)
Y,X = np.meshgrid(y,x)
X,Y = pol2cart(X,Y)
plt.subplot(2,3,4)
for i in range(2200,2201,tstep):
	a = np.transpose(np.array(f['t='+str(i)+'/psiR']))
	b = np.transpose(np.array(f['t='+str(i)+'/psiI']))
	a = a+1j*b
	k = m*np.arange(0,Lmax+1,1)
	k = k.reshape(-1,1)      
	b = 2*np.exp(1j*k*y)
	b[0,:] = b[0,:]/2
	pfunc = np.dot(a,b)
	pfunc = np.real(pfunc)
	pfunc0 = eq+.25*x*x/6
	pfunc0 = pfunc0.reshape(-1,1)
	pfunc0 = np.matlib.repmat(pfunc0,1,630)
	pfunc = pfunc+pfunc0
	plt.contourf(Y[:650,:], X[:650,:], pfunc[:650,:], 80, alpha=.99, cmap='hot')
	plt.contour(Y[:650,:], X[:650,:], pfunc[:650,:], 13, colors='gray',linewidths=0.2)
	#plt.contour(Y, X, pfunc, 30, colors='gray', linewidth=.5)
#	plt.title(r't='+str(i),fontdict=font)
	plt.text(-0.6,0.45,'II',fontdict=font2)
	plt.xlabel(r'$r/a$',fontdict=font)
	plt.ylabel(r'$z/a$',fontdict=font)
	plt.tick_params(labelsize=11)
	plt.gca().axis('square')
plt.subplot(2,3,5)
for i in range(4500,4501,tstep):
	a = np.transpose(np.array(f['t='+str(i)+'/psiR']))
	b = np.transpose(np.array(f['t='+str(i)+'/psiI']))
	a = a+1j*b
	k = m*np.arange(0,Lmax+1,1)
	k = k.reshape(-1,1)      
	b = 2*np.exp(1j*k*y)
	b[0,:] = b[0,:]/2
	pfunc = np.dot(a,b)
	pfunc = np.real(pfunc)
	pfunc0 = eq+.25*x*x/6
	pfunc0 = pfunc0.reshape(-1,1)
	pfunc0 = np.matlib.repmat(pfunc0,1,630)
	pfunc = pfunc+pfunc0
	plt.contourf(Y[:650,:], X[:650,:], pfunc[:650,:], 80, alpha=.99, cmap='hot')
	plt.contour(Y[:650,:], X[:650,:], pfunc[:650,:], 13, colors='gray',linewidths=0.2)
	#plt.contour(Y, X, pfunc, 30, colors='gray', linewidth=.5)
#	plt.title(r't='+str(i),fontdict=font)
	plt.text(-0.6,0.45,'III',fontdict=font2)
	plt.xlabel(r'$r/a$',fontdict=font)
	plt.tick_params(labelsize=11)
	plt.gca().axis('square')
plt.subplot(2,3,6)
for i in range(12000,12001,tstep):
	a = np.transpose(np.array(f['t='+str(i)+'/psiR']))
	b = np.transpose(np.array(f['t='+str(i)+'/psiI']))
	a = a+1j*b
	k = m*np.arange(0,Lmax+1,1)
	k = k.reshape(-1,1)      
	b = 2*np.exp(1j*k*y)
	b[0,:] = b[0,:]/2
	pfunc = np.dot(a,b)
	pfunc = np.real(pfunc)
	pfunc0 = eq+.25*x*x/6
	pfunc0 = pfunc0.reshape(-1,1)
	pfunc0 = np.matlib.repmat(pfunc0,1,630)
	pfunc = pfunc+pfunc0
	plt.contourf(Y[:650,:], X[:650,:], pfunc[:650,:], 80, alpha=.99, cmap='hot')
	plt.contour(Y[:650,:], X[:650,:], pfunc[:650,:], 13, colors='gray',linewidths=0.2)
	#plt.contour(Y, X, pfunc, 30, colors='gray', linewidth=.5)
#	plt.title(r't='+str(i),fontdict=font)
	plt.text(-0.6,0.45,'IV',fontdict=font2)
	plt.xlabel(r'$r/a$',fontdict=font)
	plt.tick_params(labelsize=11)
	plt.gca().axis('square')
	plt.show()
	#plt.clf()
