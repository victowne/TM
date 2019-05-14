import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np

a = np.loadtxt('-01.40625-p.txt')
x = a[:,0]
y = a[:,1]
z = a[:,2]

triang = tri.Triangulation(x, y)
mask = np.where(x[triang.triangles].max(axis=1)-x[triang.triangles].min(axis=1)>3, 1, 0)
mask2 = np.where(y[triang.triangles].max(axis=1)-y[triang.triangles].min(axis=1)>4, 1, 0)
mask = mask+mask2
triang.set_mask(mask)

plt.tricontourf(triang, z, 100, cmap='hot')
plt.axis('equal')
plt.colorbar()
plt.show()