import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(0,7,200)
y = np.ones(200)
for i in range(0,200):
        y[i] = np.sin(x[i])
plt.plot(x,y)
plt.show()
