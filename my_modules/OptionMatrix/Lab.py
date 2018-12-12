## TODO 试验区，待删除
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(25, 150, 2.5)
Y = np.arange(0.02, 0.5, 0.02)
#网格化数据
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(delta)
Z = vfunc(1.0,X,100.0,0.07,Y,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()