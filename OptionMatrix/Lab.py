## TODO 试验区，待删除
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Delta与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(25.0, 150.0, 2.5)
Y = np.arange(0.02, 0.5, 0.02)
#网格化数据
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(Delta)
Z = vfunc(0.0,X,100.0,0.07,Y,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# DdeltaDvol与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(50.0 ,150.0, 2.5)
Y = np.arange(0.0001, 1.0, 0.05)
#网格化数据
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(DdeltaDvol)
Z = vfunc(X,100.0,0.05,Y,0.2)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# DvannaDvol与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(0.01 ,200.0, 2.5)
Y = np.arange(0.0001, 1.0, 0.05)
#网格化数据
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(DvannaDvol)
Z = vfunc(X,100.0,0.07,Y,0.4)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# DdeltaDtime与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(50.0 , 150.0, 2.5)
Y = np.arange(0.0001, 0.2, 0.01)
#网格化数据
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(DdeltaDtime)
Z = vfunc(-1.0,X,100.0,0.05,Y,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# Gamma与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(50.0 , 150.0, 2.5)
Y = np.arange(0.015, 1.0, 0.02)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(Gamma)
Z = vfunc(X,100.0,0.05,Y,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(0.1 , 200.0, 2.5)
Y = np.arange(0.01, 4.0, 0.02)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(Gamma)
Z = vfunc(X,100.0,0.05,Y,0.8)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# DgammaDvol与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(50.0 , 150.0, 2.5)
Y = np.arange(0.015, 1.0, 0.02)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(DgammaDvol)
Z = vfunc(X,100.0,0.05,Y,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# DgammaDspot与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(50.0 , 150.0, 2.5)
Y = np.arange(0.015, 1.0, 0.02)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(DgammaDspot)
Z = vfunc(X,100.0,0.05,Y,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# DgammaDtime与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
X = np.arange(50.0 , 150.0, 2.5)
Y = np.arange(0.015, 0.25, 0.01)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(DgammaDtime)
Z = vfunc(X,100.0,0.05,Y,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# Vega与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
Y = np.arange(50.0 , 150.0, 2.5)
X = np.arange(0.015, 1.0, 0.02)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(Vega)
Z = vfunc(Y,100.0,0.05,X,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()
#Vega(55,60,0.105,0.75,0.3,q=0.0355)
figure = plt.figure()
ax = Axes3D(figure)
Y = np.arange(.1 , 250.0, 2.5)
X = np.arange(0.015, 2.0, 0.02)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(Vega)
Z = vfunc(Y,100.0,0.5,X,0.1,q=0.5)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# DvegaDvol与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
Y = np.arange(30.0 , 250.0, 2.5)
X = np.arange(0.03, 1.0, 0.02)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(DvegaDvol)
Z = vfunc(Y,100.0,0.05,X,0.2,0.05)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# DvommaDvol与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
Y = np.arange(50.0 , 150.0, 2.5)
X = np.arange(0.015, 1.0, 0.02)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(DvommaDvol)
Z = vfunc(Y,100.0,0.05,X,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()

# DvegaDtime与S、T的关系
figure = plt.figure()
ax = Axes3D(figure)
Y = np.arange(0.1 , 350.0, 2.5)
X = np.arange(0.015, 1.0, 0.02)
X, Y = np.meshgrid(X, Y)
vfunc=np.vectorize(DvommaDvol)
Z = vfunc(Y,100.0,0.05,X,0.3)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.show()