# Reference:
# <<Stochastic Analysis for Finance with Simulations>>, Geon Ho Choe.
#

import numpy as np
import matplotlib.pyplot as plt;
from scipy.stats import norm

# Chapter 3. the Lebesgue Integral
## Cantor set
def CantorSet():
    k=20
    M=10000
    Cantor=np.zeros(M)
    for i in range(M):
        for j in range(k):
            Cantor[i]=(Cantor[i]+2.0*np.random.randint(0,2))/3.0 #理论上，ln(BinsNumber)/ln(10)=< prec <= ln(3)/ln(10)*k
    return Cantor
Cantor=CantorSet()
plt.hist(Cantor,bins=3**7);plt.show() # 结果似乎不是很理想

# Chapter 4. Basic Probability Theory
## Simulation 4.2 (equally speced points in the unit interval)-；
N=100
X=np.random.randn(N)
Y=np.sort(X)
U=np.arange(1/(N+1),1.0,1/(N+1))
V=norm.cdf(Y)
x=np.arange(0,1.0,0.01)
plt.plot(x,x);plt.plot(V,U,'o');plt.show()
## Simulation 4.3 (Quantile-Quantile plot)
N=100
X=np.random.randn(N)
Y=np.sort(X)
U=np.arange(1/(N+1),1.0,1/(N+1))
V=norm.ppf(U)
x=np.arange(-4,4,0.02)
plt.plot(x,x);plt.plot(Y,V,'o');plt.show()
## Simulation 4.4 (correlated normal varibles)
N=5000
z1=np.random.randn(N)
z2=np.random.randn(N)
print(np.corrcoef(z1,z2))
rho=-0.9
x1=z1
x2=rho*z1+np.sqrt(1-rho*rho)*z2
print(np.corrcoef(x1,x2))
plt.figure()
plt.subplot(1,2,1)
plt.plot(z1,z2,'.')
plt.subplot(1,2,2)
plt.plot(x1,x2,'.')
plt.show()
## Simulation 4.5 (Noncentral Chi2-distribution)
from scipy.stats import ncx2
x=np.arange(.001,10.01,0.1)
for k in [2,5]:
    for lam in range(1,4):
        y=ncx2.pdf(x,k,lam)
        plt.plot(x,y)
plt.show()

# Chapter 7. Brownian motion
## Simulation7.1 (sample paths of Brownian motion)
n=500
t=20
dt=t/n
time=np.arange(0,t+dt*0.1,dt)
m=50 # number of sample paths
w=np.zeros((m,n+1))
dw=np.sqrt(dt)*np.random.randn(m,n)
for i in range(n):
    w[:,i+1]=w[:,i]+dw[:,i]
for j in range(m):
    plt.plot(time,w[j,:])
plt.show()
## Simulation 7.2 (average of geometric Brownian motion)
mu=.25
sigma=.3
n=200
t=5
tt=np.arange(.0,t+0.001,.01)
dt=t/n
time=np.arange(.0,t+dt*0.1,dt)
num=1000 # number of sample paths
w=np.zeros((num,n+1))
dw=np.zeros((num,n))
s=np.zeros((num,n+1))
s0=1
s[:,0]=s0
for i in range(num):
    dw[i,:]=np.sqrt(dt)*np.random.randn(n)
for i in range(num):
    for j in range(n):
        s[i,j+1]=s[i,j]+mu*s[i,j]*dt+sigma*s[i,j]*dw[i,j]
ave=np.zeros(n+1)
for j in range(n+1):
    ave[j]=s[:,j].mean()
plt.plot(time,ave);#plt.show()
plt.plot(tt,s0*np.exp(mu*tt));plt.show()
