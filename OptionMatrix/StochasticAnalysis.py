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
## Simulation 7.3 (Brownian motion with boundary condition)
import numpy as np
import matplotlib.pyplot as plt;
T=1
n=8
dt=T/2**n
time=np.arange(.0,T+dt*.1,dt)
M=2**n+1
w=np.zeros(M)
w[0]=0 # initial condition
w[-1]=.3 # condition on Brownian motion w at time T
for i in range(n):
    increment=2**(n-i)
    for j in range(0,2**n,increment):
        index1=j
        index2=j+increment
        t1=time[index1]
        t2=time[index2]
        w1=w[index1]
        w2=w[index2]
        ave=(w1+w2)*.5
        var=(t2-t1)*.25
        ind_mid=int((index1+index2)*.5)
        w[ind_mid]=np.random.normal(ave,np.sqrt(var))
for i in range(n+1):
    t_value=np.zeros(2**i+1)
    w_value=np.zeros(2**i+1)
    for k in range(2**i+1):
        t_value[k]=k*dt*2**(n-i)
        w_value[k]=w[k*2**(n-i)]
    plt.plot(t_value,w_value)
plt.show()
## Simulation 7.4 (LLN for Brownian motion)
N=1000
T=500
dt=T/N
time=np.arange(.0,T+dt*.1,dt)
M=50
W=np.zeros((M,N+1))
dW=np.sqrt(dt)*np.random.randn(M,N)
for i in range(N):
    W[:,i+1]=W[:,i]+dW[:,i]
X=np.zeros((M,N+1))
X[:,0]=1.0
for i in range(1,N+1):
    X[:,i]=W[:,i]/(i*dt)
for j in range(M):
    plt.plot(time,X[j,:])
plt.show()

# Chapter 10. the Ito Integral
import numpy as np
import matplotlib.pyplot as plt
## Simulation 10.1 (Ito integral)
T=3.0
N=1000
dt=T/N
t=np.arange(.0,T+dt*.1,dt)
dW=np.sqrt(dt)*np.random.randn(N)
W=np.zeros(N+1)
integral=np.zeros(N+1)
exact=np.zeros(N+1)
for i in range(N):
    W[i+1]=W[i]+dW[i]
    integral[i+1]=integral[i]+W[i]*dW[i]
    exact[i+1]=W[i+1]*W[i+1]*.5-(i+1)*dt*.5
plt.plot(t,exact,'k-')
plt.plot(t,integral,'r.')
plt.show()
## Simulation 10.2 (Ito integral)
T=3
M=10
N=500
dt=T/N
t=np.arange(.0,T+dt*.1,dt)
dW=np.sqrt(dt)*np.random.randn(M,N)
W=np.zeros((M,N+1))
W[:,1:]+=dW.cumsum(axis=1)
integral=np.zeros((M,N+1))
integral[:,1:]+=(W[:,:-1]*dW).cumsum(axis=1) # 此处用的是W[:,:-1]，而不是W[:,1:]，结果差距很大，因此要注意时间的匹配
for i in range(M):
    plt.plot(t,integral[i])
plt.plot(t,-.5*t)
plt.show()
# Simulation 10.3 (convergence)
T=3
M=100000
N=500
dt=T/N
t=np.arange(.0,T+dt*.1,dt)
dW=np.sqrt(dt)*np.random.randn(M,N)
W=np.zeros((M,N+1))
W[:,1:]+=dW.cumsum(axis=1)
integral=np.zeros((M,N+1))
integral[:,1:]+=(W[:,:-1]*dW).cumsum(axis=1) # 此处用的是W[:,:-1]，而不是W[:,1:]，结果差距很大，因此要注意时间的匹配
exact=W[:,-1]**2*.5-T*.5
error=integral[:,-1]-exact
print(np.mean(error**2))