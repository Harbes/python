# chapter 2=========================================
## example 2.1: stationary AR(1) process
from modules import simulation as simu
import numpy as np
import matplotlib.pyplot as plt
alpha=1
rho=0.95
sig=2
mu=alpha/(1-rho)
sig2=sig**2/(1-rho**2)
gamma1=sig2*rho
M=5000
y=simu(phi=rho,const=alpha,dist_sigma=sig,M=M+50).AR[50:]
mu_=np.mean(y)
sig2_=np.var(y)
gamma1_=np.mean((y-mu_)[1:]*(y-mu_)[:-1])
print(mu,mu_)
print(sig2,sig2_)
print(gamma1,gamma1_)
plt.plot(y)
plt.show()

## example 2.2: Nonlinear time series
import numpy as np
import matplotlib.pyplot as plt
M=5000
u=np.random.normal(0,2,M)
y=u[1:]*u[:-1]
print('mean:',np.mean(y))
print('var & cov',np.cov(y[1:],y[:-1])[0])
print('higher moments:',np.cov((y**2)[1:],(y**2)[:-1])[0,1])
plt.plot(y,label='y')
plt.plot(u,label='u')
plt.legend()
plt.show()

## example 2.3: autoregressive conditional heteroskedasticity
import numpy as np
import matplotlib.pyplot as plt
M=2000
z=np.random.standard_normal(M)
y=np.zeros(M)
alpha0=10
alpha1=0.9
for i in range(1,M):
    y[i]=z[i]*np.sqrt(alpha0+alpha1*(y[i-1])**2)
print('mean:',np.mean(y))
print('var & cov',np.cov(y[1:],y[:-1])[0])
print('higher moments:',np.mean(y**2))
plt.plot(y,label='y')
plt.legend()
plt.show()

## example 2.4: bilinear time series
import numpy as np
import matplotlib.pyplot as plt
M=500
sig=2
u=np.random.normal(0,sig,M)
delta=0.95
y=u[2:]+delta*u[1:-1]*u[:-2]
print('mean:',np.mean(y))
print('2nd moment:',np.mean(y**2))
print('co-2nd moment',np.mean(y[1:]*y[:-1]))
plt.plot(y)
plt.show()


## example 2.5: exponential distribution
import numpy as np
import matplotlib.pyplot as plt
M=5000
m=np.empty(M)
mu=5
r=0.2
for i in range(1,M+1):
    m[i-1]=np.mean(np.random.exponential(mu,i))
plt.plot(m)
plt.plot(np.arange(1,M+1),mu*np.ones(M),'k-')
plt.plot(np.arange(1,M+1),(mu-r)*np.ones(M),'r--')
plt.plot(np.arange(1,M+1),(mu+r)*np.ones(M),'r--')
plt.show()



































