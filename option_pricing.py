
### Pricing by Monte-Carlo
import numpy as np
import math
from numba import jit

M = 500  # nb of time steps of size dt
N = 500  # nb of stochastic  realization
L = 40  # nb of sampling point for S
K = 100  # strike price
T = 1;
sigma = 0.2;
r = 0.1;
er = np.exp(-r)
left = 0;
right = 130;  # the barrier
pi2 = 2 * math.pi
dt = T / M;
sdt = np.sqrt(dt);
eps = 1e-5


@jit
def stoch_barrier(x, m, Smin, Smax):
    if x <= Smin or x >= Smax:
        return -1000
    S = x
    gauss = np.random.standard_normal(m)
    for i in range(m):
        if S <= Smin or S >= Smax:
            return -1000
        S += S * sigma * gauss[i] * sdt + r * dt
    return S


@jit
def stoch(x, m):
    S = x
    gauss = np.random.standard_normal(m)
    for i in range(m):
        S += S * sigma * gauss[i] * sdt + r * dt
    return S


@jit
def cal_option(N, S0):
    value = 0
    for i in range(N):
        S = stoch_barrier(S0, M, left, right)
        y = 0
        if S >= 0:
            y = er * np.maximum(S - K,0)
        value += y
    return value / N
@jit
def cal_option_normal(N, S0):
    value = 0
    for i in range(N):
        S = stoch(S0, M)
        y = 0
        if S >= 0:
            y = er * np.maximum(S - K,0)
        value += y
    return value / N

import matplotlib.pyplot as plt

S0=np.linspace(2*K/L,2*K,L)
c=np.empty_like(S0)
for i in range(L):
    c[i]=cal_option_normal(5000,S0[i])
plt.plot(20*S0/K,c)
plt.show()



### Variance reduction by  control variate
import numpy as np
import math
from numba import jit
import matplotlib.pyplot as plt
from modules import euro_option as e_option
M=100  # nb of time steps of size dt
N=100000 #nb of stochastic realization
L=100  #nb of sampling point for S
K=100;sigma=0.2;r=0.1;er=np.exp(-r)
T=1
dt=1/M;sdt=np.sqrt(dt)

@jit
def cal_option_control(x):
    meanY=x/er;meanX=0;meanZ=0;barY=0;
    varX=0;varY=0;Exy=0
    for i in range(N):
        S1=x;S2=x
        gauss=np.random.standard_normal(M)
        for ii in range(M):
            y=sigma*gauss[ii]*sdt
            S1+=S1*(y+r*dt)
            S2+=S2*(-y+r*dt)
        Y=S1
        X=np.maximum(S1-K,0)*er
        Z = (X + np.maximum(S2 - K, 0) * er )/ 2
        meanX+=X;meanZ += Z; barY+=Y;
        Exy +=X*Y; varX +=X**2;varY+=Y**2
    meanX/=N;meanZ/=N;barY/=N;
    varX=varX/N-meanX**2
    varY=varY/N-2*meanY*barY+meanY**2
    Exy=Exy/N-meanX*barY
    b=Exy/varY;
    C=meanX-b*(barY-meanY)
    return meanX,meanZ,C
S0=np.arange(70,130.001,2)
c_mean=np.empty_like(S0)
c_antithetic=np.empty_like(S0)
c_control=np.empty_like(S0)
for i in range(S0.size):
    c_mean[i],c_antithetic[i],c_control[i]=cal_option_control(S0[i])
exact=e_option(S0,r=0.1).value_BSM
plt.plot(S0,c_mean-exact,label="Monte-Carlo")
plt.plot(S0,c_antithetic-exact,label="antithetic")
plt.plot(S0,c_control-exact,label="control")
plt.legend()
plt.show()