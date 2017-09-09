
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
        S += S *math.exp(sigma * gauss[i] * sdt + r * dt)
    return S


@jit
def stoch(x, m):
    S = x
    gauss = np.random.standard_normal(m)
    for i in range(m):
        S += S *math.exp(sigma * gauss[i] * sdt + r * dt)
    return S


@jit
def cal_option(N, S0):
    value = 0
    for i in range(N):
        S = stoch_barrier(S0, M, left, right)
        y = 0
        if S >= 0:
            y = er * np.maximum(K-S,0)
        value += y
    return value / N
@jit
def cal_option_normal(N, S0):
    value = 0
    for i in range(N):
        S = stoch(S0, M)
        y = 0
        if S >= 0:
            y = er * np.maximum(K-S,0)
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


# European call by binomial tree
import numpy as np
from numba import jit
import math
import matplotlib.pyplot as plt
from modules import euro_option as e_option
from scipy.special import comb
M=81  # nb of time steps of size dt
N=100000 #nb of stochastic realization
L=100  #nb of sampling point for S
K=100;sigma=0.2;r=0.1;er=np.exp(-r)
T=1
dt=T/(M-1);sdt=np.sqrt(dt)

@jit
def option_binomial(S0):
    disc=math.exp(-r*dt)
    u=(1+math.sqrt(math.exp(sigma**2*dt)-1))/disc
    d=(1-math.sqrt(math.exp(sigma**2*dt)-1))/disc
    p=0.5
    S=np.empty(M+1);S[0]=S0
    for m in range(1,M+1):
        for n in range(m,0,-1):
            S[n]=u*S[n-1]
        S[0]=d*S[0]
    C=np.empty(M+1)
    for n in range(M+1):
        C[n]=S[n]-K if S[n]>K else 0
    for m in range(M,0,-1):
        for n in range(m):
            C[n]=(p*C[n+1]+(1-p)*C[n])*disc
    return C[0]

@jit
def option_binomial_comb(St):#,K,r,T,sigma,M,otype='call',American=False):
    K = 100.0; r = 0.05; T = 1.0; sigma = 0.2; M = 80; otype = 'call'; American = False
    dt = T / M;
    sdt = math.sqrt(dt)
    disc = math.exp(-r * dt)
    u=(1+math.sqrt(math.exp(sigma**2*dt)-1))/disc
    d=(1-math.sqrt(math.exp(sigma**2*dt)-1))/disc
    p=0.5
    #alpha = (math.exp(-r * dt) + math.exp((r + sigma ** 2) * dt)) / 2.0
    #u = alpha + math.sqrt(alpha ** 2 - 1)
    #d = 1.0 / u
    #p = (1.0 / disc - d) / (u - d)
    S = np.empty(M+1);
    for i in range(M+1):
        S[i]=St*u**(i)*d**(M-i)
    C = np.zeros(M+1)
    for n in range(M+1):
        C[n] = S[n] - K if S[n] > K else 0
    #C=np.where(S-K>0,S-K,C)
    for m in range(M , 0, -1):
        for n in range(m):
            C[n] = (p * C[n + 1] + (1 - p) * C[n]) * disc
    return C[0]
%timeit option_binomial_comb(100.0)#,100.0,0.05,1.0,0.2,80)
S0=np.arange(20,130)
c_binomial=np.empty(S0.size)
%timeit option_binomial(S0[0])
%timeit for i in range(S0.size): c_binomial[i]=option_binomial(S0[i])
opt_bin=np.vectorize(option_binomial)
%timeit opt_bin(S0=S0)
%timeit exact=e_option(St=S0,r=0.1).value_BSM()
plt.plot(S0,c_binomial,label="binomial")
plt.plot(S0,exact,label='BSM')
plt.legend()
plt.show()

%timeit e_option().value_BSM()
%timeit option_binomial(100)
exact=e_option(r=0.1).value_BSM


# binary tree for an American put
import numpy as np
from numba import jit
import math
import matplotlib.pyplot as plt
from modules import ame_option as a_option,euro_option as e_option
M=50  # nb of time steps of size dt
N=100000 #nb of stochastic realization
L=100  #nb of sampling point for S
K=100;sigma=0.2;r=0.1;
T=1
dt=T/M;sdt=np.sqrt(dt)

@jit
def option_binomial(S0):
    disc=math.exp(-r*dt)
    #u=(1+math.sqrt(math.exp(sigma**2*dt)-1))/disc
    #d=(1-math.sqrt(math.exp(sigma**2*dt)-1))/disc
    #p=0.5
    alpha=(math.exp(-r*dt)+math.exp((r+sigma**2)*dt))/2
    u=alpha+math.sqrt(alpha**1-1)
    d=1/u
    p=(1/disc-d)/(u-d)
    S=np.empty(M);S[0]=S0
    um=np.empty(M);um[0]=1
    du=np.empty(M);du[0]=1
    for m in range(1,M):
        for n in range(m,0,-1):
            S[n]=u*S[n-1]
        S[0]=d*S[0]
        um[m]=u*um[m-1]
        du[m]=du[m-1]*d/u
    P=np.zeros(M)
    for n in range(M):
        P[n]=K-S[n] if K>S[n] else 0
    for m in range(M-1,0,-1):
        for n in range(m):
            P[n]=(p*P[n+1]+(1-p)*P[n])*disc
            gain=K-S0*um[m]*du[n]
            if gain>P[n]:
                P[n]=gain
    return P[0]

S0=np.arange(10,130)
p_binomial=np.empty(S0.size)
exact=np.empty(S0.size)
for i in range(S0.size):
    p_binomial[i]=option_binomial(S0[i])
    exact[i]=a_option(St=S0[i],r=0.1,otype='put').value_binomial_0(M=100)
exact1=e_option(St=S0,r=0.1,otype='put').value_BSM
plt.plot(S0,p_binomial,label="binomial")
plt.plot(S0,exact,label='for comparison')
plt.plot(S0,exact1,label="Euro put")
plt.legend()
plt.show()

%timeit e_option().value_BSM
%timeit option_binomial(100)
exact=e_option(r=0.1).value_BSM


import numpy as np
import math
from numba import jit
from modules import Ame_option_binomial as a_bin_fun
import matplotlib.pyplot as plt

@jit
def Ame_option_trinomial(St,K,r,sigma,T,M,otype='call',method='CRR'):
    mu = r - sigma * sigma / 2.0
    lam=math.sqrt(3.0)
    dt = T/ M
    disc = math.exp(-r * dt)
    v = 0 if method=='CRR' else mu
    k1 = math.exp((r - v) * dt)
    k2 = math.exp((2*(r-v)+sigma*sigma)*dt)
    u,m=math.exp(lam*sigma*math.sqrt(dt)),math.exp(v*dt)
    d=1.0/u
    pu=(k2-(d+1)*k1+d)/(u-d)/(u-1.0)
    pd=(k2-(u+1)*k1+u)/(u-d)/(1.0-d)
    pm=1.0-pu-pd
    um=(u*m)**np.arange(M+1)
    dm=(1/u)**np.arange(2*M+1)
    S = np.empty(2*M + 1);S[0]=St*(d*m)**M
    for i in range(1,2*M+1):
        S[i]=S[i-1]*u
    o_value=1.0 if otype is 'call' else -1.0
    payoff = np.zeros(2*M + 1)
    for n in range(2*M + 1):
        payoff[n] = (S[n]-K)*o_value if (S[n]-K)*o_value >0 else 0.0
    for m in range(M, 1, -1):
        for n in range(2*m-1):
            tmp = (pu* payoff[n + 2] +pm*payoff[n+1]+ pd* payoff[n]) * disc
            payoff[n] = (St *um[m-1] * dm[2*(m-1)-n]-K)*o_value if (St * um[m-1] * dm[2*(m-1)-n]-K)*o_value- tmp>0 else tmp
    return (pu* payoff[2] +pm*payoff[1]+ pd* payoff[0]) * disc

a_bin_vecfun=np.vectorize(a_bin_fun)
a_tri_vecfun=np.vectorize(Ame_option_trinomial)
S=np.arange(10,130)
p_bin=a_bin_vecfun(S,100.0,.05,.2,1.0,80,otype='put')
p_tri=a_tri_vecfun(S,100.0,.05,.2,1.0,30,otype='put',method='CRR')
plt.plot(S,p_bin,label='Bin')
plt.plot(S,p_tri,label='Tri')
plt.legend()
plt.show()




import numpy as np
import math
from numba import jit
from modules import Ame_option_binomial as a_bin_fun
import matplotlib.pyplot as plt

@jit
def Ame_option_ExDiff(St,K,r,sigma,T,M,N,otype='call'):
    mu = r - sigma * sigma / 2.0
    dt = T/ N
    dlnS=sigma*math.sqrt(1.5*dt)
    pu=sigma*sigma*dt/(2.0*dx*dx)+mu*dt/(2.0*dx)
    pd=sigma*sigma*dt/(2.0*dx*dx)-mu*dt/(2.0*dx)
    pm=1.0-pu-pd
    S = St*np.exp(np.arange(-M ,M+ 1)*dx)
    o_value=1.0 if otype is 'call' else -1.0
    payoff = np.where((S-K)*o_value >0,(S-K)*o_value,0)
    for m in range(N, 0, -1):
        for n in range(1,2*m-1):
            tmp = pu* payoff[n + 2] +pm*payoff[n+1]+ pd* payoff[n]
            payoff[n] = (S[n]-K)*o_value if (S[n]-K)*o_value- tmp>0 else tmp

    return (pu* payoff[2] +pm*payoff[1]+ pd* payoff[0]) * disc

a_bin_vecfun=np.vectorize(a_bin_fun)
a_tri_vecfun=np.vectorize(Ame_option_trinomial)
S=np.arange(10,130)
p_bin=a_bin_vecfun(S,100.0,.05,.2,1.0,80,otype='put')
p_tri=a_tri_vecfun(S,100.0,.05,.2,1.0,30,otype='put',method='CRR')
plt.plot(S,p_bin,label='Bin')
plt.plot(S,p_tri,label='Tri')
plt.legend()
plt.show()
















