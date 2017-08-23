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





