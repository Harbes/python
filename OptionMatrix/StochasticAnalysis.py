# Reference:
# <<Stochastic Analysis for Finance with Simulations>>, Geon Ho Choe.
#

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Chapter 3. the Lebesgue Integral
## Cantor set
np.set_printoptions(precision=5)
def CantorSet():
    k=6
    M=10000
    Cantor=np.zeros(M)
    for i in range(M):
        for j in range(k):
            Cantor[i]=round((Cantor[i]+2.0*np.random.uniform())/3.0,5) #理论上，ln(BinsNumber)/ln(10)=< prec <= ln(3)/ln(10)*k
    return Cantor
Cantor=CantorSet()
plt.hist(Cantor,3**10);plt.show() # 结果似乎不是很理想


# Chapter 4. Basic Probability Theory
## Simulation 4.2 (equally speced points in the unit interval)
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

