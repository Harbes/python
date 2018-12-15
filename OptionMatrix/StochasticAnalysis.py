# Reference:
# <<Stochastic Analysis for Finance with Simulations>>, Geon Ho Choe.
#

import numpy as np
import matplotlib.pyplot as plt


# Chapter 3. the Lebesgue Integral
np.set_printoptions(precision=5)
def CantorSet():
    k=6
    M=10000
    Cantor=np.zeros(M)
    for i in range(M):
        for j in range(k):
            Cantor[i]=round((Cantor[i]+2.0*np.random.uniform())/3.0,5)
    return Cantor
#from decimal import *
#getcontext().prec = 5 # ln(BinsNumber)/ln(10)=< prec <= ln(3)/ln(10)*k

Cantor=CantorSet()
plt.hist(Cantor,3**10);plt.show()

np.log(3)/np.log(10)*7