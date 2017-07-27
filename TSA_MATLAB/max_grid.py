#=====================================================
#
#  Program to find the MLEs using grid search methods
#
#=====================================================
import numpy as np
from math import sqrt
np.random.seed(123)

# simulate the model
x=np.array([1,2,4,5,8])
beta=1.0
sig2=4.0
t=len(x)
y=beta*x+sqrt(sig2)*np.random.standard_normal(t)

#grid search on gradient sig2=4










