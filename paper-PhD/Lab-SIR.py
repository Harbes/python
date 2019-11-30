import numpy as np
import pandas as pd
P=60 # num of predictors
T= 500 # num of oberservations
rho=0.1 # average correlatioin
Sigma=np.zeros((P,P))
for i in range(P-1):
    Sigma+=np.eye(P,k=i+1)*(P-1-i)/P
Sigma=Sigma+Sigma.T+np.eye(P)
X=np.random.multivariate_normal(mean=np.zeros(P),cov=Sigma,size=T)
alpha=np.arange(P,0,-1)/P
Y=alpha@X.T
Y+=np.random.randn(T)*Y.std()
