import numpy as np

########## Ch1
###MSRP###
#input
AssetNames=['equity','Tsies','IG','HY']
r_mean=np.array([7.54,2.97,4.17,6.62])*.01
r_std=np.array([15.11,4.16,6.04,10.54])*.01
corr=np.array([[1.0,-.3,.35,.74],
               [-.3,1.0,.44,-.24],
               [.35,.44,1.0,.63],
               [.74,-.24,.63,1.0]])
SIGMA=(corr*r_std).T*r_std
SIGMA_inv=np.linalg.pinv(SIGMA)
# calculate MVP
weights_mvp=SIGMA_inv.sum(axis=0)/SIGMA_inv.sum() # IG的权重为负
r_mvp=(weights_mvp*r_mean).sum()
std_mvp=np.sqrt(1/SIGMA_inv.sum())
# constrained minimization
from scipy.optimize import minimize,Bounds,LinearConstraint
bounds=Bounds(np.zeros(4),np.ones(4))
matrix_LC=np.array([np.ones(4),r_mean])
linear_constraint=LinearConstraint(matrix_LC,[1.0,0.03],[1.0,0.05])
def func(W):
    return (W*SIGMA*W).sum()
r_p=0.04
W0=np.ones(4)/4
res=minimize(func,W0,method='trust-constr',constraints=linear_constraint,bounds=bounds)