from scipy.optimize import minimize,Bounds,LinearConstraint,SR1
import numpy as np

AssetNames=['equity','Tsies','IG','HY']
r_mean=np.array([7.54,2.97,4.17,6.62])*.01
r_std=np.array([15.11,4.16,6.04,10.54])*.01
corr=np.array([[1.0,-.3,.35,.74],
               [-.3,1.0,.44,-.24],
               [.35,.44,1.0,.63],
               [.74,-.24,.63,1.0]])
SIGMA=(corr*r_std).T*r_std
SIGMA_inv=np.linalg.pinv(SIGMA)

bounds=Bounds(np.zeros(4),np.ones(4))
matrix_LC=np.ones(4)
p_r=0.04
linear_constraint=LinearConstraint(matrix_LC, [1.0], [1.0])
def func(W):
    return -W.dot(r_std)/np.sqrt(W.dot(SIGMA).dot(W))
def func_der(W):
    return (SIGMA.dot(W)-W.dot(SIGMA).dot(W)*r_std)/(W.dot(SIGMA).dot(W))**1.5
W0=np.ones(4)/4
res=minimize(func,W0,method='trust-constr',jac='2-point',hess=SR1(),constraints=linear_constraint,bounds=bounds)# 似乎要设置jac矩阵
print('Weights:',res.x,'\np_r:',res.x.dot(r_mean),'\np_std:',np.sqrt(res.x.dot(SIGMA).dot(res.x)),'\nSR:',res.x.dot(r_mean)/np.sqrt(res.x.dot(SIGMA).dot(res.x)))
