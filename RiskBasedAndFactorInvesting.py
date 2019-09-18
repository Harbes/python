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
linear_constraint=LinearConstraint(matrix_LC,[1.0,r_p],[1.0,r_p]) #两个r_p换成一个range后，minimize报错
def func(W):
    return (W*SIGMA*W).sum()
r_p=0.04
W0=np.ones(4)/4
res=minimize(func,W0,method='trust-constr',constraints=linear_constraint,bounds=bounds)# 权重总是不满足约束条件，不知为何
# 使用cvxopt包
from cvxopt import matrix,solvers

P=matrix(SIGMA,tc='d') #二次项
q=matrix([0]*4,tc='d') #一次项
G=matrix(np.diag([-1]*4),tc='d')#不等式约束(全部改写成<=)
h=matrix([0]*4,tc='d')
A=matrix(np.vstack(([1]*4,r_mean)),tc='d')#等式约束
n=100
sigmas=np.zeros(n)
r_p=0.03+np.arange(n)*.0003
weights=np.zeros((4,n))
for i in range(n):
    b=matrix([1.0,r_p[i]],tc='d')
    sol=solvers.qp(P,q,G,h,A,b)
    sigmas[i]=np.sqrt(sol['primal objective']*2.0) # 注意，目标函数有1/2，故计算方差时要还原
    weights[:,i]=np.array(sol['x']).reshape(-1)
import matplotlib.pyplot as plt
plt.plot(sigmas*100,r_p*100);plt.show()
i_min_sigma=sigmas.argmin()
sigma_min=\
    sigmas.min()
i_max_SR=(r_p/sigmas).argmax()
SR_max=\
    (r_p/sigmas).max()
#===> MVP
weights[:,23]
r_p[23]
sigmas[23]
#===> MSRP
weights[:,33]
r_p[33]
sigmas[33]