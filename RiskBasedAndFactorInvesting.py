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
# calculate MVP without nonnegative constraint
weights_mvp=SIGMA_inv.sum(axis=0)/SIGMA_inv.sum() # IG的权重为负
r_mvp=(weights_mvp*r_mean).sum()
std_mvp=np.sqrt(1/SIGMA_inv.sum())
# constrained minimization# 不知为何，报错
from scipy.optimize import minimize,Bounds,LinearConstraint
bounds=Bounds(np.zeros(4),np.ones(4))
matrix_LC=np.array([np.ones(4),r_mean])
p_r=0.04
linear_constraint=LinearConstraint(matrix_LC, [1.0, p_r], [1.0, p_r]) #两个r_p换成一个range后，minimize报错
def func(W):
    return SIGMA.dot(W).dot(W)
W0=np.ones(4)/4
res=minimize(func,W0,method='trust-constr',constraints=linear_constraint,bounds=bounds)# 似乎要设置jac矩阵
# 使用cvxopt包，计算权重非负约束下的efficient frontier
from cvxopt import matrix,solvers
P=matrix(SIGMA,tc='d') #二次项
q=matrix([0]*4,tc='d') #一次项
G=matrix(np.diag([-1]*4),tc='d')#不等式约束(全部改写成<=)
h=matrix([0]*4,tc='d')
A=matrix(np.vstack(([1]*4,r_mean)),tc='d')#等式约束
n=100
p_std=np.zeros(n)
p_r= 0.03 + np.arange(n) * .0003
p_weights=np.zeros((4, n))
for i in range(n):
    b=matrix([1.0, p_r[i]], tc='d')
    sol=solvers.qp(P,q,G,h,A,b)
    p_std[i]=np.sqrt(sol['primal objective'] * 2.0) # 注意，目标函数有1/2，故计算方差时要还原
    p_weights[:, i]=np.array(sol['x']).reshape(-1)
import matplotlib.pyplot as plt
plt.plot(p_std * 100, p_r * 100);plt.show() # 需要绘制详细一些


SR_max=\
    (p_r / p_std).max()
#===> MVP
i_std_min=p_std.argmin()
p_weights[:, i_std_min]
p_r[i_std_min]
p_std[i_std_min]
SR_std_min=p_r[i_std_min]/p_std[i_std_min];SR_std_min
#===> MSRP
i_SR_max=(p_r / p_std).argmax()
p_weights[:, i_SR_max]
p_r[i_SR_max]
p_std[i_SR_max]
(p_r/p_std)[i_SR_max]
# 1/N or equally-weighted
p_ew_weight=np.ones(4)/4
p_ew_r=r_mean.mean();p_ew_r
p_ew_std=np.sqrt(SIGMA.dot(p_weight).dot(p_weight));p_ew_std
p_ew_SR=p_ew_r/p_ew_std;p_ew_SR