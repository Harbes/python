import numpy as np
from scipy.optimize import minimize,Bounds,LinearConstraint,NonlinearConstraint,BFGS,SR1
import matplotlib.pyplot as plt

########## Ch1
#input
AssetNames=['equity','Tsies','IG','HY']
r_mean=np.array([7.54,2.97,4.17,6.62])*.01
r_std=np.array([15.11,4.16,6.04,10.54])*.01
corr=np.array([[1.0,-.3,.35,.74],
               [-.3,1.0,.44,-.24],
               [.35,.44,1.0,.63],
               [.74,-.24,.63,1.0]])
SIGMA=(corr*r_std).T*r_std

### MVP without nonnegative constraint
SIGMA_inv=np.linalg.pinv(SIGMA)
weights_mvp=SIGMA_inv.sum(axis=0)/SIGMA_inv.sum() # IG的权重为负
r_mvp=(weights_mvp*r_mean).sum()
std_mvp=np.sqrt(1/SIGMA_inv.sum())
### MVP with non-negative constraint, given expected portfolio return
bounds=Bounds(np.zeros(4),np.ones(4))
matrix_LC=np.array([np.ones(4),r_mean])
p_r=0.04
linear_constraint=LinearConstraint(matrix_LC, [1.0, p_r], [1.0, p_r])
def func(W):
    return SIGMA.dot(W).dot(W)
def func_der(W):
    return SIGMA.dot(W)*2.0
def func_hess(W):
    return SIGMA*2.0
W0=np.ones(4)/4
res=minimize(func,W0,method='trust-constr',jac=func_der,hess=func_hess,constraints=linear_constraint,bounds=bounds)# 似乎要设置jac矩阵
print('Weights:',res.x,'\nMinimum std:',np.sqrt(res.fun))
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

plt.plot(p_std * 100, p_r * 100);plt.show() # 最小方差前沿
## MVP，方法1
i_std_min=p_std.argmin()
p_weights[:, i_std_min]
p_r[i_std_min]
p_std[i_std_min]
SR_std_min=p_r[i_std_min]/p_std[i_std_min];SR_std_min
## MSRP，方法1
i_SR_max=(p_r / p_std).argmax()
p_weights[:, i_SR_max]
p_r[i_SR_max]
p_std[i_SR_max]
(p_r/p_std)[i_SR_max]

### MSRP with non-negative constraint，方法2
W0=np.ones(4)/4
bounds=Bounds(np.zeros(4),[np.inf]*4)
matrix_LC=np.ones(4)
linear_constraint=LinearConstraint(matrix_LC, 1.0,1.0)
def func(W):
    return -W.dot(r_mean)/(W.dot(SIGMA).dot(W))
def func_der(W):
    return (W.dot(r_mean)*SIGMA.dot(W)-W.dot(SIGMA).dot(W)*r_mean)/(W.dot(SIGMA).dot(W))*1.5
res=minimize(func,W0,method='trust-constr',jac=func_der,hess=SR1(),bounds=bounds,constraints=linear_constraint)
print('Weights:',res.x,'\np_r:',res.x.dot(r_mean),'\np_std:',np.sqrt(res.x.dot(SIGMA).dot(res.x)),'\nSR:',res.x.dot(r_mean)/np.sqrt(res.x.dot(SIGMA).dot(res.x)))

### MVP with non-negative constraint, 方法2
W0=np.ones(4)/4
bounds=Bounds(np.zeros(4),[np.inf]*4)
matrix_LC=np.ones(4)
linear_constraint=LinearConstraint(matrix_LC, 1.0,1.0)
def func(W):
    return W.dot(SIGMA).dot(W)
def func_der(W):
    return W.dot(SIGMA)*2.0
res=minimize(func,W0,method='trust-constr',jac=func_der,hess='2-point',bounds=bounds,constraints=linear_constraint,tol=1e-8)
print('Weights:',res.x,'\np_r:',res.x.dot(r_mean),'\np_std:',np.sqrt(res.x.dot(SIGMA).dot(res.x)),'\nSR:',res.x.dot(r_mean)/np.sqrt(res.x.dot(SIGMA).dot(res.x)))

### 1/N or equally-weighted
p_ew_weight=np.ones(4)/4
p_ew_r=r_mean.mean();p_ew_r
p_ew_std=np.sqrt(SIGMA.dot(p_ew_weight).dot(p_ew_weight));p_ew_std
p_ew_SR=p_ew_r/p_ew_std;p_ew_SR

### MDP with non-negative constraint
bounds=Bounds(np.zeros(4),np.ones(4))
matrix_LC=np.ones(4)
p_r=0.04
linear_constraint=LinearConstraint(matrix_LC, [1.0], [1.0])
def func(W):
    return -W.dot(r_std)/np.sqrt(W.dot(SIGMA).dot(W))
def func_der(W):
    return (W.dot(r_std)*SIGMA.dot(W)-W.dot(SIGMA).dot(W)*r_std)/(W.dot(SIGMA).dot(W))**1.5
W0=np.ones(4)/4
res=minimize(func,W0,method='trust-constr',jac=func_der,hess=SR1(),constraints=linear_constraint,bounds=bounds)# hess=SR1()是近似计算hess矩阵，另外，可以使用jac='2-point'近似计算
print('Weights:',res.x,'\np_r:',res.x.dot(r_mean),'\np_std:',np.sqrt(res.x.dot(SIGMA).dot(res.x)),'\nSR:',res.x.dot(r_mean)/np.sqrt(res.x.dot(SIGMA).dot(res.x)))

### risk budgeting，方法1【不需要leverage=1的约束，只需要最后scale一下即可】
budgeting=np.arange(1,5)/np.arange(1,5).sum()
c=budgeting.dot(np.log(budgeting))*1.001
def func(W):
    return SIGMA.dot(W).dot(W)
def func_der(W):
    return SIGMA.dot(W)*2.0
def func_hess(W):
    return SIGMA*2.0
def cons_f(W):
    return c-budgeting.dot(np.log(W))
def cons_jac(W):
    return -budgeting/W
from scipy.optimize import BFGS
nonlinear_constraint=NonlinearConstraint(cons_f,-np.inf,0.0,jac=cons_jac,hess=BFGS())
res=minimize(func,budgeting,method='trust-constr',jac=func_der,hess=func_hess,constraints=nonlinear_constraint)
print('Weights:',res.x/res.x.sum(),'\np_r:',res.x.dot(r_mean)/res.x.sum(),'\np_std:',np.sqrt(res.x.dot(SIGMA).dot(res.x))/res.x.sum(),'\nSR:',res.x.dot(r_mean)/np.sqrt(res.x.dot(SIGMA).dot(res.x)))
W=res.x/res.x.sum();SIGMA.dot(W)*W
### risk budgeting，方法2【不需要额外w>0的假设，因为ln(w)就隐含了w>0】
budgeting=np.arange(2,6);budgeting=budgeting/budgeting.sum()
matrix_LC=np.ones(4)
linear_constraint=LinearConstraint(matrix_LC, [1.0], [1.0])
def func(W):
    CR_b=SIGMA.dot(W)*W/budgeting
    v=0.0
    for i in range(1,len(CR_b)):
        v+=((CR_b[:-i]-CR_b[i:])**2.0).sum()
    return v
res=minimize(func,budgeting,method='trust-constr',jac='2-point',hess=BFGS(),constraints=linear_constraint)
print('Weights:',res.x,'\np_r:',res.x.dot(r_mean),'\np_std:',np.sqrt(res.x.dot(SIGMA).dot(res.x)),'\nSR:',res.x.dot(r_mean)/np.sqrt(res.x.dot(SIGMA).dot(res.x)))
W=res.x;SIGMA.dot(W)*W/(SIGMA.dot(W)*W).sum()

### IVP
p_ivp_weight=1.0/r_std/(1.0/r_std).sum()
p_ivp_r=p_ivp_weight.dot(r_mean);p_ivp_r
p_ivp_std=np.sqrt(p_ivp_weight.dot(SIGMA).dot(p_ivp_weight));p_ivp_std
p_ivp_SR=p_ivp_r/p_ivp_std;p_ivp_SR
