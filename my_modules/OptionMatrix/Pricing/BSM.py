import numpy as np
from scipy.stats import norm
from math import log,exp,sqrt

def cnd(x,method='default'): # 正态分布累计分布函数
    if method=='default':
        return norm.cdf(x)
def BS(call_condition,S,K,r,T,sigma):
    '''
    Black-Scholes(1973),without dividends, European option
    :param call_condition:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :return:
    '''
    assert (T>0.0) & (sigma>0.0)
    d1=(log(S/K)+(r+sigma**2.0*0.5)*T)/sigma*sqrt(T)
    d2=d1-sigma*sqrt(T)
    if call_condition>=0.0: # call
        return S*cnd(d1)-K*exp(-r*T)*cnd(d2)
    else: # put
        return K*exp(-r*T)*cnd(-d2)-S*cnd(-d1)

BS(-1.0,100.0,100.0,0.05,1.0,0.25)
T=1.0;sigma=0.25
T>0.0
True & (sigma>0.0)