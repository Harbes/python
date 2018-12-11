#import numpy as np
from scipy.stats import norm
from math import log,exp,sqrt

def cnd(x,method='default'): # cumulative normal distribution function
    '''
    <The Complete Guide to Option Pricing Formulas>, Chapter 13
    '''
    if method=='default':
        return norm.cdf(x)
def BSM(call_condition,S,K,r,T,sigma,q=0.0):
    '''
    European option, Black-Scholes(1973), allowing for a continuous dividend yield
    '''
    assert (T>0.0) & (sigma>0.0)
    d1=(log(S/K)+(r-q+sigma**2.0*0.5)*T)/sigma/sqrt(T)
    d2=d1-sigma*sqrt(T)
    if call_condition>=0.0: # call
        return S*exp(-q*T)*cnd(d1)-K*exp(-r*T)*cnd(d2)
    else: # put
        return K*exp(-r*T)*cnd(-d2)-S*exp(-q*T)*cnd(-d1)
def GeneralizedBSM(call_condition,S,K,r,T,sigma,q=0.0):
    '''
    European option, BSM model be "generalized" by incorporating a cost-of-carry rate b;
    c=S*exp((b-r)T)*N(d1)-K*exp(-rT)*N(d2)
    p=K*exp(-rT)*N(-d2)-S*exp((b-r)T)*N(-d1)
    d1=(ln(S/K)+(b+sigma**2/2)T)/(sigma*sart(T))
    d2=d1-sigma*sqrt(T)
    ====================
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    '''
    assert (T>0.0) & (sigma>0.0)
    b=r
    d1=(log(S/K)+(b+sigma**2.0/2.0)*T)/sigma/sqrt(T)
    d2=d1-sigma*sqrt(T)
    if call_condition>=0.0: # call
        return S*exp((b-r)*T)*cnd(d1)-K*exp(-r*T)*cnd(d2)
    else: # put
        return K*exp(-r*T)*cnd(-d2)-S*exp((b-r)*T)*cnd(-d1)

## TODO 试验区，待删除
BSM(-1.0,100.0,100.0,0.05,1.5,0.25)
GeneralizedBSM(-1.0,100.0,100.0,0.05,1.5,0.25)
