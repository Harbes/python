from math import exp,sqrt,pi
from scipy.stats import norm

def HartAlgo(x):
    y=abs(x)
    k=exp(-y**2.0*0.5)
    if y>37.0:
        p=0.0
    elif y>=7.07106781186547:
        p=k/(y+1.0/(y+2.0/(y+3.0/(y+4.0/(y+0.65)))))/2.506628274631
    else:
        A=(((((0.0352624965998911*y+0.700383064443688)*y+6.37396220353165)*y+33.912866078383)*y+112.079291497871)*y+
           221.213596169931)*y+220.206867912376
        B=((((((0.0883883476483184*y+1.75566716318264)*y+16.064177579207)*y+86.7807322029461)*y+296.564248779674)*y+
           637.333633378831)*y+793.826512519948)*y+440.413735824752
        p=k*A/B
    if x>0.0:
        p=1.0-p
    return p

def PolyApprox4(x):
    '''
    the following approximation of N(x) produces values to within 4-decimal-place accuracy
    :param x:
    :return:
    '''
    y=abs(x)
    k=1.0/(1.0+0.33267*y)
    p=exp(-y**2.0*0.5)/sqrt(2.0*pi)*(0.4361836*k-0.1201676*k**2.0+0.9372980*k**3.0)
    if x>=0.0:
        return 1.0-p
    else:
        return p

def PolyApprox6(x):
    '''
    the following approximation of N(x), described by Abromowitz and Stegun(1974), produces values to within
    6 decimal places of the true value.
    :param x:
    :return:
    '''
    y=abs(x)
    k=1.0/(1.0+0.2316419*y)
    p=exp(-y**2.0*0.5)/sqrt(2.0*pi)*\
      (0.31938153*k-0.356563782*k**2.0+1.781477937*k**3.0-1.821255978*k**4.0+1.330274429*k**5.0)
    if x>=0.0:
        return 1.0-p
    else:
        return p

def cnd(x,method='default'): # cumulative normal distribution function
    '''
    <The Complete Guide to Option Pricing Formulas>, Chapter 13
    注意：Python中norm.cdf()使用的似乎就是Hart算法
    '''
    if method=='default':
        return norm.cdf(x)
    elif method=='HartAlgo':
        return HartAlgo(x)
    elif method=='PolyApprox6':
        return PolyApprox6(x)
    else:
        return PolyApprox4(x)