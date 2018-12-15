# Reference:
# <<The Complete Guide to Option Pricing Formulas>> Espen G. Haug


from math import log,exp,sqrt,pi
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

def Delta(flag, S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    delta(call/put)=D_.../D_S (也叫"spot delta", 区别于"strike delta")
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma ** 2.0 / 2.0) * T) / sigma / sqrt(T)
    if flag >0.0: # call
        return cnd(d1)*exp((b-r)*T)
    else: # put, value<0
        return (cnd(d1)-1)*exp((b-r)*T)
def Gamma(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    gamma(call/put)=D_(D_.../D_S)/D_S
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma *sigma / 2.0) * T) / sigma / sqrt(T)
    return exp((b-r)*T-d1*d1*0.5)/S/sigma/sqrt(2.0*pi*T) # gamma(call)=gamma(put)
def Vega(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    vega(call/put)=D_.../D_sigma
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r-q
    d1 = (log(S / K) + (b + sigma * sigma/ 2.0) * T) / sigma / sqrt(T)
    return exp((b - r) * T - d1 * d1 * 0.5) *S*sqrt(0.5*T /pi)  # vega(call)=vega(put)
def Theta(flag, S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    theta(call/put)= -D_.../D_T=D_.../D_t
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma ** 2.0 / 2.0) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    A=-S*sigma*exp((b - r) * T - d1 * d1 * 0.5)/sqrt(8.0*pi*T)
    if flag>0.0:# call
        return A-(b-r)*S*cnd(d1)*exp((b-r)*T)-r*K*cnd(d2)*exp(-r*T)
    else: # put
        return A+(b-r)*S*cnd(-d1)*exp((b-r)*T)+r*K*cnd(-d2)*exp(-r*T)
def Rho(flag, S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    rho(call/put)= D_.../D_r
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d2 = (log(S / K) + (b - sigma * sigma / 2.0) * T) / sigma / sqrt(T)
    if flag>0.0:
        return T*K*exp(-r*T)*cnd(d2)
    else:
        return -T*K*exp(-r*T)*cnd(-d2)
def DdeltaDvol(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    DdeltaDvol(call/put)=vanna(call/put)= D_delta/D_sigma=D_vega/D_S
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma * sigma / 2.0) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    return -exp((b-r)*T-d1*d1*0.5)*d2/sigma/sqrt(2.0*pi)
def DvannaDvol(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    DvannaDvol(call/put)= D_vanna/D_sigma
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma ** 2.0 / 2.0) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    return -exp((b-r)*T-d1*d1*0.5)*d2/sigma/sigma/sqrt(2.0*pi)*(d1*d2-d1/d2-1.0)
def DdeltaDtime(flag, S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    DdeltaDtime(call/put)=- D_delta/D_T=D_delta/D_t
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma *sigma / 2.0) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    A=exp(-0.5*d1*d1)/sqrt(2.0*pi)*(b/sigma/sqrt(T)-d2*0.5/T)
    if flag>0.0:
        return -exp((b-r)*T)*(A+(b-r)*cnd(d1))
    else:
        return -exp((b-r)*T)*(A-(b-r)*cnd(-d1))
def DgammaDvol(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    DgammaDvol(call/put)=zomma(call/put)=D_gamma/D_sigma
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma *sigma / 2.0) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    return exp((b-r)*T-d1*d1*0.5)/S/sigma/sigma/sqrt(2.0*pi*T)*(d1*d2-1.0)
def DgammaDspot(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    DgammaDspot(call/put)=D_gamma/D_S
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma *sigma / 2.0) * T) / sigma / sqrt(T)
    return -exp((b-r)*T-d1*d1*0.5)/S/S/sigma/sqrt(2.0*pi*T)*(1.0+d1/sigma/sqrt(T))
def DgammaDtime(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    DgammaDtime(call/put)=zomma(call/put)=-D_gamma/D_T=D_gamma/D_t
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma *sigma / 2.0) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    return exp((b-r)*T-d1*d1*0.5)/S/sigma/sqrt(2.0*pi*T)*(r-b+b*d1/sigma/sqrt(T)+(1.0-d1*d2)*0.5/T)
def DvegaDvol(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    DvegaDvol(call/put)=vomma(call/put)=D_vega/D_sigma
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma *sigma / 2.0) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    return exp((b - r) * T - d1 * d1 * 0.5) *S*sqrt(0.5*T /pi)*d1*d2/sigma
def DvommaDvol(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    DvommaDvol(call/put)=ultima(call/put)=D_vomma/D_sigma
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma *sigma / 2.0) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    return exp((b - r) * T - d1 * d1 * 0.5) *S*sqrt(0.5*T /pi)*d1*d2/sigma/sigma*(d1*d2-d1/d2-d2/d1-1.0)
def DvegaDtime(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    DvegaDtime(call/put)=-D_vega/D_T=D_vega/D_t
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma * sigma / 2.0) * T) / sigma / sqrt(T)
    d2 = d1 - sigma * sqrt(T)
    return exp((b - r) * T - d1 * d1 * 0.5) *S*sqrt(0.5*T /pi)*(r-b+b*d1/sigma/sqrt(T)-(1.0+d1*d2)/T*0.5)
def VarianceVega(S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    VarianceVega(call/put)=D_.../D_variance (variance=sigma^2)
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d1 = (log(S / K) + (b + sigma *sigma / 2.0) * T) / sigma / sqrt(T)
    return S*exp((b-r)*T-0.5*d1*d1)*sqrt(T/8.0/pi)/sigma
def StrikeDelta(flag, S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    strike_delta(call/put)= D_.../D_K
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d2 = (log(S / K) + (b - sigma * sigma / 2.0) * T) / sigma / sqrt(T)
    if flag>0.0:
        return -exp(-r*T)*cnd(d2)
    else:
        return exp(-r*T)*cnd(-d2)
def StrikeGamma(flag, S, K, r, T, sigma, q=0.0):
    '''
    b=r          gives the BS(1973) stock option model;
    b=r-q        gives the Merton(1973) stock option model with continuous dividend yield q;
    b=0          gives the Black(1976) futures option model;
    b=0, r=0     gives the Asay(1982) margined futures option model;
    b=r-rf       gives the Garman and Kohlhagen(1983) currency option model
    ===================================================
    strike_delta(call/put)= D_(strike_delta)/D_K
    :param flag:
    :param S:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param q:
    :return:
    '''
    assert (T > 0.0) & (sigma > 0.0)
    b = r - q
    d2 = (log(S / K) + (b - sigma * sigma / 2.0) * T) / sigma / sqrt(T)
    return exp(-r*T-d2*d2*0.5)/X/sigma/sqrt(2.0*pi*T)

