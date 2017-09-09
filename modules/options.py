#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 16:16:44 2017

@author: harbes
usage
============
from modules import options
temp=options.euro_option()
temp.value_BSM()
== or ==
from modules import euro_option
temp=euro_option(otype='put')
a.value_BSM()
"""

import numpy as np
import scipy.stats as scis
from scipy.special import comb
from numba import jit
import math
norm=scis.norm
# TODO 尝试利用broadcasting将所有methods向量化；或者丢弃class，全都写成函数的形式，将参数以字典形式储存
class euro_option:
    def __init__(self,St=100.0,K=100.0,t=0.0,T=1.0,r=.05,sigma=.2,otype='call'):
        # (self,St,K,t,T,r,sigma,otype):#
        self.St = St
        self.K = K
        self.t = t
        self.T = T
        self.r = r
        self.sigma = sigma
        self.otype = otype

    def d1f(self):
        return (np.log(self.St / self.K) + (self.r + self.sigma ** 2 / 2) * (self.T-self.t)) / (self.sigma \
                                                                                                * np.sqrt(self.T-self.t))

    def value_BSM(self):
        '''
        using the BSM model to calculate the value of option
   
        return: 
        ======
            the value of option 
        '''
        #norm = scis.norm
        d1 = self.d1f()
        d2 = d1-self.sigma * np.sqrt(self.T-self.t)
        if self.otype == 'call':
            return self.St * norm(0, 1).cdf(d1) - np.exp(-self.r * (self.T-self.t)) * self.K * norm(0, 1).cdf(d2)
        else:
            return np.exp(-self.r * (self.T-self.t)) * self.K * norm(0, 1).cdf(-d2) - self.St * norm(0, 1).cdf(-d1)

#TODO 尝试向量化
    @jit
    def value_CRR(self, M=80):
        '''
        Cox-Ross-Rubinstein European option valuation（binomial）

        return
        ======
        value of option   
        '''
        #import numpy as np
        dt = self.T / M
        df = np.exp(-self.r * dt)  # discount per interval
        # binomial parameters
        u = np.exp(self.sigma * np.sqrt(dt))
        d = 1 / u
        p = (np.exp(self.r * dt) - d) / (u - d)

        S = np.empty(M)
        S[0] = self.St
        for m in range(1, M):
            for n in range(m, 0, -1):
                S[n] = u * S[n - 1]
            S[0] = d * S[0]
        opt = np.empty(M)
        if self.otype=="call":
            for n in range(M):
                opt[n] = S[n] - self.K if S[n] > self.K else 0
        else:
            for n in range(M):
                opt[n] = self.K-S[n] if self.K>S[n] else 0
        for m in range(M - 1, 0, -1):
            for n in range(m):
                opt[n] = (p * opt[n + 1] + (1 - p) * opt[n]) * df
        return opt[0]

    @property
    def delta(self):
        '''
        calculate the delta of the given option
            call: 𝜕C/𝜕S
            put：𝜕P／𝜕S
        :return: 
            the value of delta of the given option
        '''
        #norm = scis.norm
        if self.otype=='call':
            return norm(0, 1).cdf(self.d1f())
        else:
            return norm(0,1).cdf(self.d1f())-1
    @property
    def gamma(self):
        '''
        calculate the gamma of the given option
            call:𝜕 delta／𝜕 S
            put: 𝜕 delta／𝜕 S
        :return: 
        '''
        #norm = scis.norm
        return norm(0,1).pdf(self.d1f())/(self.St*self.sigma*np.sqrt(self.T-self.t))
    @property
    def theta(self):
        '''
        calculate the theta of given option
            call:-𝜕C／𝜕(T-t)
            put: -𝜕P/𝜕(T-t)
        :return: 
        '''
        #norm = scis.norm
        d1=self.d1f()
        d2 = d1 - self.sigma * np.sqrt(self.T - self.t)
        if self.otype=='call':
            return -self.St*norm(0,1).pdf(d1)*self.sigma/2/np.sqrt(self.T-self.t)\
                   -self.r*self.K*np.exp(-self.r*(self.T-self.t))*norm(0,1).cdf(d2)
        else:
            return -self.St * norm(0, 1).pdf(d1) * self.sigma / 2 / np.sqrt(self.T - self.t) \
                   + self.r * self.K * np.exp(-self.r * (self.T - self.t)) * norm(0, 1).cdf(-d2)
    @property
    def rho(self):
        '''
        calculate the rho of the given option
            call:𝜕C/𝜕r
            put: 𝜕P/𝜕r
        :return: 
        '''
        #norm = scis.norm
        d2 = self.d1f() - self.sigma * np.sqrt(self.T - self.t)
        if self.otype=='call':
            return self.K * (self.T - self.t) * np.exp(-self.r * (self.T - self.t)) * norm(0, 1).cdf(d2)
        else:
            return -self.K * (self.T - self.t) * np.exp(-self.r * (self.T - self.t)) * norm(0, 1).cdf(-d2)

    @property
    def vega(self):
        '''
        calculate the vega of the given option
            call:𝜕C/𝜕σ
            put: 𝜕P/𝜕σ
        :return: 
        '''
        #norm=scis.norm
        return self.St*norm(0,1).pdf(self.d1f())*np.sqrt(self.T-self.t)

#TODO American call 好像不对，需要修改
class ame_option():
    def __init__(self,St=100.0,K=100.0,t=0.0,T=1.0,r=.05,sigma=.2,otype='call'):
        self.St = St
        self.K = K
        self.t = t
        self.T = T
        self.r = r
        self.sigma = sigma
        self.otype = otype
    @jit
    def inner_value(self,S,K):
        if self.otype=='call':
            return np.maximum(S-K,0)
        else:
            return np.maximum(K-S,0)


    @jit
    def value_binomial_0(self, M=80):
        dt = (self.T-self.t) / M
        df = np.exp(-self.r * dt)  # discount per interval
        # binomial parameters
        u = np.exp(self.sigma * np.sqrt(dt))
        d = 1 / u
        p = (np.exp(self.r * dt) - d) / (u - d)

        mu = np.arange(M + 1)
        mu = np.resize(mu, (M + 1, M + 1))
        md = np.transpose(mu)
        mu = u ** (mu - md)
        md = d ** md
        S = self.St * mu * md

        h = self.inner_value(S,self.K)
        # h=np.maximum(S-self.K,0) if self.otype=='call' else np.maximum(self.K-S,0)
        V = self.inner_value(S,self.K)
        # V = np.maximum(S - self.K, 0) if self.otype == 'call' else np.maximum(self.K - S, 0)

        z = 0
        for i in range(M - 1, -1, -1):
            C_temp = (p * V[0:M - z, i + 1] + (1 - p) * V[1:M - z + 1, i + 1]) * df
            V[0:M - z, i] = np.where(h[0:M - z, i] > C_temp, h[0:M - z, i], C_temp)
            z += 1
        return V[0, 0]
    @jit
    def value_binomial_1(self,M=50):
        lenS = np.shape(self.St)[0] if np.shape(self.St) else 1
        lenK = np.shape(self.K)[0] if np.shape(self.K) else 1
        lent = np.shape(self.t)[0] if np.shape(self.t) else 1
        lenT = np.shape(self.T)[0] if np.shape(self.T) else 1
        lenr = np.shape(self.r)[0] if np.shape(self.r) else 1
        lensigma = np.shape(self.sigma)[0] if np.shape(self.sigma) else 1
        parameter_tuple=(lenS,lenK,lent,lenT,lenr,lensigma)
        max_len=np.max(parameter_tuple)
        for i,para in enumerate(parameter_tuple):
            if para != 1 and para != max_len:
                print(para)
                raise ValueError('There will be errors in broadcasting, please check the length of every parameter!')
                break
        St = np.resize(self.St, (M + 1, M + 1, lenS))
        K = np.resize(self.K, (M + 1, M + 1, lenK))
        t = np.resize(self.t, (M + 1, M + 1, lent))
        T = np.resize(self.T, (M + 1, M + 1, lenT))
        r = np.resize(self.r, (M + 1, M + 1, lenr))
        sigma = np.resize(self.sigma, (M + 1, M + 1, lensigma))

        dt = (T - t) / M
        df = np.exp(-r * dt)
        u = np.exp(sigma * np.sqrt(dt))
        d = 1 / u
        p = (np.exp(r * dt) - d) / (u - d)

        mu = np.arange(M + 1)
        mu = np.resize(mu, (M + 1, M + 1))
        md = mu.T
        u_d = np.resize(mu - md, (max_len, M + 1, M + 1)).transpose((1, 2, 0))

        mu = u ** (u_d)
        md = d ** np.resize(md, (max_len, M + 1, M + 1)).transpose((1, 2, 0))
        S = St * mu * md
        h=self.inner_value(S,K)
        #h = np.maximum(S - K, 0)  # if self.otype=='call' else np.maximum(self.K-S,0)
        V=self.inner_value(S,K)
        #V = np.maximum(S - K, 0)  # if self.otype == 'call' else np.maximum(self.K - S, 0)
        z = 0
        for i in range(M - 1, -1, -1):
            C_temp = (p[0:M - z, i + 1, :] * V[0:M - z, i + 1, :] + (1 - p[1:M - z + 1, i + 1, :]) \
                      * V[1:M - z + 1,i + 1, :]) * df[0:M - z,i + 1,:]
            V[0:M - z, i, :] = np.where(h[0:M - z, i, :] > C_temp, h[0:M - z, i, :], C_temp)
            z += 1
        return V[0, 0, :]

    @jit
    def value_binomial_2(self, M=50):
        lenS = np.shape(self.St)[0] if np.shape(self.St) else 1
        lenK = np.shape(self.K)[0] if np.shape(self.K) else 1
        lent = np.shape(self.t)[0] if np.shape(self.t) else 1
        lenT = np.shape(self.T)[0] if np.shape(self.T) else 1
        lenr = np.shape(self.r)[0] if np.shape(self.r) else 1
        lensigma = np.shape(self.sigma)[0] if np.shape(self.sigma) else 1
        parameter_tuple = (lenS, lenK, lent, lenT, lenr, lensigma)
        max_len = np.max(parameter_tuple)

        for i, para in enumerate(parameter_tuple):
            if para != 1 and para != max_len:
                print(para)
                raise ValueError('There will be errors in broadcasting, please check the length of every parameter!')
                break

        St = np.resize(self.St, (M + 1, M + 1, lenS))
        K = np.resize(self.K, (M + 1, M + 1, lenK))
        t = np.resize(self.t, (M + 1, M + 1, lent))
        T = np.resize(self.T, (M + 1, M + 1, lenT))
        r = np.resize(self.r, (M + 1, M + 1, lenr))
        sigma = np.resize(self.sigma, (M + 1, M + 1, lensigma))

        dt = (T - t) / M
        df = np.exp(-r * dt)
        u = np.exp(sigma * np.sqrt(dt))
        d = 1 / u
        p = (np.exp(r * dt) - d) / (u - d)

#TODO 尝试使用“深度组合”改写这一段，看看能不能提速,例如np.dstack((tuple(Mu),)*3)

        mu = np.arange(M + 1)
        mu = np.resize(mu, (M + 1, M + 1))
        md = mu.T
        mu = u ** (np.resize(mu - md, (max_len, M + 1, M + 1)).transpose((1, 2, 0)))
        md = d ** np.resize(md, (max_len, M + 1, M + 1)).transpose((1, 2, 0))
        S = St * mu * md
        h = self.inner_value(S, K)
        # h = np.maximum(S - K, 0)  # if self.otype=='call' else np.maximum(self.K-S,0)
        V = self.inner_value(S, K)
        # V = np.maximum(S - K, 0)  # if self.otype == 'call' else np.maximum(self.K - S, 0)
        z = 0
        V_temp = V[0:M - z + 1, M, :]
        for i in range(M - 1, -1, -1):
            C_temp = (p[0:M - z, i + 1, :] * V_temp[0:M - z, :] + (1 - p[1:M - z + 1, i + 1, :]) \
                      * V_temp[1:M - z + 1,:]) * df[0:M - z,i + 1, :]
            V_temp = np.where(h[0:M - z, i, :] > C_temp, h[0:M - z, i, :], C_temp)
            z += 1
        return V_temp[0]
    def value_LSMC(self,I=25000,M=50):
        '''
        此方法似乎不稳健
        :param I: 
        :param M: 
        :return: 
        '''
        dt=self.T/M
        df=np.exp(-self.r*dt)

        #stock price path
        S=self.St*np.exp(np.cumsum((self.r-0.5*self.sigma**2)*dt+self.sigma*np.sqrt(dt)*np.random.standard_normal((M+1,I)),axis=0))
        S[0]=self.St

        h=self.inner_value(S)

        #present value vector
        V=h[-1]
        if self.otype=='call':
            for t in range(M-1,0,-1):
                rg=np.polyfit(S[t],V*df,5)
                C=np.polyval(rg,S[t])
                V=np.where(h[t]<C,h[t],V*df)
        else:
            for t in range(M - 1, 0, -1):
                rg = np.polyfit(S[t], V * df, 5)
                C = np.polyval(rg, S[t])
                V = np.where(h[t] > C, h[t], V * df)
        return df * np.sum(V) / I
    # TODO 期权的动态复制策略
    def dynamicReplication(self):
        pass
    # TODO
    def MonteCarloAntitheticDeltaGammaControlVariates(self):
        pass
    # TODO
    def MCEuroSpreadOption(self):
        pass
    # TODO
    def MCAsianOption(self):
        pass
def udp_binomial(M,mu,dt,sigma,method='CRR'):
    if method =='CRR':
        v=0
    else:
        v=mu
    return math.exp(v*dt+sigma*math.sqrt(dt)),math.exp(v*dt-sigma*math.sqrt(dt)),0.5+(mu-v)/sigma*math.sqrt(dt)/2.0

@jit
def Ame_option_binomial(St,K,r,sigma,T,M,otype='call',method='CRR'):
    mu = r - sigma * sigma / 2.0
    dt = T/ M
    disc = math.exp(-r * dt)
    v = 0 if method=='CRR' else mu
    u,d=math.exp(v*dt+sigma*math.sqrt(dt)),math.exp(v*dt-sigma*math.sqrt(dt))
    p = (1.0/disc- d) / (u - d)
    um=u**np.arange(M+1)
    du=(d/u)**np.arange(M+1)
    S = np.empty(M + 1);S[0]=St*d**M
    for i in range(1,M+1):
        S[i]=S[i-1]*u/d
    o_value=1.0 if otype is 'call' else -1.0
    payoff = np.zeros(M + 1)
    for n in range(M + 1):
        payoff[n] = (S[n]-K)*o_value if (S[n]-K)*o_value >0 else 0.0
    for m in range(M, 1, -1):
        for n in range(m):
            tmp = (p * payoff[n + 1] + (1 - p) * payoff[n]) * disc
            payoff[n] = (St *um[m-1] * du[m-1-n]-K)*o_value if (St * um[m-1] * du[m-1-n]-K)*o_value- tmp>0 else tmp
    return (p*payoff[1]+(1-p)*payoff[0])*disc
@jit
def Ame_option_trinomial(St,K,r,sigma,T,M,otype='call',method='CRR'):
    mu = r - sigma * sigma / 2.0
    lam=math.sqrt(1.5)
    dt = T/ M
    disc = math.exp(-r * dt)
    v = 0 if method=='CRR' else mu
    k1 = math.exp((r - v) * dt)
    k2 = math.exp((2*(r-v)+sigma*sigma)*dt)
    u,m=math.exp(lam*sigma*math.sqrt(dt)),math.exp(v*dt)
    d=1.0/u
    pu=(k2-(d+1)*k1+d)/(u-d)/(u-1.0)
    pd=(k2-(u+1)*k1+u)/(u-d)/(1.0-d)
    pm=1.0-pu-pd
    um=(u*m)**np.arange(M+1)
    dm=(1/u)**np.arange(2*M+1)
    S = np.empty(2*M + 1);S[0]=St*(d*m)**M
    for i in range(1,2*M+1):
        S[i]=S[i-1]*u
    o_value=1.0 if otype is 'call' else -1.0
    payoff = np.zeros(2*M + 1)
    for n in range(2*M + 1):
        payoff[n] = (S[n]-K)*o_value if (S[n]-K)*o_value >0 else 0.0
    for m in range(M, 1, -1):
        for n in range(2*m-1):
            tmp = (pu* payoff[n + 2] +pm*payoff[n+1]+ pd* payoff[n]) * disc
            payoff[n] = (St *um[m-1] * dm[2*(m-1)-n]-K)*o_value if (St * um[m-1] * dm[2*(m-1)-n]-K)*o_value- tmp>0 else tmp
    return (pu* payoff[2] +pm*payoff[1]+ pd* payoff[0]) * disc


@jit
def Ame_option_ExDiff(St, K, r, sigma, T, M, N, otype='call'):
    '''

    :param St:
    :param K:
    :param r:
    :param sigma:
    :param T:
    :param M: 50
    :param N: 100
    :param otype:
    :return:
    '''
    mu = r - sigma * sigma / 2.0
    dt = T / N
    disc = math.exp(-r * dt)
    dx = sigma * math.sqrt(1.5 * dt)
    pu = sigma * sigma * dt / (2.0 * dx * dx) + mu * dt / (2.0 * dx)
    pd = sigma * sigma * dt / (2.0 * dx * dx) - mu * dt / (2.0 * dx)
    pm = 1.0 - pu - pd
    S = St * np.exp(np.arange(-M, M + 1) * dx)
    o_value = 1.0 if otype is 'call' else -1.0
    payoff = np.where((S - K) * o_value > 0, (S - K) * o_value, 0)
    f = np.empty(2 * M + 1)
    f[0] = K * 0.5 * (1 - o_value)
    f[-1] = (S[-1] - K) * 0.5 * (1 + o_value)
    for m in range(N, 1, -1):
        for n in range(1, 2 * M):
            tmp = (pu * payoff[n + 1] + pm * payoff[n] + pd * payoff[n - 1]) * disc
            f[n] = (S[n] - K) * o_value if (S[n] - K) * o_value - tmp > 0 else tmp
        payoff[1:-1] = f[1:-1]
    return (pu * payoff[M + 1] + pm * payoff[M] + pd * payoff[M - 1]) * disc
@jit
def Euro_option_binomial(St,K,r,sigma,T,M=80,otype='call',method='CRR'):
    mu = r - sigma * sigma / 2.0
    dt = T/ M
    disc = math.exp(-r * dt)
    v = 0 if method=='CRR' else mu
    u,d,p=math.exp(v*dt+sigma*math.sqrt(dt)),math.exp(v*dt-sigma*math.sqrt(dt)),0.5-(mu-v)/sigma*math.sqrt(dt)/2.0
    S = np.empty(M + 1);
    S[0] = St * d ** M
    for i in range(1, M + 1):
        S[i] = S[i - 1] * u / d
    o_value = 1 if otype is 'call' else -1
    payoff = np.zeros(M + 1)
    for n in range(M + 1):
        payoff[n] = (S[n] - K) * o_value if (S[n] - K) * o_value > 0 else 0
    for m in range(M, 0, -1):
        for n in range(m):
            payoff[n] = (p * payoff[n + 1] + (1 - p) * payoff[n]) * disc
    return payoff[0]
@jit
def Euro_option_ExDiff(St, K, r, sigma, T, M, N, otype='call'):
    '''

    :param St:
    :param K:
    :param r:
    :param sigma:
    :param T:
    :param M: 50
    :param N: 100
    :param otype:
    :return:
    '''
    mu = r - sigma * sigma / 2.0
    dt = T / N
    disc = math.exp(-r * dt)
    dx = sigma * math.sqrt(1.5 * dt)
    pu = sigma * sigma * dt / (2.0 * dx * dx) + mu * dt / (2.0 * dx)
    pd = sigma * sigma * dt / (2.0 * dx * dx) - mu * dt / (2.0 * dx)
    pm = 1.0 - pu - pd
    S = St * np.exp(np.arange(-M, M + 1) * dx)
    o_value = 1.0 if otype is 'call' else -1.0
    payoff = np.where((S - K) * o_value > 0, (S - K) * o_value, 0)
    f = np.empty(2 * M + 1)
    f[0] = K * 0.5 * (1 - o_value)
    f[-1] = (S[-1] - K) * 0.5 * (1 + o_value)
    for m in range(N, 1, -1):
        for n in range(1, 2 * M):
            f[n]= (pu * payoff[n + 1] + pm * payoff[n] + pd * payoff[n - 1]) * disc
        payoff[1:-1] = f[1:-1]
    return (pu * payoff[M + 1] + pm * payoff[M] + pd * payoff[M - 1]) * disc
# TODO
@jit
def option_binomial(St=100.0,K=100.0,r=0.05,T=1.0,sigma=0.2,M=100,otype='call',American=False):
    dt = T / M;
    sdt = math.sqrt(dt)
    disc = math.exp(-r * dt)
    # u=(1+math.sqrt(math.exp(sigma**2*dt)-1))/disc
    # d=(1-math.sqrt(math.exp(sigma**2*dt)-1))/disc
    # p=0.5
    alpha = (math.exp(-r * dt) + math.exp((r + sigma ** 2) * dt)) / 2.0
    u = alpha + math.sqrt(alpha ** 2 - 1)
    d = 1.0 / u
    p = (1.0 / disc - d) / (u - d)
    S = np.empty(M+1);
    S[0] = St
    if American:
        um = np.empty(M + 1);
        um[0] = 1
        du = np.empty(M + 1);
        du[0] = 1
    for m in range(1, M + 1):
        for n in range(m, 0, -1):
            S[n] = u * S[n - 1]
        S[0] = d * S[0]
        if American:
            um[m] = u * um[m - 1]
            du[m] = du[m - 1] * d / u
    payoff= np.zeros(M+1)
    if otype=='call':
        for n in range(M+1):
            payoff[n] = S[n]-K if S[n]>K else 0
    else:
        for n in range(M+1):
            payoff[n] = K - S[n] if K > S[n] else 0
    for m in range(M, 0, -1):
        for n in range(m):
            payoff[n] = (p * payoff[n + 1] + (1 - p) * payoff[n]) * disc
            if American:
                gain = K - St * um[m] * du[n]
                if gain > payoff[n]:
                    payoff[n] = gain
    return payoff[0]

@jit
def option_binomial_comb(St=100.0,K=100.0,r=0.05,T=1.0,sigma=0.2,M=80,otype='call',American=False):
    dt = T / M;
    sdt = math.sqrt(dt)
    disc = math.exp(-r * dt)
    u=(1+math.sqrt(math.exp(sigma**2*dt)-1))/disc
    d=(1-math.sqrt(math.exp(sigma**2*dt)-1))/disc
    p=0.5
    S = np.empty(M+1);
    for i in range(M+1):
        S[i]=St*u**(i)*d**(M-i)
    payoff = np.zeros(M+1)
    if otype == 'call':
        for n in range(M+1):
            payoff[n] = S[n] - K if S[n] > K else 0
    else:
        for n in range(M+1):
            payoff[n] = K - S[n] if K > S[n] else 0
    if American:
        for m in range(M , 0, -1):
            for n in range(m):
                payoff[n] = (p * payoff[n + 1] + (1 - p) * payoff[n]) * disc
                gain = K - St * u**n * d**(m-1-n)
                if gain > payoff[n]:
                    payoff[n] = gain
        return payoff[0]
    else:
        for m in range(M+1):
            payoff[m] *= comb(M,m)*p**m*(1-p)**(M-m)
        return np.sum(payoff)*math.exp(-r*T)


def option_MC(St=100.0,K=100.0,r=0.05,T=1.0,sigma=0.2,M=100,N=100,JumpLambda=False,JumpKappa=False,otype='call',American=False):
    '''

    :param St:
    :param K:
    :param r:
    :param T:
    :param sigma:
    :param M:
    :param N:
    :param JumpLambda: Bool,number
    :param JumpKappa: Bool,number
    :param otype:
    :param American:
    :return:
    '''
    dt=T/N
    e=np.random.standard_normal((M,N))
    if JumpLambda:
        jump=np.random.poisson(JumpLambda,(M,N))
        mu=r-JumpLambda*JumpKappa
        tmp = mu * dt + sigma * math.sqrt(dt) * e+math.sqrt(dt)*jump
    else:
        mu=r
        tmp = mu * dt + sigma * math.sqrt(dt) * e
    payoff=np.maximum(St*np.exp(np.cumsum(tmp,axis=1)[:,-1])-K,0)
    return math.exp(-r*T)*np.mean(payoff)



class GenRelatedNormal:
    def __init__(self,mu,mat):
        if type(mu) is not np.ndarray or type(mat) is not np.ndarray:
            raise ValueError('请检查并保证输入的数据是数组！')
        self.Mu=mu
        self.VarCov=mat
    def byNumpy(self,nums=100000,seed=None):
        np.random.seed(seed=seed)
        return np.random.multivariate_normal(self.Mu,self.VarCov,nums)
    def byCholesky(self,nums=100000):
        lb=np.linalg.cholesky(self.VarCov)
        m,_=np.shape(self.VarCov)
        mu_=np.zeros(m)
        sigma_=np.eye(m)
        return self.Mu+np.random.multivariate_normal(mu_,sigma_,nums)@lb.T
# TODO 其他算法
@jit
def BrownianBridge(final=None,NumOfsteps=100,T=1):
    delta_t=math.sqrt(T/NumOfsteps)
    e=np.random.standard_normal(NumOfsteps)
    Z=np.zeros(NumOfsteps);
    if final is None:
        Z[-1]=e[-1]*math.sqrt(NumOfsteps)
    else:
        Z[-1]=final
    k = NumOfsteps - 1;
    j = 0;
    for i in range(1,NumOfsteps-1):
        gamma=(i-j)/(k-j) # gamma=1/(k-i+1)
        Z[i]=(1-gamma)*Z[j]+gamma*Z[k]+e[i]*math.sqrt(gamma*(1-gamma)*(k-j)*delta_t)
        j=i
    return Z

