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
from math import log,sqrt,exp,pi
assert isinstance(pi,object )
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
class ame_option:
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

#HestonPrice('put',5.0,.05,0.0,-.8,.5,.5,100.0,100.0,.03,.02,.05,1,0.00001,50,0.001)
class HestonOptions:
    def __init__(self,St=100.0,K=100.0,tau=1.0,r=.05,q=.02,otype='call',kappa=5.0,theta=.05,lam=0.0,rho=-.8,sigma=.5,v0=.05):
        '''         
            % Heston (1993) price of a European option.
            % Uses the original formulation by Heston
            % Heston parameters:
            %    kappa  = volatility mean reversion speed parameter
            %    theta  = volatility mean reversion level parameter
            %    lamb = risk parameter
            %    rho    = correlation between two Brownian motions
            %    sigma  = volatility of variance
            %    v0     = initial variance
            % Option features.
            %    otype = 'C'all or 'P'ut
            %    K = strike price
            %    St = spot price
            %    r = risk free rate
            %    q = dividend yield
            %    T = maturity
            % Integration features
            %    L = lower limit
            %    U = upper limit
            %    dphi = integration increment
            % 此处主要是用基于梯形法则的数值积分
            % example: HestonPrice('put',5.0,.05,0.0,-.8,.5,.5,100.0,100.0,.03,.02,.05,1,0.00001,50,0.001)
        '''
        self.St=St
        self.K=K
        self.tau=tau
        self.r=r
        self.q=q
        self.otype=otype
        self.kappa=kappa
        self.theta=theta
        self.lam=lam
        self.rho=rho
        self.sigma=sigma
        self.v0=v0

    def HestonIntegrand_Trap(self,phi,Pnum=0):
        '''
        phi = integration variable
        '''
        a = self.kappa * self.theta
        u = np.array([0.5, -.5])[Pnum]
        b = np.array([self.kappa + self.lam - self.rho * self.sigma, self.kappa + self.lam])[Pnum]
        d = np.sqrt((self.rho * self.sigma * 1.0j * phi - b) ** 2 - self.sigma * self.sigma * (2 * u * 1.0j * phi - phi * phi))
        g = (b - self.rho * self.sigma * 1.0j * phi + d) / (b - self.rho * self.sigma * 1.0j * phi - d)
        c = 1.0 / g
        D = (b - self.rho * self.sigma * 1.0j * phi - d) / self.sigma / self.sigma * (1 - np.exp(-d * self.tau)) / (1 - c * np.exp(-d * self.tau))
        G = (1 - c * np.exp(-d * self.tau)) / (1 - c)
        C = (self.r - self.q) * 1.0j * phi * self.tau + a / self.sigma / self.sigma * ((b - self.rho * self.sigma * 1.0j * phi - d) * self.tau - 2 * np.log(G))
        f = np.exp(C + D * self.v0 + 1.0j * phi * log(self.St))
        return (np.exp(-1.0j * phi * log(self.K)) * f / 1.0j / phi).real
    def HestonIntegrand_Orig(self,phi,Pnum=0):
        a = self.kappa * self.theta
        u = np.array([0.5, -.5])[Pnum]
        b = np.array([self.kappa + self.lam - self.rho * self.sigma, self.kappa + self.lam])[Pnum]
        d = np.sqrt((self.rho * self.sigma * 1.0j * phi - b) ** 2 - self.sigma * self.sigma * (2 * u * 1.0j * phi - phi * phi))
        g = (b - self.rho * self.sigma * 1.0j * phi + d) / (b - self.rho * self.sigma * 1.0j * phi - d)
        G = (1 - g * np.exp(d * self.tau)) / (1 - g)
        C = (self.r - self.q) * 1.0j * phi * self.tau + a / self.sigma / self.sigma * ((b - self.rho * self.sigma * 1.0j * phi + d) * self.tau - 2 * np.log(G))
        D = (b - self.rho * self.sigma * 1.0j * phi + d) / self.sigma / self.sigma * (1 - np.exp(d * self.tau)) / (1 - g * np.exp(d * self.tau))
        f = np.exp(C + D * self.v0 + 1.0j * phi * log(self.St))
        return (np.exp(-1.0j * phi * log(self.K)) * f / 1.0j / phi).real
    def HestonCharacteristicFunction(self, phi):
        a = self.kappa * self.theta
        u = -0.5
        b = self.kappa + self.lam
        d = np.sqrt(
            (self.rho * self.sigma * 1.0j * phi - b) ** 2 - self.sigma * self.sigma * (2 * u * 1.0j * phi - phi * phi))
        g = (b - self.rho * self.sigma * 1.0j * phi + d) / (b - self.rho * self.sigma * 1.0j * phi - d)
        G = (1 - g * np.exp(d * self.tau)) / (1 - g)
        C = (self.r - self.q) * 1.0j * phi * self.tau + a / self.sigma / self.sigma * (
        (b - self.rho * self.sigma * 1.0j * phi + d) * self.tau - 2 * np.log(G))
        D = (b - self.rho * self.sigma * 1.0j * phi + d) / self.sigma / self.sigma * (1 - np.exp(d * self.tau)) / (
        1 - g * np.exp(d * self.tau))
        return np.exp(C + D * self.v0 + 1.0j * phi * log(self.St))
    @jit
    def HestonPrice(self,Lphi=0.0001,Uphi=50.0,dphi=0.1,Trap=1):
        '''
            % Heston (1993) price of a European option.
            % Uses the original formulation by Heston
            % Heston parameters:
            %    kappa  = volatility mean reversion speed parameter
            %    theta  = volatility mean reversion level parameter
            %    lambda = risk parameter
            %    rho    = correlation between two Brownian motions
            %    sigma  = volatility of variance
            %    v0     = initial variance
            % Option features.
            %    PutCall = 'C'all or 'P'ut
            %    K = strike price
            %    S = spot price
            %    r = risk free rate
            %    q = dividend yield
            %    T = maturity
            % Integration features
            %    L = lower limit
            %    U = upper limit
            %    dphi = integration increment
            % 此处主要是用基于梯形法则的数值积分
            '''
        phi = np.arange(Lphi, Uphi + .001, dphi)
        N = len(phi)
        int1 = np.empty(N)
        int2 = np.empty(N)
        if Trap==1:
            Heston_Inte=self.HestonIntegrand_Trap
        else:
            Heston_Inte=self.HestonIntegrand_Orig
        int1 = Heston_Inte(phi, Pnum=0)
        int2 = Heston_Inte(phi, Pnum=1)
        I1 = np.trapz(int1) * dphi
        I2 = np.trapz(int2) * dphi
        P1 = 0.5 + 1.0 / np.pi * I1
        P2 = 0.5 + 1.0 / np.pi * I2
        if self.otype == 'call':
            return self.St * np.exp(-self.q * self.tau) * P1 - self.K * np.exp(-self.r * self.tau) * P2
        else:
            return self.K * np.exp(-self.r * self.tau) * (1.0 - P2) - self.St * np.exp(-self.q * self.tau) * (1.0 - P1)
    def HestonPriceTrapezoidal(self,Lphi=0.00001,Uphi=50.0,dphi=0.1):
        phi=np.arange(Lphi,Uphi,dphi)
        n=len(phi)
        h=(phi[n-1]-phi[0])/n
        weights=np.ones(n);weights[0]=0.5;weights[n-1]=0.5
        int1=h*weights*self.HestonIntegrand_Trap(phi,Pnum=0)
        int2=h*weights*self.HestonIntegrand_Trap(phi,Pnum=1)
        I1=np.sum(int1)
        I2=np.sum(int2)
        P1=0.5+1.0/pi*I1
        P2=0.5+1.0/pi*I2
        if self.otype == 'call':
            return self.St * np.exp(-self.q * self.tau) * P1 - self.K * np.exp(-self.r * self.tau) * P2
        else:
            return self.K * np.exp(-self.r * self.tau) * (1.0 - P2) - self.St * np.exp(-self.q * self.tau) * (1.0 - P1)

    def HestonInteConsol_Trap(self,phi):
        '''
        phi = integration variable
        Pnum = 1 or 2 (for the probabilities)
        Trap = 1 "Little Trap" formulation
               0  Original Heston formulation

        '''
        a = self.kappa * self.theta
        u = 0.5
        b = self.kappa + self.lam - self.rho * self.sigma
        d = np.sqrt((self.rho * self.sigma * 1.0j * phi - b) ** 2 - self.sigma * self.sigma * (2 * u * 1.0j * phi - phi * phi))
        g = (b - self.rho * self.sigma * 1.0j * phi + d) / (b - self.rho * self.sigma * 1.0j * phi - d)
        c = 1.0 / g
        D = (b - self.rho * self.sigma * 1.0j * phi - d) / self.sigma / self.sigma * (1 - np.exp(-d * self.tau)) / (1 - c * np.exp(-d * self.tau))
        G = (1 - c * np.exp(-d * self.tau)) / (1 - c)
        C = (self.r - self.q) * 1.0j * phi * self.tau + a / self.sigma / self.sigma * ((b - self.rho * self.sigma * 1.0j * phi - d) * self.tau - 2 * np.log(G))
        f1 = np.exp(C + D * self.v0 + 1.0j * phi * log(self.St))
        u =  -.5
        b =  self.kappa + self.lam
        d = np.sqrt(
            (self.rho * self.sigma * 1.0j * phi - b) ** 2 - self.sigma * self.sigma * (2 * u * 1.0j * phi - phi * phi))
        g = (b - self.rho * self.sigma * 1.0j * phi + d) / (b - self.rho * self.sigma * 1.0j * phi - d)
        c = 1.0 / g
        D = (b - self.rho * self.sigma * 1.0j * phi - d) / self.sigma / self.sigma * (1 - np.exp(-d * self.tau)) / (
        1 - c * np.exp(-d * self.tau))
        G = (1 - c * np.exp(-d * self.tau)) / (1 - c)
        C = (self.r - self.q) * 1.0j * phi * self.tau + a / self.sigma / self.sigma * (
        (b - self.rho * self.sigma * 1.0j * phi - d) * self.tau - 2 * np.log(G))
        f2 = np.exp(C + D * self.v0 + 1.0j * phi * log(self.St))
        return (np.exp(-1.0j * phi * log(self.K)) / 1.0j / phi*(self.St*exp(-self.q*self.tau)*f1-self.K*exp(-self.r*self.tau)*f2)).real
    def HestonInteConsol_Orig(self,phi):
        a = self.kappa * self.theta
        u = 0.5
        b = self.kappa + self.lam - self.rho * self.sigma, self.kappa + self.lam
        d = np.sqrt((self.rho * self.sigma * 1.0j * phi - b) ** 2 - self.sigma * self.sigma * (2 * u * 1.0j * phi - phi * phi))
        g = (b - self.rho * self.sigma * 1.0j * phi + d) / (b - self.rho * self.sigma * 1.0j * phi - d)
        G = (1 - g * np.exp(d * self.tau)) / (1 - g)
        C = (self.r - self.q) * 1.0j * phi * self.tau + a / self.sigma / self.sigma * ((b - self.rho * self.sigma * 1.0j * phi + d) * self.tau - 2 * np.log(G))
        D = (b - self.rho * self.sigma * 1.0j * phi + d) / self.sigma / self.sigma * (1 - np.exp(d * self.tau)) / (1 - g * np.exp(d * self.tau))
        f1 = np.exp(C + D * self.v0 + 1.0j * phi * log(self.St))
        u = -0.5
        b =self.kappa + self.lam
        d = np.sqrt(
            (self.rho * self.sigma * 1.0j * phi - b) ** 2 - self.sigma * self.sigma * (2 * u * 1.0j * phi - phi * phi))
        g = (b - self.rho * self.sigma * 1.0j * phi + d) / (b - self.rho * self.sigma * 1.0j * phi - d)
        G = (1 - g * np.exp(d * self.tau)) / (1 - g)
        C = (self.r - self.q) * 1.0j * phi * self.tau + a / self.sigma / self.sigma * (
        (b - self.rho * self.sigma * 1.0j * phi + d) * self.tau - 2 * np.log(G))
        D = (b - self.rho * self.sigma * 1.0j * phi + d) / self.sigma / self.sigma * (1 - np.exp(d * self.tau)) / (
        1 - g * np.exp(d * self.tau))
        f2 = np.exp(C + D * self.v0 + 1.0j * phi * log(self.St))
        return (np.exp(-1.0j * phi * log(self.K)) / 1.0j / phi*(self.St*exp(-self.q*self.tau)*f1-self.K*exp(-self.r*self.tau)*f2)).real
    def HestonPriceConsol(self,Lphi=0.0001,Uphi=50.0,dphi=0.1,Trap=1):
        '''
            % Heston (1993) price of a European option.
            % Uses the original formulation by Heston
            % Heston parameters:
            %    kappa  = volatility mean reversion speed parameter
            %    theta  = volatility mean reversion level parameter
            %    lambda = risk parameter
            %    rho    = correlation between two Brownian motions
            %    sigma  = volatility of variance
            %    v0     = initial variance
            % Option features.
            %    PutCall = 'C'all or 'P'ut
            %    K = strike price
            %    S = spot price
            %    r = risk free rate
            %    q = dividend yield
            %    T = maturity
            % Integration features
            %    L = lower limit
            %    U = upper limit
            %    dphi = integration increment
            % 此处主要是用基于梯形法则的数值积分
             '''
        phi = np.arange(Lphi, Uphi + .001, dphi)
        n = len(phi)
        inte = np.empty(n)
        if Trap==1:
            Heston_Inte=self.HestonInteConsol_Trap
        else:
            Heston_Inte=self.HestonInteConsol_Orig
        inte = Heston_Inte(phi)
        I = np.trapz(inte) * dphi
        if self.otype == 'call':
            return self.St * np.exp(-self.q * self.tau) *0.5 - self.K * np.exp(-self.r * self.tau) * 0.5 + I/pi
        else:
            return self.K * np.exp(-self.r * self.tau) *0.5 - self.St * np.exp(-self.q * self.tau) *0.5+I/pi

    def ExampleIntegrandDiscontinuity():
        phi=np.arange(0.00001,10.01,0.05)
        n=len(phi)
        A=np.empty((n,2))
        A[:,0]=HestonIntegrand_Orig(phi,10.0,0.05,0.0,-0.9,0.75,3.0,75.0,100.0,0.0,0.0,0.05,0)
        A[:, 1] =HestonIntegrand_Orig(phi,10.0,0.05,0.0,-0.9,0.09,1.0,75.0,100.0,0.0,0.0,0.05,0)
        plt.plot(phi,A[:,0],label='sigma=0.75,tau=3.0')
        plt.plot(phi,A[:,1],label='sigma=0.09,tau=1.0')
        plt.legend()
        plt.show()
    def ExampleIntegrandOscillation():
        phi = np.arange(0.00001, 100.01, 0.5)
        n = len(phi)
        A = np.empty((n, 2))
        A[:, 0] =HestonIntegrand_Orig(phi,10.0,0.01,0.0,-0.9,0.175,1.0/52.0,10.0,7.0,0.0,0.0,0.01,0)
        A[:, 1] =HestonIntegrand_Orig(phi,10.0,0.07,0.0,-0.9,0.09, 1.0,     10.0,10.0,0.0,0.0,0.01,0)
        plt.plot(phi, A[:, 0], label='sigma=0.175,tau=1/52,v0=0.01,theta=0.01,St=7')
        plt.plot(phi, A[:, 1], label='sigma=0.09,tau=1.0,v0=0.07,theta=0.07,St=10')
        plt.legend()
        plt.show()
    def ExampleAlbrecher_vs_Heston(Pnum=0):
        HOpt=HestonOptions(St=100.0,K=100.0,tau=5.0,r=0.035,q=0.0,kappa=1.5768,sigma=0.5751,rho=-0.5711,theta=0.0398,v0=0.0175,lam=0.0)
        phi=np.arange(0.0001,10.01,0.1)
        n=len(phi)
        A=np.empty((n,2))
        A[:,0]=HestonIntegrand_Orig(phi,1.5768,0.0398,0.0,0.5711,0.5751,5.0,100.0,100.0,0.035,0.0,0.0175,Pnum)
        A[:,1]=HestonIntegrand_Trap(phi,1.5768,0.0398,0.0,0.5711,0.5751,5.0,100.0,100.0,0.035,0.0,0.0175,Pnum)
        plt.plot(phi,A[:,0],label='Heston Integrand')
        plt.plot(phi,A[:,1],label='Albrecher Integrand')
        plt.legend()
        plt.show()
    def ExampleEffect_of_Correlation_on_Terminal_Price_Density(export=False):
        dphi=0.05
        phi=np.arange(0.00001,100.01,dphi)
        x=np.arange(4.3,4.9,0.005)
        paras={'St':100.0,'K':100.0,'tau':0.5,'r':0.0,'q':0.0,'kappa':2.0,'theta':0.01,'sigma':0.1,'v0':0.01}
        f1=np.empty(len(x))
        f2=np.empty(len(x))
        f3=np.empty(len(x))
        for i in range(len(x)):
            intneg=(np.exp(-1.0j*phi*x[i])*HestonCharacteristicFunction(phi,2.0,0.01,0.0,-0.8,0.1,0.5,100.0,100.0,0.0,0.0,0.01)).real
            intzero=(np.exp(-1.0j*phi*x[i])*HestonCharacteristicFunction(phi,2.0,0.01,0.0,0.0,0.1,0.5,100.0,100.0,0.0,0.0,0.01)).real
            intpos=(np.exp(-1.0j*phi*x[i])*HestonCharacteristicFunction(phi,2.0,0.01,0.0,0.8,0.1,0.5,100.0,100.0,0.0,0.0,0.01)).real
            f1[i]=np.trapz(intneg)/math.pi*dphi
            f2[i]=np.trapz(intzero)/math.pi*dphi
            f3[i]=np.trapz(intpos)/math.pi*dphi
        plt.plot(x,f1,label='rho=-0.8')
        plt.plot(x,f2,label='rho=0')
        plt.plot(x,f3,label='rho=0.8')
        plt.legend()
        plt.show()
        if export:
            return x,f1,f2,f3
    def ExampleEffect_of_Vol_of_Variance_on_Terminal_Price_Density(export=False):
        dphi=0.05
        phi=np.arange(0.00001,100.01,dphi)
        x=np.arange(4.3,4.9,0.005)
        paras={'St':100.0,'K':100.0,'tau':0.5,'r':0.0,'q':0.0,'kappa':2.0,'theta':0.01,'v0':0.01,'rho':0.0}
        f1=np.empty(len(x))
        f2=np.empty(len(x))
        f3=np.empty(len(x))
        for i in range(len(x)):
            intzero=(np.exp(-1.0j*phi*x[i])*HestonCharacteristicFunction(phi,2.0,0.01,0.0,0.0,0.00001,0.5,100.0,100.0,0.0,0.0,0.01)).real
            intlow=(np.exp(-1.0j*phi*x[i])*HestonCharacteristicFunction(phi,2.0,0.01,0.0,0.0,0.2,0.5,100.0,100.0,0.0,0.0,0.01)).real
            inthigh=(np.exp(-1.0j*phi*x[i])*HestonCharacteristicFunction(phi,2.0,0.01,0.0,0.0,0.4,0.5,100.0,100.0,0.0,0.0,0.01)).real
            f1[i]=np.trapz(intzero)/math.pi*dphi
            f2[i]=np.trapz(intlow)/math.pi*dphi
            f3[i]=np.trapz(inthigh)/math.pi*dphi
        plt.plot(x,f1,label='sigma=0')
        plt.plot(x,f2,label='sigma=0.2')
        plt.plot(x,f3,label='sigma=0.4')
        plt.legend()
        plt.show()
        if export:
            return x,f1,f2,f3
    def ExampleEffect_of_Correlation_on_Heston_and_BlackScholes_Prices():
        '''
        一个值得深入思考的现象：
import numpy as np
from modules import HestonOptions as Hoptions
import matplotlib.pyplot as plt
x,f1,f2,f3=Hoptions.ExampleEffectOfCorrelation(True)

plt.plot(x,np.cumsum(f3-f2).cumsum(),label='survive diff : positive-normal')
plt.legend()
plt.show()
        '''
        S=np.arange(70,141,0.05)
        volneg=0.071037274323352*sqrt(2)
        volpos=0.070386797400082*sqrt(2)
        BScallpos=euro_option(S,r=0.0,sigma=volpos,T=0.5).value_BSM()
        BScallneg=euro_option(S,r=0.0,sigma=volneg,T=0.5).value_BSM()
        Hcallpos=np.empty_like(S)
        Hcallneg=np.empty_like(S)
        for i in range(len(S)):
            Hcallpos[i]=HestonPriceConsol(2.0,0.01,0.0,0.5,0.1,0.5,100.0,S[i],0.0,0.0,0.01,1,0.00001,75.001,0.1)
            Hcallneg[i]=HestonPriceConsol(2.0,0.01,0.0,-0.5,0.1,0.5,100.0,S[i],0.0,0.0,0.01,1,0.00001,75.001,0.1)
        plt.plot(S,Hcallpos-BScallpos,label='rho=+0.5')
        plt.plot(S,Hcallneg-BScallneg,label='rho=-0.5')
        plt.legend()
        plt.show()
    def ExampleEffect_of_Sigma_on_Heston_and_BlackScholes_Prices():
        S = np.arange(70, 141, 0.05)
        volzero = 0.070712338973920 * sqrt(2)
        BScall = euro_option(S, r=0.0, sigma=volzero, T=0.5).value_BSM()
        Hcallhigh = np.empty_like(S)
        Hcalllow= np.empty_like(S)
        for i in range(len(S)):
            Hcallhigh[i] = HestonPriceConsol(2.0,0.01,0.0,0.0,0.2,0.5,100.0,S[i],0.0,0.0,0.01,1,0.00001,90.001,0.1)
            Hcalllow[i] = HestonPriceConsol(2.0,0.01,0.0,0.0,0.1,0.5,100.0,S[i],0.0,0.0,0.01,1,0.00001,90.001,0.1)
        plt.plot(S, Hcallhigh - BScall, label='sigma=0.2')
        plt.plot(S, Hcalllow - BScall, label='sigma=0.1')
        plt.legend()
        plt.show()
    def ExampleEffect_of_Heston_Parameters_on_Implied_Volatility():
        #H_paras={'r':0.05,'q':0.0,'tau':0.25,'rho':0.0,'kappa':2.0,'theta':0.01}
        K=np.arange(95,110.001,0.5)
        n=len(K)
        KK=np.resize(K,(3,n)).T
        Uphi=50.01
        dphi=0.1
        sigma=(0.35,0.55,0.8)
        Hcall=np.empty((n,3))
        for i in range(3):
            for j in range(n):
                Hcall[j,i]=HestonPriceConsol(2.0,0.01,0.0,0.0,sigma[i],0.25,K[j],100.0,0.05,0.0,0.05,1,0.00001,Uphi,dphi)
        IVsigma=imp_vol_Formula(Hcall,100.0,KK,0.05,0.25,0.5,'call')

        rho = (-0.25, 0, 0.25)
        for i in range(3):
            for j in range(n):
                Hcall[j, i] = HestonPriceConsol(2.0, 0.01, 0.0, rho[i],0.25, 0.25, K[j], 100.0, 0.05, 0.0, 0.05,
                                                1, 0.00001,Uphi,dphi)
        IVrho = imp_vol_Formula(Hcall, 100.0, KK, 0.05, 0.25, 0.5, 'call')

        v0 = (0.02, 0.015, 0.01)
        for i in range(3):
            for j in range(n):
                Hcall[j, i] = HestonPriceConsol(2.0, 0.01, 0.0, 0.0, 0.5, 0.25, K[j], 100.0, 0.05, 0.0, v0[i],
                                                1, 0.00001,Uphi,dphi)
        IVv0 = imp_vol_Formula(Hcall, 100.0, KK, 0.05, 0.25, 0.5, 'call')

        kappa=(5.0,2.0,1.0)
        for i in range(3):
            for j in range(n):
                Hcall[j, i] = HestonPriceConsol(kappa[i], 0.01, 0.0, 0.0, 0.2, 0.25, K[j], 100.0, 0.05, 0.0, 0.01,
                                                1, 0.00001,Uphi,dphi)
        IVkappa = imp_vol_Formula(Hcall, 100.0, KK, 0.05, 0.25, 0.5, 'call')

        theta=(0.02,0.015,0.01)
        for i in range(3):
            for j in range(n):
                Hcall[j, i] = HestonPriceConsol(2.0, theta[i], 0.0, 0.0, 0.5, 0.25, K[j], 100.0, 0.05, 0.0, 0.05,
                                                1, 0.00001,Uphi,dphi)
        IVtheta = imp_vol_Formula(Hcall, 100.0, KK, 0.05, 0.25, 0.5, 'call')
        line=('k-','r-','g-')
        fig=plt.figure()
        ax=fig.add_subplot(3,2,1)
        for i in range(3):
            ax.plot(K,IVsigma[:,i],line[i],label='sigma={}'.format(sigma[i]))
        ax.legend()
        ax = fig.add_subplot(3, 2, 2)
        for i in range(3):
            ax.plot(K,IVrho[:,i],line[i],label='rho={}'.format(rho[i]))
        ax.legend()
        ax = fig.add_subplot(3, 2, 3)
        for i in range(3):
            ax.plot(K, IVv0[:, i], line[i],label='v0={}'.format(v0[i]))
        ax.legend()
        ax = fig.add_subplot(3, 2, 4)
        for i in range(3):
            ax.plot(K, IVkappa[:, i], line[i],label='kappa={}'.format(kappa[i]))
        ax.legend()
        ax = fig.add_subplot(3, 2, 5)
        for i in range(3):
            ax.plot(K, IVtheta[:, i],line[i], label='theta={}'.format(theta[i]))
        ax.legend()
        plt.show()
    def ExamplePlot_ImpliedVol_and_LocalVol_surface():
        v0=0.4*0.4
        theta=0.2*0.2
        rho=-0.6
        kappa=1.5
        sigma=5.1
        lam=0
        S=100.0
        r=0.0
        q=0.0
        Lphi=0.00001
        Uphi=50.001
        dphi=0.1

        K=np.arange(95,105.01,0.5)
        tau=np.arange(0.25,1.0301,0.04)
        Nk=len(K)
        Nt=len(tau)
        tt,KK= np.meshgrid(tau,K)
        price=v_HestonPriceConsol(kappa, theta, lam, rho, sigma, tt, KK, S, r, q, v0, 1, Lphi, Uphi,dphi)
        LV = LocalVariance_FD(kappa, theta, lam, rho, sigma, tt, KK, S, r, q, v0, 1, Lphi, Uphi, dphi, 0.001,0.0125)
        IV=imp_vol_Formula(price,S,KK,r,tt,0.4,'call')
        fig=plt.figure()
        ax=Axes3D(fig)
        ax.plot_surface(KK,tt,IV,rstride=1,cstride=1,cmap='rainbow')
        ax.plot_surface(KK,tt,LV)
        plt.show()
    def ExampleHeston_Local_Vol_Gatheral_Example():
        v0=0.57942364171876
        theta=0.00701250929769
        kappa=0.00000645051426
        sigma=0.63913150446556
        rho=-0.80541286329551
        lam=0.0
        T=222.0/365.0
        r=0.0
        q=0.0
        Lphi=1e-5
        Uphi=50.01
        dphi=0.1
        Trap=1
        PutCall = 'put'

        K=np.arange(15,65.01)
        S=30.67
        xT=np.log(K/S)
        # local vol
        dt=1e-5
        dK=1e-1
        LVFD=LocalVariance_FD(kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,Trap,Lphi,Uphi,dphi,dt,dK)
        LVAN=v_LocalVariance_Analytic(kappa,theta,lam,rho,sigma,T,K,S,r,v0,Trap,Lphi,Uphi,dphi)
        LVAP=LocalVariance_Approx(kappa,theta,rho,sigma,T,K,S,v0)

        # implied vol
        price=v_HestonPriceConsol(kappa, theta, lam, rho, sigma, T, K, S, r, q, v0, Trap, Lphi, Uphi,dphi,PutCall)
        IV=imp_vol_Formula(price,S,K,r,T,0.5,PutCall)

        plt.plot(K,LVAP,'r-',label='Approx')
        plt.plot(K,LVAN,'p-',label='Analytic')
        plt.plot(K,LVFD,'kx',label='Finite Diff')
        plt.plot(K,IV,'rx-',label='Implied Vol')
        plt.legend()
        plt.show()
    def ExampleHeston_Price_Using_Heston_and_Attari():
        S=30.0
        K=20.0
        T=1.0/12.0
        r=0.01
        q=0.0
        kappa=1.4
        theta=0.05
        sigma=0.3
        v0=0.05
        rho=-0.8
        lam=0.0
        Trap=1
        PutCall='call'
        if Trap==1:
            HestonIntegrand=HestonIntegrand_Trap
            AttariIntegrand=AttariIntegrand_Trap
        else:
            HestonIntegrand = HestonIntegrand_Orig
            AttariIntegrand = AttariIntegrand_Orig
        phi=np.arange(1e-6,20.0001,0.01)
        Heston=HestonIntegrand(phi,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,1)
        Attari=AttariIntegrand(phi,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0)
        zero=np.zeros(len(phi))
        plt.plot(phi,Heston,label='Heston')
        plt.plot(phi,Attari,label='Attari')
        plt.plot(phi,zero,'k-')
        plt.legend()
        plt.show()
    def ExampleComparison_of_Heston_and_Carr_Madan_Calls(PutCall='call'):
        S=100.0
        K=100.0
        T=0.5
        r=0.1
        q=0.07
        kappa=2.0
        theta=0.06
        sigma=0.1
        v0=0.06
        rho=-0.7
        lam=0.0
        Trap=1
        alpha=1.75
        Lphi=1e-5
        Uphi=50.001
        dphi=0.1
        callHeston=HestonPriceConsol(kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,Trap,Lphi,Uphi,dphi,PutCall)
        priceCarrMadan=CarrMadanPrice(alpha,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,Lphi,Uphi,dphi,PutCall)
        print('Heston call price:{}'.format(callHeston))
        print('Carr&Madan call price:{}'.format(priceCarrMadan))
    def ExampleCarrMadan_OTM_Damped_Undamped_Integrands_Figure():
        S=1.0
        K=0.96
        tau=1.0/52.0
        r=0.03
        q=0.0
        kappa=2.0
        theta=0.25
        sigma=0.3
        lam=0
        v0=0.05
        rho=-0.8
        Trap=1
        alpha=1.1

        U=np.arange(-170,170.001)
        T=np.arange(1,5)/52.0
        UU,TT=np.meshgrid(U,T)
        undamped=CarrMadanIntegrandOTM(UU,kappa,theta,lam,rho,sigma,TT,K,S,r,q,v0)
        a=CarrMadanIntegrandOTM(UU-1.0j*alpha,kappa,theta,lam,rho,sigma,TT,K,S,r,q,v0)
        b=CarrMadanIntegrandOTM(UU+1.0j*alpha,kappa,theta,lam,rho,sigma,TT,K,S,r,q,v0)
        damped=(a-b)*0.5
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot_surface(UU, TT, undamped, rstride=1, cstride=1, cmap='rainbow')
        ax.plot_surface(UU, TT, damped)
        plt.show()






def HestonIntegrand_Trap(phi,kappa,theta,lam,rho,sigma,tau,K,S,r,q,v0,Pnum=0):
    '''
    % Returns the integrand for the risk neutral probabilities P1 and P2.
    % phi = integration variable
    XXXXXX Pnum = 1 or 2 (for the probabilities)XXXXX
    % Heston parameters:
    %    kappa  = volatility mean reversion speed parameter
    %    theta  = volatility mean reversion level parameter
    %    lambda = risk parameter
    %    rho    = correlation between two Brownian motions
    %    sigma  = volatility of variance
    %    v      = initial variance
    % Option features.
    %    PutCall = 'C'all or 'P'ut
    %    K = strike price
    %    S = spot price
    %    r = risk free rate
    %    q = dividend yield
    %    Trap = 1 "Little Trap" formulation
    %           0  Original Heston formulation
    '''
    a=kappa*theta
    u=np.array([0.5,-.5])[Pnum]
    b=np.array([kappa+lam-rho*sigma,kappa+lam])[Pnum]
    d=np.sqrt((rho*sigma*1.0j*phi-b)**2-sigma*sigma*(2*u*1.0j*phi-phi*phi))
    g=(b-rho*sigma*1.0j*phi+d)/(b-rho*sigma*1.0j*phi-d)
    c=1.0/g
    D=(b-rho*sigma*1.0j*phi-d)/sigma/sigma*(1-np.exp(-d*tau))/(1-c*np.exp(-d*tau))
    G=(1-c*np.exp(-d*tau))/(1-c)
    C=(r-q)*1.0j*phi*tau+a/sigma/sigma*((b-rho*sigma*1.0j*phi-d)*tau-2*np.log(G))
    f=np.exp(C+D*v0+1.0j*phi*log(S))
    return (np.exp(-1.0j*phi*np.log(K))*f/1.0j/phi).real
def HestonIntegrand_Orig(phi,kappa,theta,lam,rho,sigma,tau,K,S,r,q,v0,Pnum=0):
    '''
    % Returns the integrand for the risk neutral probabilities P1 and P2.
    % phi = integration variable
    XXXXXX Pnum = 1 or 2 (for the probabilities)XXXXX
    % Heston parameters:
    %    kappa  = volatility mean reversion speed parameter
    %    theta  = volatility mean reversion level parameter
    %    lambda = risk parameter
    %    rho    = correlation between two Brownian motions
    %    sigma  = volatility of variance
    %    v      = initial variance
    % Option features.
    %    PutCall = 'C'all or 'P'ut
    %    K = strike price
    %    S = spot price
    %    r = risk free rate
    %    q = dividend yield
    %    Trap = 1 "Little Trap" formulation
    %           0  Original Heston formulation
    '''
    a=kappa*theta
    u=np.array([0.5,-.5])[Pnum]
    b=np.array([kappa+lam-rho*sigma,kappa+lam])[Pnum]
    d=np.sqrt((rho*sigma*1.0j*phi-b)**2-sigma*sigma*(2*u*1.0j*phi-phi*phi))
    g=(b-rho*sigma*1.0j*phi+d)/(b-rho*sigma*1.0j*phi-d)
    G=(1-g*np.exp(d*tau))/(1-g)
    C=(r-q)*1.0j*phi*tau+a/sigma/sigma*((b-rho*sigma*1.0j*phi+d)*tau-2*np.log(G))
    D=(b-rho*sigma*1.0j*phi+d)/sigma/sigma*(1-np.exp(d*tau))/(1-g*np.exp(d*tau))
    f=np.exp(C+D*v0+1.0j*phi*log(S))
    return (np.exp(-1.0j*phi*np.log(K))*f/1.0j/phi).real
def HestonPrice(kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,Trap,Lphi,Uphi,dphi,PutCall='call'):
    '''
    % Heston (1993) price of a European option.
    % Uses the original formulation by Heston
    % Heston parameters:
    %    kappa  = volatility mean reversion speed parameter
    %    theta  = volatility mean reversion level parameter
    %    lambda = risk parameter
    %    rho    = correlation between two Brownian motions
    %    sigma  = volatility of variance
    %    v0     = initial variance
    % Option features.
    %    PutCall = 'C'all or 'P'ut
    %    K = strike price
    %    S = spot price
    %    r = risk free rate
    %    q = dividend yield
    %    T = maturity
    % Integration features
    %    L = lower limit
    %    U = upper limit
    %    dphi = integration increment
    % 此处主要是用基于梯形法则的数值积分
    % example: HestonPrice('put',5.0,.05,0.0,-.8,.5,.5,100.0,100.0,.03,.02,.05,1,0.00001,50,0.001)
    '''
    phi=np.arange(Lphi,Uphi+.001,dphi)
    N=len(phi)
    int1=np.empty(N)
    int2=np.empty(N)
    if Trap == 1:
        Heston_Inte = HestonIntegrand_Trap
    else:
        Heston_Inte = HestonIntegrand_Orig
    int1=Heston_Inte(phi,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,0)
    int2=Heston_Inte(phi,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,1)
    I1=np.trapz(int1)*dphi
    I2=np.trapz(int2)*dphi
    P1=0.5+1.0/pi*I1
    P2=0.5+1.0/pi*I2
    if PutCall=='call':
        return S*np.exp(-q*T)*P1-K*np.exp(-r*T)*P2
    else:
        return K*np.exp(-r*T)*(1.0-P2)-S*np.exp(-q*T)*(1-P1)
v_HestonPrice=np.vectorize(HestonPrice)
def HestonCharacteristicFunction(phi,kappa,theta,lam,rho,sigma,tau,S,r,q,v0):
    a = kappa * theta
    u = -.5
    b = kappa + lam
    d = np.sqrt((rho * sigma * 1.0j * phi - b) ** 2 - sigma * sigma * (2 * u * 1.0j * phi - phi * phi))
    g = (b - rho * sigma * 1.0j * phi + d) / (b - rho * sigma * 1.0j * phi - d)
    G = (1 - g * np.exp(d * tau)) / (1 - g)
    C = (r - q) * 1.0j * phi * tau + a / sigma / sigma * ((b - rho * sigma * 1.0j * phi + d) * tau - 2 * np.log(G))
    D = (b - rho * sigma * 1.0j * phi + d) / sigma / sigma * (1 - np.exp(d * tau)) / (1 - g * np.exp(d * tau))
    return np.exp(C + D * v0 + 1.0j * phi * np.log(S))

def HestonInteConsol_Trap(phi,kappa,theta,lam,rho,sigma,tau,K,S,r,q,v0):
    '''
    % Returns the integrand for the risk neutral probabilities P1 and P2.
    % phi = integration variable
    XXXXXX Pnum = 1 or 2 (for the probabilities)XXXXX
    % Heston parameters:
    %    kappa  = volatility mean reversion speed parameter
    %    theta  = volatility mean reversion level parameter
    %    lambda = risk parameter
    %    rho    = correlation between two Brownian motions
    %    sigma  = volatility of variance
    %    v      = initial variance
    % Option features.
    %    PutCall = 'C'all or 'P'ut
    %    K = strike price
    %    S = spot price
    %    r = risk free rate
    %    q = dividend yield
    %    Trap = 1 "Little Trap" formulation
    %           0  Original Heston formulation
    '''
    a=kappa*theta
    u=0.5
    b=kappa+lam-rho*sigma
    d=np.sqrt((rho*sigma*1.0j*phi-b)**2-sigma*sigma*(2*u*1.0j*phi-phi*phi))
    g=(b-rho*sigma*1.0j*phi+d)/(b-rho*sigma*1.0j*phi-d)
    c=1.0/g
    D=(b-rho*sigma*1.0j*phi-d)/sigma/sigma*(1-np.exp(-d*tau))/(1-c*np.exp(-d*tau))
    G=(1-c*np.exp(-d*tau))/(1-c)
    C=(r-q)*1.0j*phi*tau+a/sigma/sigma*((b-rho*sigma*1.0j*phi-d)*tau-2*np.log(G))
    f1=np.exp(C+D*v0+1.0j*phi*log(S))

    u = -.5
    b = kappa + lam
    d = np.sqrt((rho * sigma * 1.0j * phi - b) ** 2 - sigma * sigma * (2 * u * 1.0j * phi - phi * phi))
    g = (b - rho * sigma * 1.0j * phi + d) / (b - rho * sigma * 1.0j * phi - d)
    c = 1.0 / g
    D = (b - rho * sigma * 1.0j * phi - d) / sigma / sigma * (1 - np.exp(-d * tau)) / (1 - c * np.exp(-d * tau))
    G = (1 - c * np.exp(-d * tau)) / (1 - c)
    C = (r - q) * 1.0j * phi * tau + a / sigma / sigma * ((b - rho * sigma * 1.0j * phi - d) * tau - 2 * np.log(G))
    f2= np.exp(C + D * v0 + 1.0j * phi * log(S))

    return (np.exp(-1.0j*phi*log(K))/1.0j/phi*(S*exp(-q*tau)*f1-K*exp(-r*tau)*f2)).real
def HestonInteConsol_Orig(phi,kappa,theta,lam,rho,sigma,tau,K,S,r,q,v0):
    '''
    % Returns the integrand for the risk neutral probabilities P1 and P2.
    % phi = integration variable
    XXXXXX Pnum = 1 or 2 (for the probabilities)XXXXX
    % Heston parameters:
    %    kappa  = volatility mean reversion speed parameter
    %    theta  = volatility mean reversion level parameter
    %    lambda = risk parameter
    %    rho    = correlation between two Brownian motions
    %    sigma  = volatility of variance
    %    v      = initial variance
    % Option features.
    %    PutCall = 'C'all or 'P'ut
    %    K = strike price
    %    S = spot price
    %    r = risk free rate
    %    q = dividend yield
    %    Trap = 1 "Little Trap" formulation
    %           0  Original Heston formulation
    '''
    a=kappa*theta
    u=0.5
    b=kappa+lam-rho*sigma
    d=np.sqrt((rho*sigma*1.0j*phi-b)**2-sigma*sigma*(2*u*1.0j*phi-phi*phi))
    g=(b-rho*sigma*1.0j*phi+d)/(b-rho*sigma*1.0j*phi-d)
    G=(1-g*np.exp(d*tau))/(1-g)
    C=(r-q)*1.0j*phi*tau+a/sigma/sigma*((b-rho*sigma*1.0j*phi+d)*tau-2*np.log(G))
    D=(b-rho*sigma*1.0j*phi+d)/sigma/sigma*(1-np.exp(d*tau))/(1-g*np.exp(d*tau))
    f1=np.exp(C+D*v0+1.0j*phi*log(S))

    u = -.5
    b = kappa + lam
    d = np.sqrt((rho * sigma * 1.0j * phi - b) ** 2 - sigma * sigma * (2 * u * 1.0j * phi - phi * phi))
    g = (b - rho * sigma * 1.0j * phi + d) / (b - rho * sigma * 1.0j * phi - d)
    G = (1 - g * np.exp(d * tau)) / (1 - g)
    C = (r - q) * 1.0j * phi * tau + a / sigma / sigma * ((b - rho * sigma * 1.0j * phi + d) * tau - 2 * np.log(G))
    D = (b - rho * sigma * 1.0j * phi + d) / sigma / sigma * (1 - np.exp(d * tau)) / (1 - g * np.exp(d * tau))
    f2 = np.exp(C + D * v0 + 1.0j * phi * log(S))
    return (np.exp(-1.0j*phi*log(K))/1.0j/phi*(S*exp(-q*tau)*f1-K*exp(-r*tau)*f2)).real
def HestonPriceConsol(kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,Trap,Lphi,Uphi,dphi,PutCall='call'):
    phi=np.arange(Lphi,Uphi,dphi)
    N=len(phi)
    inte=np.empty(N)
    if Trap == 1:
        Heston_Inte = HestonInteConsol_Trap
    else:
        Heston_Inte = HestonInteConsol_Orig
    inte=Heston_Inte(phi,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0)
    I=np.trapz(inte)*dphi
    if PutCall=='call':
        return S*np.exp(-q*T)*0.5-K*np.exp(-r*T)*0.5+I/pi
    else:
        return K*np.exp(-r*T)*0.5-S*np.exp(-q*T)*0.5+I/pi
v_HestonPriceConsol=np.vectorize(HestonPriceConsol)

def LocalVariance_FD(kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,Trap,Lphi,Uphi,dphi,dt,dK,PutCall='call'):
    '''obtain local volatility in the Heston model using the central difference approximation to the derivatives'''
    CT =v_HestonPriceConsol(kappa,theta,lam,rho,sigma,T+dt,K,S,r,q,v0,Trap,Lphi,Uphi,dphi,PutCall)
    CT_=v_HestonPriceConsol(kappa,theta,lam,rho,sigma,T-dt,K,S,r,q,v0,Trap,Lphi,Uphi,dphi,PutCall)
    dCdT=(CT-CT_)/2/dt
    CK=v_HestonPriceConsol(kappa,theta,lam,rho,sigma,T,K+dK,S,r,q,v0,Trap,Lphi,Uphi,dphi,PutCall)
    CK0=v_HestonPriceConsol(kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,Trap,Lphi,Uphi,dphi,PutCall)
    CK_=v_HestonPriceConsol(kappa,theta,lam,rho,sigma,T,K-dK,S,r,q,v0,Trap,Lphi,Uphi,dphi,PutCall)
    d2CdK2=(CK+CK_-2*CK0)/dK/dK
    return np.sqrt(2*dCdT/K/K/d2CdK2)
def dPjdT_Trap(phi,kappa,theta,lam,rho,sigma,T,K,S,r,v0,Pnum):
    u = np.array([0.5, -.5])[Pnum]
    b = np.array([kappa + lam - rho * sigma, kappa + lam])[Pnum]
    d = np.sqrt((rho * sigma * 1.0j * phi - b) ** 2 - sigma * sigma * (2 * u * 1.0j * phi - phi * phi))
    g = (b - rho * sigma * 1.0j * phi + d) / (b - rho * sigma * 1.0j * phi - d)
    edT=np.exp(-d*T)
    c=1.0/g
    dDdT=(b-rho*sigma*1.0j*phi-d)/sigma/sigma*d*edT*(1.0-c)/(1.0-c*edT)**2
    dCdT=r*1.0j*phi+kappa*theta/sigma/sigma*((b-rho*sigma*1.0j*phi-d)-2.0*c*d*edT/(1-c*edT))
    G=(1-c*edT)/(1-c)
    D=(b-rho*sigma*1.0j*phi-d)/sigma/sigma*(1-edT)/(1-c*edT)
    C=r*1.0j*phi*T+kappa*theta/sigma/sigma*((b-rho*sigma*1.0j*phi-d)*T-2.0*np.log(G))
    dfdT=np.exp(C+D*v0+1.0j*phi*log(S))*(dCdT+dDdT*v0)
    return (K**(-1.0j*phi)*dfdT/1.0j/phi).real
def dPjdT_Orig(phi,kappa,theta,lam,rho,sigma,T,K,S,r,v0,Pnum):
    u = np.array([0.5, -.5])[Pnum]
    b = np.array([kappa + lam - rho * sigma, kappa + lam])[Pnum]
    d = np.sqrt((rho * sigma * 1.0j * phi - b) ** 2 - sigma * sigma * (2 * u * 1.0j * phi - phi * phi))
    g = (b - rho * sigma * 1.0j * phi + d) / (b - rho * sigma * 1.0j * phi - d)
    edT=np.exp(d*T)
    dDdT=(b-rho*sigma*1.0j*phi-d)/sigma/sigma*d*edT*(g-1.0)/(1.0-g*edT)**2
    dCdT=r*1.0j*phi+kappa*theta/sigma/sigma*((b-rho*sigma*1.0j*phi+d)-2.0*g*d*edT/(1-g*edT))
    G=(1-g*edT)/(1-g)
    D=(b-rho*sigma*1.0j*phi+d)/sigma/sigma*(1-edT)/(1-g*edT)
    C=r*1.0j*phi*T+kappa*theta/sigma/sigma*((b-rho*sigma*1.0j*phi+d)*T-2.0*np.log(G))
    dfdT=np.exp(C+D*v0+1.0j*phi*log(S))*(dCdT+dDdT*v0)
    return (K**(-1.0j*phi)*dfdT/1.0j/phi).real
def d2P1dK2(phi,kappa,theta,lam,rho,sigma,T,K,S,r,v0):
    u=0.5
    b=kappa+lam-rho*sigma
    d = np.sqrt((rho * sigma * 1.0j * phi - b) ** 2 - sigma * sigma * (2 * u * 1.0j * phi - phi * phi))
    g = (b - rho * sigma * 1.0j * phi + d) / (b - rho * sigma * 1.0j * phi - d)
    edT = np.exp(d * T)
    G = (1 - g * edT) / (1 - g)
    D = (b - rho * sigma * 1.0j * phi + d) / sigma / sigma * (1 - edT) / (1 - g * edT)
    C = r * 1.0j * phi * T + kappa * theta / sigma / sigma * ((b - rho * sigma * 1.0j * phi + d) * T - 2.0 * np.log(G))
    return ((1.0j*phi+1)*K**(-1.0j*phi-2.0)*np.exp(C+D*v0+1.0j*phi*log(S))).real
def d2P2dK2(phi,kappa,theta,lam,rho,sigma,T,K,S,r,v0):
    u=-0.5
    b=kappa+lam
    d = np.sqrt((rho * sigma * 1.0j * phi - b) ** 2 - sigma * sigma * (2 * u * 1.0j * phi - phi * phi))
    g = (b - rho * sigma * 1.0j * phi + d) / (b - rho * sigma * 1.0j * phi - d)
    edT = np.exp(d * T)
    G = (1 - g * edT) / (1 - g)
    D = (b- rho * sigma * 1.0j * phi + d) / sigma / sigma * (1 - edT) / (1 - g * edT)
    C = r * 1.0j * phi * T + kappa * theta / sigma / sigma * ((b - rho * sigma * 1.0j * phi + d) * T - 2.0 * np.log(G))
    return ((1.0j*phi-1.0)*K**(-1.0j*phi-1.0)*np.exp(C+D*v0+1.0j*phi*log(S))).real
def LocalVariance_Analytic(kappa,theta,lam,rho,sigma,T,K,S,r,v0,Trap,Lphi,Uphi,dphi):
    phi=np.arange(Lphi,Uphi,dphi)
    if Trap==1:
        dPjdT=dPjdT_Trap
        HestonIntegrand=HestonIntegrand_Trap
    else:
        dPjdT=dPjdT_Orig
        HestonIntegrand = HestonIntegrand_Orig
    int1=dPjdT(phi,kappa,theta,lam,rho,sigma,T,K,S,r,v0,0)
    int2 = dPjdT(phi, kappa, theta, lam, rho, sigma, T, K, S, r, v0, 1)
    int3=HestonIntegrand(phi,kappa,theta,lam,rho,sigma,T,K,S,r,0.0,v0,1)
    int4=d2P1dK2(phi,kappa,theta,lam,rho,sigma,T,K,S,r,v0)
    int5 = d2P2dK2(phi, kappa, theta, lam, rho, sigma, T, K, S, r, v0)

    dP1dT=np.trapz(int1)*dphi/pi
    dP2dT = np.trapz(int2) * dphi / pi
    P2=0.5+np.trapz(int3)*dphi/pi
    dCdT=S*dP1dT-K*np.exp(-r*T)*(-r*P2+dP2dT)
    dP1dK2=np.trapz(int4)*dphi/pi
    TwodP2dK2=np.trapz(int5)*dphi/pi
    d2CdK2=S*dP1dK2-np.exp(-r*T)*TwodP2dK2
    return np.sqrt(2.0*dCdT/K/K/d2CdK2)
v_LocalVariance_Analytic=np.vectorize(LocalVariance_Analytic)
def LocalVariance_Approx(kappa,theta,rho,sigma,T,K,S,v0):
    kappa_=kappa-rho*sigma*0.5
    theta_=theta*kappa/kappa_
    xT=np.log(K/S)
    wT=(v0-theta)*(1-np.exp(-kappa*T))/kappa+theta*T
    vT=(v0-theta_)*np.exp(-kappa_*T)+theta_
    F1=(v0-theta)/(kappa_-kappa)
    E1=np.exp((kappa_-kappa)*T)-1.0
    F2=theta/kappa_
    E2=np.exp(kappa_*T)-1.0
    Integral=np.exp(-kappa_*T)*(F1*E1+F2*E2)
    return np.sqrt(vT+rho*sigma*xT/wT*Integral)

def AttariIntegrand_Trap(phi,kappa,theta,lam,rho,sigma,tau,K,S,r,q,v0):
    a=kappa*theta
    u=-0.5
    b=kappa+lam
    d = np.sqrt((rho * sigma * 1.0j * phi - b) ** 2 - sigma * sigma * (2 * u * 1.0j * phi - phi * phi))
    g = (b - rho * sigma * 1.0j * phi + d) / (b - rho * sigma * 1.0j * phi - d)
    c = 1.0 / g
    D = (b - rho * sigma * 1.0j * phi - d) / sigma / sigma * (1 - np.exp(-d * tau)) / (1 - c * np.exp(-d * tau))
    G = (1 - c * np.exp(-d * tau)) / (1 - c)
    C = (r - q) * 1.0j * phi * tau + a / sigma / sigma * ((b - rho * sigma * 1.0j * phi - d) * tau - 2.0 * np.log(G))
    f=np.exp(C+D*v0-1.0j*phi*r*tau)
    L=-r*tau+np.log(K/S)
    return ((f.real+f.imag/phi)*np.cos(L*phi)+(f.imag-f.real/phi)*np.sin(L*phi))/(1.0+phi*phi)
def AttariIntegrand_Orig(phi,kappa,theta,lam,rho,sigma,tau,K,S,r,q,v0):
    a=kappa*theta
    u=-0.5
    b=kappa+lam
    d = np.sqrt((rho * sigma * 1.0j * phi - b) ** 2 - sigma * sigma * (2 * u * 1.0j * phi - phi * phi))
    g = (b - rho * sigma * 1.0j * phi + d) / (b - rho * sigma * 1.0j * phi - d)

    D = (b - rho * sigma * 1.0j * phi + d) / sigma / sigma * (1 - np.exp(d * tau)) / (1 - g * np.exp(d * tau))
    G = (1 - g * np.exp(d * tau)) / (1 - g)
    C = (r - q) * 1.0j * phi * tau + a / sigma / sigma * ((b - rho * sigma * 1.0j * phi + d) * tau - 2.0 * np.log(G))
    f=np.exp(C+D*v0-1.0j*phi*r*tau)
    L=-r*tau+np.log(K/S)
    return ((f.real+f.imag/phi)*np.cos(L*phi)+(f.imag-f.real/phi)*np.sin(L*phi))/(1.0+phi*phi)
def AttariPrice(kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,Trap,Lphi,Uphi,dphi,PutCall='call'):
    phi=np.arange(Lphi,Uphi,dphi)
    if Trap==1:
        AttariIntegrand=AttariIntegrand_Trap
    else:
        AttariIntegrand=AttariIntegrand_Orig
    int1=AttariIntegrand(phi,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0)
    if PutCall=='call':
        return S*np.exp(-q*T)-K*np.exp(-r*T)*(0.5+np.trapz(int1)*dphi/pi)
    else:
        return K*np.exp(-r*T)*(0.5-np.trapz(int1)*dphi/pi)

def CarrMadanIntegrand(phi,alpha,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,PutCall):
    S *= exp(-q * T)
    iscall=1.0 if PutCall=='call' else -1.0
    return (np.exp(-1.0j * phi * log(K) - r * T) * HestonCharacteristicFunction(phi - (iscall*alpha + 1.0) * 1.0j, kappa, theta, lam,
                                                                        rho, sigma, T, S, r, 0.0,v0) / (
    alpha * alpha + iscall*alpha - phi * phi + 1.0j * (2.0*iscall * alpha + 1.0) * phi)).real
def CarrMadanIntegrandOTM(phi,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0):
    CF=HestonCharacteristicFunction(phi-1.0j,kappa,theta,lam,rho,sigma,T,S,r,q,v0)
    return (np.exp(-1.0j * phi * np.log(K) - r * T)*(S**(1.0j*phi+1.0)/(1.0+1.0j*phi)-np.exp(r*T)*S**(1.0j*phi+1.0)/(1.0j*phi)-CF/(phi*phi-1.0j*phi))).real
def CarrMadanPrice(alpha,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,Lphi,Uphi,dphi,PutCall='call'):
    phi=np.arange(Lphi,Uphi,dphi)
    int1=CarrMadanIntegrand(phi,alpha,kappa,theta,lam,rho,sigma,T,K,S,r,q,v0,PutCall)
    return exp(-alpha * log(K)) * np.trapz(int1) * dphi / pi


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
def Ame_option_ImDiff(St, K, r,q, sigma, T, M, N, otype='call'):
    '''

    :param St:
    :param K:
    :param r:
    :param q:
    :param sigma:
    :param T:
    :param M:  20
    :param N:  30
    :param otype:
    :return:
    '''
    dt = T / N
    dS=St/M
    S = dS * np.arange(1,2*M+1)
    o_value = 1.0 if otype is 'call' else -1.0
    payoff = np.where((S - K) * o_value > 0, (S - K) * o_value, 0)
    f=np.copy(payoff)
    A=np.zeros((2*M,2*M))
    for i in range(2*M):
        alpha=0.5*sigma*sigma*S[i]*S[i]*dt/dS/dS
        betha=0.5*(r-q)*S[i]*dt/dS
        Lph=-alpha+betha
        Dph=1+r*dt+2*alpha
        Uph=-alpha-betha
        if i==0:
            A[i,i]=Dph+2*Lph
            A[i,i+1]=Uph-Lph
        elif i<2*M-1:
            A[i,i-1]=Lph
            A[i,i]=Dph
            A[i,i+1]=Uph
        else:
            A[i,i-1]=Lph-Uph
            A[i,i]=Dph+2*Uph
    A_inv=np.linalg.pinv(A)
    for m in range(N-1):
        f=A_inv@f
        f=np.where(payoff>f,payoff,f)
    f = A_inv @ f
    return f[M-1]
def Euro_option_BSM(St,K,r,sigma,T,otype='call'):
    d1=(np.log(St / K) + (r + sigma *sigma / 2) * T) / (sigma*np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    if otype == 'call':
        return St * norm(0, 1).cdf(d1) - np.exp(-r * T) * K * norm(0, 1).cdf(d2)
    else:
        return np.exp(-r * T) * K * norm(0, 1).cdf(-d2) - St * norm(0, 1).cdf(-d1)
def vega_euro(St,K,r,sigma,T,otype='call'):
    d1 = (np.log(St / K) + (r + sigma *sigma / 2) * T) / (sigma * np.sqrt(T))
    return St * norm(0, 1).pdf(d1) * np.sqrt(T)
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
def Euro_option_trinomial(St,K,r,sigma,T,M,otype='call',method='CRR'):
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
            payoff[n] = (pu* payoff[n + 2] +pm*payoff[n+1]+ pd* payoff[n]) * disc
    return payff[0]
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
@jit
def Euro_option_ImDiff(St, K, r,q, sigma, T, M, N, otype='call'):
    '''

    :param St:
    :param K:
    :param r:
    :param q:
    :param sigma:
    :param T:
    :param M: 15
    :param N: 30
    :param otype:
    :return:
    '''
    dt = T / N
    dS=St/M
    S = dS * np.arange(1,2*M+1)
    o_value = 1.0 if otype is 'call' else -1.0
    payoff = np.where((S - K) * o_value > 0, (S - K) * o_value, 0)
    A=np.zeros((2*M,2*M))
    for i in range(2*M):
        alpha=0.5*sigma*sigma*S[i]*S[i]*dt/dS/dS
        betha=0.5*(r-q)*S[i]*dt/dS
        Lph=-alpha+betha
        Dph=1+r*dt+2*alpha
        Uph=-alpha-betha
        #A[i,i-1]=Lph-(i==2*M-1)*Uph
        #A[i,i]=Dph+2*Lph*(i==0)+2*Uph*(i==2*M-1)
        #A[i,i+1]=Uph-Lph*(i==0)
        if i==0:
            A[i,i]=Dph+2*Lph
            A[i,i+1]=Uph-Lph
        elif i<2*M-1:
            A[i,i-1]=Lph
            A[i,i]=Dph
            A[i,i+1]=Uph
        else:
            A[i,i-1]=Lph-Uph
            A[i,i]=Dph+2*Uph
    A_inv=np.linalg.pinv(A)
    for m in range(N):
        payoff=A_inv@payoff
    return payoff[M-1]
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
def imp_vol(Ct,St,K,r,T,sigma_est,otype):
    def diff(sigma):
        return Euro_option_BSM(St,K,r,sigma,T,otype)-Ct
    return fsolve(diff,sigma_est)[0]
def imp_vol_Formula(Ct,St,K,r,T,sigma_0,otype,vol=0.000001):
    iter_max_num=10
    sigma=sigma_0
    dC = Euro_option_BSM(St, K, r, sigma, T, otype) - Ct
    for i in range(iter_max_num):
        sigma=sigma-dC/vega_euro(St,K,r,sigma,T,otype)
        dC = Euro_option_BSM(St, K, r, sigma, T, otype) - Ct
        if np.percentile(abs(dC),10)<vol:
            break
    return sigma
class Norm:
    def __init__(self,mu,sigma):
        self.mu=mu
        self.sigma=sigma
    def cdf(self,d):
        d=(d-self.mu)/self.sigma
        a1=0.319381530
        a2=-0.356563782
        a3=1.781477937
        a4=-1.821255978
        a5=1.330274429
        gamma=0.2316419
        normalprime=math.exp(-d*d/2.0)/math.sqrt(2*math.pi)
        if d>0:
            k1 = 1.0 / (1.0 + gamma * d)
            return 1-normalprime*(a1*k1+a2*k1**2+a3*k1**3+a4*k1**4+a5*k1**5)
        else:
            k1 = 1.0 / (1.0 - gamma * d)
            return normalprime * (a1 * k1 + a2 * k1 ** 2 + a3 * k1 ** 3 + a4 * k1 ** 4 + a5 * k1 ** 5)

    def pdf(self, d):
        return np.exp(-(d - self.mu) ** 2 * 0.5 / self.sigma / self.sigma) / sqrt(2 * pi) / self.sigma
