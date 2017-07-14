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
from numba import jit
norm=scis.norm
class euro_option:
    def __init__(self,St=100,K=100,t=0,T=1,r=.05,sigma=.2,otype='call'):
        # (self,St,K,t,T,r,sigma,otype):#
        self.St = St
        self.K = K
        self.t = t
        self.T = T
        self.r = r
        self.sigma = sigma
        self.otype = otype

    def d1f(self):
        return (np.log(self.St / self.K) + (self.r + self.sigma ** 2 / 2) * (self.T-self.t)) / (self.sigma * np.sqrt(self.T-self.t))
    @property
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

#TODO 使用CRR运行太慢，待改进
    @jit
    def value_CRR(self, M=1000):
        '''
        Cox-Ross-Rubinstein European option valuation

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

        mu = np.arange(M + 1)
        mu = np.resize(mu, (M + 1, M + 1))
        md = np.transpose(mu)
        mu = u ** (mu - md)
        md = d ** md
        S = self.St * mu * md

        # inner value
        if self.otype == 'call':
            V = np.maximum(S - self.K, 0)
        else:
            V = np.maximum(self.K - S, 0)
        z = 0
        for t in range(M - 1, -1, -1):  # backwards iteration
            V[0:M - z, t] = (p * V[0:M - z, t + 1] + (1 - p) * V[1:M - z + 1, t + 1]) * df
            z += 1
        return V[0, 0]
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
            return -self.St*norm(0,1).pdf(d1)*self.sigma/2/np.sqrt(self.T-self.t)-self.r*self.K*np.exp(-self.r*(self.T-self.t))*norm(0,1).cdf(d2)
        else:
            return -self.St * norm(0, 1).pdf(d1) * self.sigma / 2 / np.sqrt(self.T - self.t) + self.r * self.K * np.exp(
                -self.r * (self.T - self.t)) * norm(0, 1).cdf(-d2)
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
    def __init__(self,St=100,K=100,t=0,T=1,r=.05,sigma=.2,otype='call'):
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
    def value_binomial_0(self, M=50):
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
            C_temp = (p[0:M - z, i + 1, :] * V[0:M - z, i + 1, :] + (1 - p[1:M - z + 1, i + 1, :]) * V[1:M - z + 1,i + 1, :]) * df[0:M - z,i + 1,:]
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
            C_temp = (p[0:M - z, i + 1, :] * V_temp[0:M - z, :] + (1 - p[1:M - z + 1, i + 1, :]) * V_temp[1:M - z + 1,
                                                                                                   :]) * df[0:M - z,
                                                                                                         i + 1, :]
            V_temp = np.where(h[0:M - z, i, :] > C_temp, h[0:M - z, i, :], C_temp)
            z += 1
        return V_temp
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
