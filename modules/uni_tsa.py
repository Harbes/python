import numpy as np
from scipy.stats import norm, chi2
from scipy.linalg import toeplitz
from math import sqrt
import matplotlib.pyplot as plt
from numba import jit
class uni_tsa:
    def __init__(self,arr):
        if type(arr) is not np.ndarray:
            self.arr=np.array(arr)
        else:
            self.arr=arr
    def unit_test(self,lag=1,ct='c'):
        from functools import reduce
        delta_x=self.arr[1:]-self.arr[:-1]
        tmp=np.array([delta_x[lag-i:-i] for i in range(1,lag+1)])
        y = delta_x[lag:]
        if ct == 'c':
            if lag>0:
                x=reduce(lambda x, y: np.append(x, y, axis=0), ([np.ones(len(delta_x)-lag)],[self.arr[lag:-1]],tmp))
            else:
                x = reduce(lambda x, y: np.append(x, y, axis=0),
                           ([np.ones(len(delta_x) - lag)], [self.arr[lag:-1]]))
            beta = y @ x.T @ np.linalg.pinv(x @ x.T)
            sigma = np.std(y - beta @ x)
            return beta[1]/sigma/sqrt(np.linalg.pinv(x @ x.T)[1,1])
        elif ct =='ct':
            if lag>0:
                x = reduce(lambda x, y: np.append(x, y, axis=0), ([np.ones(len(delta_x) - lag)],\
                                        [np.array(range(1,len(delta_x) - lag + 1))], [self.arr[lag:-1]], tmp))
            else:
                x = reduce(lambda x, y: np.append(x, y, axis=0), ([np.ones(len(delta_x) - lag)],\
                                                                  [np.array(range(1, len(delta_x) - lag + 1))],\
                                                                  [self.arr[lag:-1]]))
            beta=y @ x.T @ np.linalg.pinv(x @ x.T)
            sigma=np.std(y-beta@x)
            return beta[2]/sigma/sqrt(np.linalg.pinv(x @ x.T)[2,2])
        elif ct=='ctt':
            if lag>0:
                x = reduce(lambda x, y: np.append(x, y, axis=0), ([np.ones(len(delta_x) - lag)], \
                    [np.array(range(1, len(delta_x) - lag + 1))],[(np.array(range(1,len(delta_x) - lag + 1)))**2],\
                                                                  [self.arr[lag:-1]], tmp))
            else:
                x = reduce(lambda x, y: np.append(x, y, axis=0), ([np.ones(len(delta_x) - lag)], \
                    [np.array(range(1, len(delta_x) - lag + 1))],[(np.array(range(1, len(delta_x) - lag + 1))) ** 2],\
                                                                  [self.arr[lag:-1]]))
            beta=y @ x.T @ np.linalg.pinv(x @ x.T)
            sigma=np.std(y-beta@x)
            return beta[3]/sigma/sqrt(np.linalg.pinv(x @ x.T)[3,3])
        else:
            if lag>0:
                x=np.append([self.arr[lag:-1]],tmp, axis=0)
            else:
                x = self.arr[lag:-1]
            beta = y @ x.T @ np.linalg.pinv(x @ x.T)
            sigma = np.std(y - beta @ x)
            return beta[0]/sigma/sqrt(np.linalg.pinv(x @ x.T)[0,0])

    def cov(self,lag=0):
        if lag >= self.arr.shape[0]:
            raise ValueError("Not enough observations to compute autocovariance")
        else:
            x1 = self.arr[lag:len(self.arr)]
            return np.sum((x1 - np.mean(self.arr)) * (self.arr[0:len(self.arr) - lag] - np.mean(self.arr))) / len(x1)
    def acf(self,lag=0):
        return self.cov(lag=lag)/self.cov(lag=0)
    def _acfs(self,max_lag):
        acfs = np.array([self.acf(lag) for lag in range(1,max_lag+1)])
        acf_interval = np.ones(max_lag)
        acf_interval[1:] = 2 * acfs[:-1] ** 2
        acf_interval = np.sqrt(np.cumsum(acf_interval) / len(self.arr))
        return np.r_[1,acfs],np.r_[0,acf_interval]

    def acf_plot(self, max_lag=24,alpha=0.05):
        s=norm.ppf(1-alpha/2)
        acfs,acf_interval=self._acfs(max_lag)
        plt.figure(figsize=(15, 5))
        ax = plt.subplot(111)
        plt.stem(range(max_lag+1),acfs)
        plt.xlabel("Lag")
        plt.ylabel("Autocorrelation")
        plt.title("ACF")
        plt.plot(range(max_lag+1), s *acf_interval, linestyle="-.", color="red")
        plt.axhline(0, linestyle='-', color='black')
        plt.plot(range(max_lag+1), -s *acf_interval, linestyle="-.", color="red")
        plt.show()
    def acf_test(self,lag=None,test_type='individual'):
        '''
        test the acf
        :param lag: 
            lag=lag if test_type='individual'
            lag=max_lag(a list or tuple) if test_type='portmanteau'/'joint'
        :param test_type: 
            'individual' or 'portmanteau/joint'
        :return: 
            the test result
        '''

        if test_type == 'individual':
            if lag is None:
                lag=1
            elif lag<=0:
                raise ValueError('lag must be greater than zero')
            t_ratio=self.acf(lag=lag)/sqrt((1+2*np.sum([self.acf(lag=i)**2 for i in range(lag-1)]))/len(self.arr))
            p_value=2*norm(0,1).cdf(-np.abs(t_ratio))
            print('===========================')
            print('Null hypotheses is rho_l=0 ')
            print('===========================')
            print('   t statistics is %12.3f' %t_ratio)
            print('   p-value is %17.3f' %p_value)
        else:
            if lag is None:
                max_lags=[5,10]
            elif (np.array(lag)<=0).any():
                raise ValueError('lag must be greater than zero')
            else:
                max_lags=lag
            T=len(self.arr)
            print('=============================================')
            print('Test for autocorrelation: Ljung and Box(1978)')
            print('Null hypotheses is rho_1=rho_2=...=0')
            print('=============================================')
            for m in range(len(max_lags)):
                Q = T * (T + 2) * np.sum([self.acf(lag=i) ** 2 / (T - i) for i in range(1, max_lags[m])])
                p_value = 1-chi2.cdf(Q, max_lags[m])
                print('Q({}):'.format(max_lags[m]))
                print('   chi_square statistics is %12.3f' % Q)
                print('   p-value is %26.3f' % p_value)

    def pacf(self,max_lag=40,method='Yule-Walker'):

        if method=='Yule-Walker':
            pacfs = [1]
            for j in range(1,max_lag+1):
                gamma=[self.cov(lag=i) for i in range(j+1)]
                R=toeplitz(gamma[:-1])
                pacfs.append(np.linalg.solve(R, gamma[1:])[-1])
            return pacfs
        else:
            return  np.r_[1,[self.pacf_ols(lag=i) for i in range(1,max_lag+1)]]
    def pacf_ols(self,lag=1):
        tmp=np.array([self.arr[lag-i:-i] for i in range(1,lag+1)])
        y=self.arr[lag:]
        x=np.append([np.ones(len(y))],tmp,axis=0)
        return y@x.T@np.linalg.pinv(x@x.T)[lag]
    def pacf_plot(self,max_lag=24,alpha=0.05,method='Yule-Walker'):
        s = norm.ppf(1 - alpha / 2)
        pacfs=self.pacf(max_lag=max_lag,method=method)
        plt.figure(figsize=(15, 5))
        ax = plt.subplot(111)
        plt.stem(range(max_lag+1), pacfs)
        plt.xlabel("Lag")
        plt.ylabel("Partial Autocorrelation")
        plt.title("PACF")
        plt.axhline(s / np.sqrt(len(self.arr)), linestyle="-.", color="red")
        plt.axhline(0, linestyle='-', color='black')
        plt.axhline(- s / np.sqrt(len(self.arr)), linestyle="-.", color="red")
        plt.show()
    def acf_pacf_plot(self,max_lag=24,alpha=0.05,pacf_method='Yule-Walker'):
        s = norm.ppf(1 - alpha / 2)
        acfs, acf_interval = self._acfs(max_lag)
        pacfs=self.pacf(max_lag=max_lag,method=pacf_method)
        fig=plt.figure(figsize=(15, 5))
        ax1= fig.add_subplot(211)
        ax1.stem(range(max_lag+1), acfs)
        plt.ylabel("ACF")
        plt.plot(range( max_lag+1),s *acf_interval, linestyle="-.", color="red")
        plt.axhline(0,linestyle='-',color='black')
        plt.plot(range(max_lag+1),-s*acf_interval, linestyle="-.", color="red")

        ax2 = fig.add_subplot(212)
        ax2.stem(range(max_lag+1), pacfs)
        plt.xlabel("Lag")
        plt.ylabel("PACF")
        plt.axhline(s / np.sqrt(len(self.arr)), linestyle="-.", color="red")
        plt.axhline(0, linestyle='-', color='black')
        plt.axhline(- s / np.sqrt(len(self.arr)), linestyle="-.", color="red")
        plt.show()
    def model(self):
        # TODO 时间序列模型参数估计，主要针对 ARIMA 模型
        pass
    def garch(self,p=1,o=0,q=1,beta_init=[.1,.08,.9]):
        # TODO 需要提高GARCH运行速度
        from scipy.optimize import fmin
        from numba import jit
        #,fmin_powell,fmin_cg,fmin_bfgs,fmin_ncg,fmin_l_bfgs_b,fmin_cobyla
        rtn=self.arr-np.mean(self.arr)
        @jit
        def neg_log_likelihood(beta):
            w=beta[0]
            a=beta[1]
            b=beta[2]
            rtn2=rtn**2
            T=len(rtn2)
            ht=np.empty(T)
            ht[0]=sum(rtn2)/T
            for i in range(1,T):
                ht[i]=w+a*rtn2[i-1]+b*ht[i-1]
            sqrt_ht=np.sqrt(ht)
            x=rtn/sqrt_ht
            punishment=10000 if a+b>=1 else 0
            return (T-1)/2*np.log(2*np.pi)+.5*np.sum(np.log(ht[1:])+x[1:]**2)+punishment
        def constraint(beta):
            return 1-(beta[1]+beta[2])
        theta_estimate=fmin(func=neg_log_likelihood,x0=beta_init)
        #theta_estimate = fmin_cobyla(neg_log_likelihood, x0=beta_init,cons=[constraint,])

        return theta_estimate

    def Granger(self,arr2):
        # TODO Granger因果检验
        pass

class simulation:
    def __init__(self,phi=None,theta=None,const=0,y_init=None,dist='normal',dist_sigma=1,M=1000,seed=None):
        self.phi=phi
        self.theta=theta
        self.const=const
        self.y_init=y_init
        self.dist=dist
        self.dist_sigma=dist_sigma
        self.M = M
        self.seed=seed
        if self.theta is None:
            self.theta=0
            self.q=0
        else:
            self.theta=np.array(self.theta)
            self.q=(self.theta).size
        if self.phi is None:
            self.phi=0
            self.p=0
        else:
            self.phi=np.array(self.phi)
            self.p=(self.phi).size
        self.pq=np.maximum(self.p,self.q)
    @property
    def generate_random(self):
        if self.seed:
            np.random.seed(self.seed)
        if self.dist=='normal':
            epsilon=np.random.normal(0,scale=self.dist_sigma,size=self.M+self.pq)
        return epsilon
    @property
    def MA(self):
        epsilon=self.generate_random
        e_array=(np.vstack([epsilon[self.pq-i-1:self.M+self.pq-1-i] for i in range(self.q)])).T
        return self.const+(self.theta*e_array).sum(axis=1)+epsilon[self.pq:]
    @property
    def AR(self):
        epsilon = self.generate_random
        y_t=np.empty(self.M+self.pq)
        if self.y_init is None:
            y_t[:self.pq] = 0
        else:
            y_t[self.pq-self.p:self.pq] = np.array(self.y_init)
        for i in range(self.pq, self.M + self.pq):
            y_t[i] = self.const + (self.phi * y_t[i - self.p:i]).sum() + epsilon[i]
        return y_t[self.pq:]
    @property
    def ARMA(self):
        epsilon = self.generate_random
        e_array = (np.vstack([epsilon[self.pq - 1 - i:self.M + self.pq - 1 - i] for i in range(self.q)])).T
        y_t=np.empty(self.M+self.pq)
        if self.y_init is None:
            y_t[:self.pq] = 0
        else:
            y_t[self.pq-self.p:self.pq] = np.array(self.y_init)
        for i in range(self.pq, self.M + self.pq):
            y_t[i] = self.const + (self.phi * y_t[i - self.p:i]).sum() +epsilon[i]+(self.theta*e_array[i-self.pq]).sum()
        return y_t[self.pq:]





    





