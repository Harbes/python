import numpy as np
from scipy.stats import norm, chi2
from scipy.linalg import toeplitz
from math import sqrt
import statsmodels.tsa.stattools as tsa_tools
import matplotlib.pyplot as plt
from modules import simulation as sim



def autoCov(arr,lag):
    '''
    支持 1-D和2-D
    :param arr:
    :param lag:
    :return:
    '''
    if lag >= arr.shape[0] | lag<0:
        raise ValueError("Not enough observations to compute autocovariance, or lag<0")
    else:
        x1 = arr[lag:len(arr)]
        len_x1=len(x1)
        return np.sum((x1 - np.mean(arr)) * (arr[0:len_x1] - np.mean(arr))) / len_x1
def acf(arr,lag):
    return autoCov(arr,lag)/autoCov(arr,0)
def acfs(arr, max_lag):
    '''
    用于计算lag=1：max_lag的acf，同时计算置信区间？（没有用到置信水平，如何算出置信区间？？？）
    :param arr:
    :param max_lag:
    :return:
    '''
    acfs_ = np.array([acf(arr,lag) for lag in range(1, max_lag + 1)])
    acf_interval = np.ones(max_lag)
    acf_interval[1:] = 2 * acfs_[:-1] ** 2
    acf_interval = np.sqrt(np.cumsum(acf_interval) / len(arr))
    return np.r_[1, acfs_], np.r_[0, acf_interval]
def acfs_plot(arr, max_lag, alpha=0.05):
    '''
    画出给定最大滞后阶数的acf图形，包括置信区间
    :param arr:
    :param max_lag:
    :param alpha:
    :return:
    '''
    s = norm.ppf(1 - alpha / 2)
    acfs_, acf_interval = acfs(arr,max_lag)
    plt.figure(figsize=(15, 5))
    ax = plt.subplot(111)
    plt.stem(range(max_lag + 1), acfs_)
    plt.xlabel("Lag")
    plt.ylabel("Autocorrelation")
    plt.title("ACF")
    plt.plot(range(max_lag + 1), s * acf_interval, linestyle="-.", color="red")
    plt.axhline(0, linestyle='-', color='black')
    plt.plot(range(max_lag + 1), -s * acf_interval, linestyle="-.", color="red")
    plt.show()
def acf_test(arr,lag,test_type='individual'):
    '''
    test the acf 保持怀疑
    :param lag:
        lag=lag if test_type='individual'
        lag=max_lag(a list or tuple) if test_type='portmanteau'/'joint'
    :param test_type:
        'individual' or 'portmanteau/joint'
    :return:
        the test result
    '''

    if test_type == 'individual':
        if lag<=0:
            raise ValueError('lag must be greater than zero')
        t_ratio=acf(arr,lag)/sqrt((1+2*np.sum([acf(arr,i)**2 for i in range(lag-1)]))/len(arr)) #?
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

def pacfs(arr,max_lag,method='ols'):
    if method=='Yule-Walker':
        pacfs = [1]
        for j in range(1,max_lag+1):
            gamma_=[autoCov(arr,i) for i in range(j+1)]
            R=toeplitz(gamma_[:-1])
            pacfs.append(np.linalg.solve(R,gamma_[1:])[-1])
        return np.array(pacfs)
    else:
        return  np.r_[1,[pacf_ols(arr,i) for i in range(1,max_lag+1)]]
def pacf_ols(arr,lag):
    if lag<=0:
        raise ValueError('lag不能小于1')
    tmp=np.array([arr[lag-i:-i] for i in range(1,lag+1)])
    y=arr[lag:]
    x=np.append([np.ones(len(y))],tmp,axis=0)
    return y@x.T@np.linalg.pinv(x@x.T)[lag]
def pacf_plot(arr,max_lag,alpha=0.05,method='ols'):
    s = norm.ppf(1 - alpha / 2)
    pacfs_=pacfs(arr,max_lag,method=method)
    plt.figure(figsize=(15, 5))
    ax = plt.subplot(111)
    plt.stem(range(max_lag+1), pacfs_)
    plt.xlabel("Lag")
    plt.ylabel("Partial Autocorrelation")
    plt.title("PACF")
    plt.axhline(s / np.sqrt(len(arr)), linestyle="-.", color="red")
    plt.axhline(0, linestyle='-', color='black')
    plt.axhline(- s / np.sqrt(len(arr)), linestyle="-.", color="red")
    plt.show()
def acf_pacf_plot(arr,max_lag,alpha=0.05,pacf_method='ols'):
    s = norm.ppf(1 - alpha / 2)
    acfs_, acf_interval = acfs(arr,max_lag)
    pacfs_=pacfs(arr,max_lag,method=pacf_method)
    fig=plt.figure(figsize=(15, 5))
    ax1= fig.add_subplot(211)
    ax1.stem(range(max_lag+1), acfs_)
    plt.ylabel("ACF")
    plt.plot(range( max_lag+1),s *acf_interval, linestyle="-.", color="red")
    plt.axhline(0,linestyle='-',color='black')
    plt.plot(range(max_lag+1),-s*acf_interval, linestyle="-.", color="red")

    ax2 = fig.add_subplot(212)
    ax2.stem(range(max_lag+1), pacfs_)
    plt.xlabel("Lag")
    plt.ylabel("PACF")
    plt.axhline(s / np.sqrt(len(arr)), linestyle="-.", color="red")
    plt.axhline(0, linestyle='-', color='black')
    plt.axhline(- s / np.sqrt(len(arr)), linestyle="-.", color="red")
    plt.show()
acf_pacf_plot(a,30)