import numpy as np
from scipy.stats import norm, chi2
from scipy.linalg import toeplitz
from math import sqrt
import statsmodels.tsa.api as tsa
import statsmodels.api as sm
import matplotlib.pyplot as plt
from modules import simulation as sim


def acovfs(arr):
    '''
    计算auto-cov；
    :param arr: 1-D array
    :return: lag=0~len(arr)
    '''
    return tsa.stattools.acovf(arr)
def acovf(arr,lag):
    '''
    支持 1-D和2-D；unbiased
    :param arr:
    :param lag:
    :return:
    '''
    return (arr[lag:len(arr)] - np.mean(arr)) @ (arr[0:len(arr)-lag] - np.mean(arr))/(len(arr)-lag)
def acfs(arr,lags,LBq=False,alpha=None):
    return tsa.stattools.acf(arr,nlags=lags,qstat=LBq,alpha=alpha)
def pacfs(arr,lags,alpha=None):
    return tsa.stattools.pacf(arr,nlags=lags,method='ols',alpha=alpha)
def plot_acfs_pacfs(arr,lags,method='statsmodels'):
    if method=='statsmodels':
        fig = plt.figure(figsize=(12, 8))
        ax1 = fig.add_subplot(211)
        fig = sm.graphics.tsa.plot_acf(arr, lags=lags, ax=ax1)
        ax2 = fig.add_subplot(212)
        fig = sm.graphics.tsa.plot_pacf(arr, lags=lags, ax=ax2)
        plt.show()
    else:
        acf_, conf_acf = acfs(arr, lags, alpha=0.05)
        pacf_, conf_pacf = pacfs(arr, lags, alpha=0.05)
        conf_acf -=acf_[:,None]
        conf_pacf -=pacf_[:,None]
        fig = plt.figure(figsize=(15, 5))
        ax1 = fig.add_subplot(211)
        ax1.stem(range(1, len(acf_)), acf_[1:])
        plt.ylabel("ACF")
        plt.plot(range(1, len(acf_)), conf_acf[1:, 1], linestyle="-.", color="red")
        plt.axhline(0, linestyle='-', color='black')
        plt.plot(range(1, len(acf_)), conf_acf[1:, 0], linestyle="-.", color="red")

        ax2 = fig.add_subplot(212)
        ax2.stem(range(1, len(pacf_)), pacf_[1:])
        plt.xlabel("Lag")
        plt.ylabel("PACF")
        plt.plot(range(1, len(pacf_)), conf_pacf[1:, 1], linestyle="-.", color="red")
        plt.axhline(0, linestyle='-', color='black')
        plt.plot(range(1, len(pacf_)), conf_pacf[1:, 0], linestyle="-.", color="red")
        plt.show()

def ccovfs(x,y):
    '''
    计算交叉协方差
    :param x:
    :param y:
    :return:
    '''
    return tsa.stattools.ccovf(x,y)
def ccfs(x,y):
    '''
    cross correlation
    :param x:
    :param y:
    :return:
    '''
    return tsa.stattools.ccf(x,y)
def periodogram():
    # TODO 待学习 periodogram
    pass
def adfuller(x):
    '''
    Augmented Dickey-Fuller unit root test
    :param x:
    :return:
    '''
    adf,pvalue,usedlag,nobs,cvalues,_=tsa.stattools.adfuller(x) # 可以改变 maxlag参数
    print('            adf:',adf)
    print('         pvalue:', pvalue)
    print('critical values:','\n',
           '               ',' 1%:', cvalues['1%'])
    print('                ',' 5%:',cvalues['5%'])
    print('                ','10%:', cvalues['10%'])
    print('        usedlag:', usedlag)
    print('           nobs:', nobs)
def GrangerCausalityTests(arr2,maxlag):
    '''

    :param arr2: 2D-array(nobs.2)
    :param maxlag:
    :return:
    '''
    return tsa.stattools.grangercausalitytests(arr2,maxlag)
def levinson_durbin(arr,nlags=10):
    # TODO
    return tsa.stattools.levinson_durbin(arr,nlags=nlags)
def arma_order_select_ic():
    # TODO
    pass
def ARMA_fit(endogenous,order,exog=None,trend='c'):
    '''
    例如：ARMA_fit(y,[3,0,1])
    使用results.summary()打印结果，或者使用results.params或者 .resid/.pvalues
    :param endogenous:
    :param order:
    :param trend:
    :return:
    '''
    return tsa.ARMA(endogenous,order,exog=exog).fit(trend=trend)
def ARIMA_fit(endog,order,exog=None,trend='c',method='css-mle'):
    return tsa.ARIMA(endog,order,exog=exog).fit(trend=trend,method=method)
def KalmanFilter():
    # TODO
    pass
def VecAutoReg_fit(endog,maxlag=None,method='ols',trend='c'):
    return tsa.VAR(endog).fit(maxlags=maxlag,method=method,trend=trend)
def DynamicVAR():
    #TODO
    pass
def SVAR():
    #TODO
    pass



matfn = 'E:/pythonD/TSA_MATLAB/sims_data.mat'  # the path of .mat data
data = sio.loadmat(matfn)
r=data['ytdata'][:,0]
mpo=np.log(data['ytdata'][:,[3,4,5]])
y=np.c_[r[12:],(mpo[12:]-mpo[:-12])*100]
VecAutoReg_fit(y,maxlag=2).summary()
