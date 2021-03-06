import numpy as np
import statsmodels.api as sm
import pandas as pd

def OLS(Y,X,constant_included=True,print_res='params'):
    '''
    # Fit and summarize OLS model # 具有相同的解释变量
    spector_data = sm.datasets.spector.load()
    spector_data.exog = sm.add_constant(spector_data.exog, prepend=True)
    mod = sm.OLS(spector_data.endog, spector_data.exog)
    res = mod.fit()
    print(res.summary())
    def OLS_R2_adj(Y): # all Y share the X
        return OLS(Y,spector_data.exog,print_res='rsquared_adj')
    n=3000
    rand=np.random.randn(n)
    a=np.kron(spector_data.endog,rand).reshape(32,n)
    pd.DataFrame(a).apply(OLS_R2_adj) # 对于多个参数的func，可以使用pipe
    :param Y:
    :param X:
    :param constant_included:
    :param print_res: 'params'、'tvalues'、'pvalues'、'rsquared'、'rsquared_adj'
    :return:
    '''
    if constant_included is not True:
        X=sm.add_constant(X,prepend=True)
    res=sm.OLS(Y, X).fit()
    if print_res=='params':
        return res.params
    elif print_res=='tvalues':
        return res.tvalues
    elif print_res=='pvalues':
        return res.pvalues
    elif print_res=='rsquared':
        return res.rsquared
    elif print_res=='rsquared_adj':
        return res.rsquared_adj
    else:
        return res.params

def GLS(Y,X,sigma,constant_included=True,print_res='params'):
    '''
    def GLS_params(Y):
        return GLS(Y,data.exog,sigma,print_res='params')
    data = sm.datasets.longley.load()
    data.exog = sm.add_constant(data.exog)
    ols_resid = sm.OLS(data.endog, data.exog).fit().resid
    res_fit = sm.OLS(ols_resid[1:], ols_resid[:-1]).fit()
    rho = res_fit.params;rho
    from scipy.linalg import toeplitz
    order = toeplitz(np.arange(16));order
    sigma = rho**order
    n=3000
    rand=np.random.randn(n)
    a=np.kron(data.endog,rand).reshape(16,n)
    pd.DataFrame(a).apply(GLS_params)
    :param Y:
    :param X:
    :param constant_included:
    :param print_res: 'params'、'tvalues'、'pvalues'、'rsquared'、'rsquared_adj'
    :return:
    '''
    if constant_included is not True:
        X=sm.add_constant(X,prepend=True)
    res=sm.GLS(Y, X, sigma=sigma).fit()
    if print_res=='params':
        return res.params
    elif print_res=='tvalues':
        return res.tvalues
    elif print_res=='pvalues':
        return res.pvalues
    elif print_res=='rsquared':
        return res.rsquared
    elif print_res=='rsquared_adj':
        return res.rsquared_adj
    else:
        return res.params
