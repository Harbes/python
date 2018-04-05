import pandas as pd
import numpy as np
from scipy.stats import norm


# start = lambda eqd: eqd.index[0]
# end = lambda eqd: eqd.index[-1]
# days = lambda eqd: (eqd.index[-1] - eqd.index[0]).days
# trades_per_month = lambda eqd: eqd.groupby(
#     lambda x: (x.year, x.month)
# ).apply(lambda x: x[x != 0].count()).mean()
# profit = lambda eqd: eqd.sum()
# average = lambda eqd: eqd[eqd != 0].mean()
# average_gain = lambda eqd: eqd[eqd > 0].mean()
# average_loss = lambda eqd: eqd[eqd < 0].mean()
# winrate = lambda eqd: float(sum(eqd > 0)) / len(eqd)
# payoff = lambda eqd: eqd[eqd > 0].mean() / -eqd[eqd < 0].mean()
# pf = PF = lambda eqd: abs(eqd[eqd > 0].sum() / eqd[eqd < 0].sum())
# maxdd = lambda eqd: (eqd.cumsum().expanding().max() - eqd.cumsum()).max()
# rf = RF = lambda eqd: eqd.sum() / maxdd(eqd)
# trades = lambda eqd: len(eqd[eqd != 0])
# _days = lambda eqd: eqd.resample('D').sum().dropna() #意图是什么

def kwarges_Options():
    options={
        'periods_per_year':252.0,
        'MARR':None
    }
    return options
def Sharpe(df,Rf=0.0):
    ''' annualised Sharpe ratio '''
    d=df.sub(Rf,axis=0)
    return d.mean() / d.std() * (options['periods_per_year']**0.5)
def DownsideRisk(df,method='full',**kwargs):
    '''

    :param df:
    :param MARR:
        minimum acceptable rate of return
    :param method:
    :return:
    '''
    MARR=kwargs['MARR']
    if MARR is None:
        MARR=df.mean()
    d=df.sub(MARR,axis=1)
    if method == 'full':
        return np.sqrt((d**2.0)[d<0.0].sum()/len(d))
    elif method=='subset':
        return np.sqrt((d**2.0)[d<0.0].mean())
def Sortino(df,Rf=0.0,**kwargs):
    return df.sub(Rf,axis=0).mean()/DownsideRisk(df,**kwargs)
def VaR(df,alpha=0.05,method='history'):
    if method=='history':
        return df.quantile(q=alpha)
    elif method=='norm':
        return -norm.ppf(alpha)*df.std()-df.mean()
def ExpectedShortfall(df,alpha,method='history'):
    if method=='history':
        return df[df<=df.quantile(q=alpha)].mean()
    elif method=='norm':
        return df[df<=norm.ppf(alpha)*df.std()+df.mean()].mean()
def MaxDrawdown(df):
    df_CumRet=(df+1.0).cumprod()
    df_CumMax=df_CumRet.cummax()
    return (df_CumMax-df_CumRet).cummax()/df_CumMax*100.0
def ulcer(df):
    '''The ulcer index is a stock market risk measure or technical analysis indicator devised by Peter Martin in 1987'''
    d=(df+1.0).cumprod()
    d_cummax=d.cummax()
    return ((((d - d_cummax)/d_cummax) ** 2.0).mean()) ** 0.5


def UPI(df, Rf=0.0):
    return df.sub(Rf,axis=0) / ulcer(df)


#def mpi(eqd):
#    """ Modified UPI, with enumerator resampled to months (to be able to
#    compare short- to medium-term strategies with different trade frequencies. """
#    return eqd.resample('M').sum().mean() / ulcer(eqd)
#MPI = mpi


#def mcmdd(eqd, runs=100, quantile=0.99, array=False):
#    maxdds = [maxdd(eqd.take(np.random.permutation(len(eqd)))) for i in range(runs)]
#    if not array:
#        return pd.Series(maxdds).quantile(quantile)
#    else:
#        return maxdds


# def performance_summary(equity_diffs, quantile=0.99, precision=4):
#     def _format_out(v, precision=4):
#         if isinstance(v, dict):
#             return {k: _format_out(v) for k, v in list(v.items())}
#         if isinstance(v, (float, np.float)):
#             v = round(v, precision)
#         if isinstance(v, np.generic):
#             return np.asscalar(v)
#         return v
#
#     def force_quantile(series, q):
#         return sorted(series.values)[int(len(series) * q)]
#
#     eqd = equity_diffs[equity_diffs != 0]
#     if getattr(eqd.index, 'tz', None) is not None:
#         eqd = eqd.tz_convert(None)
#     if len(eqd) == 0:
#         return {}
#     hold = holding_periods(equity_diffs)
#
#     return _format_out({
#         'backtest': {
#             'from': str(start(eqd)),
#             'to': str(end(eqd)),
#             'days': days(eqd),
#             'trades': len(eqd),
#             },
#         'performance': {
#             'profit': eqd.sum(),
#             'averages': {
#                 'trade': average(eqd),
#                 'gain': average_gain(eqd),
#                 'loss': average_loss(eqd),
#                 },
#             'winrate': winrate(eqd),
#             'payoff': payoff(eqd),
#             'PF': PF(eqd),
#             'RF': RF(eqd),
#             },
#         'risk/return profile': {
#             'sharpe': sharpe(eqd),
#             'sortino': sortino(eqd),
#             'maxdd': maxdd(eqd),
#             'WCDD (monte-carlo {} quantile)'.format(quantile): mcmdd(eqd, quantile=quantile),
#             'UPI': UPI(eqd),
#             'MPI': MPI(eqd),
#             }
#         })
