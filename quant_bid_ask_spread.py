# estimate the (effective) bid-ask spread, and compare it with illiquidity measure(Amihud,2002), turnover ...
# paper: 【RFS.2017】A Simple Estimation of Bid-Ask Spreads from Daily Close, High, and Low Prices
# Critical Review：中国股票的换手率太高(流动性较好)，而上述方法更适合非流动性股票（可以尝试按size-liquidity进行分组，再分别计算不同测度之间的相关性）

import pandas as pd
from pandas import DataFrame
import numpy as np
from pandas.tseries.offsets import YearEnd
from datetime import datetime


def import_data():
    global data_path, amount, size_free, size_tot, adj_open, adj_close, adj_high, adj_low, low, high, close
    data_path = '/Users/harbes/data/xccdata'
    # data_path='F:/data/xccdata'
    data_set = 'PV_datetime'
    Vars = ['size_tot', 'size_free', 'adj_close', 'adj_open', 'adj_high', 'adj_low', 'amount', 'low', 'high', 'clsprc']
    data = pd.read_pickle(data_path + '/' + data_set)[Vars]
    amount = data['amount'].unstack() * 0.1  # 保持单位与size一致
    size_free = data['size_free'].unstack()
    size_tot = data['size_tot'].unstack()
    adj_open = data['adj_open'].unstack()
    adj_close = data['adj_close'].unstack()
    adj_high = data['adj_high'].unstack()
    adj_low = data['adj_low'].unstack()
    close = data['clsprc'].unstack()
    low = data['low'].unstack()
    high = data['high'].unstack()


def calculate_liquidity_measures():
    global spread, illiquidity, turnover
    # spread = 2*np.sqrt(np.maximum((np.log(close) - np.log(high * low) * 0.5) * (np.log(close) - np.log((high * low).shift(-1)) * 0.5),0))
    spread = 2 * np.sqrt(np.maximum((np.log(adj_close) - np.log(adj_high * adj_low) * 0.5) * (
    np.log(adj_close) - np.log((adj_high * adj_low).shift(-1)) * 0.5), 0))
    illiquidity = np.abs((adj_close - adj_open) / adj_open) * 10000 / amount
    turnover = amount / size_free


def resample_data():
    global spread, illiquidity, turnover
    key = lambda x: x.year * 100 + x.month
    spread = spread.groupby(key).mean() * 100
    illiquidity = illiquidity.groupby(key).mean()
    turnover = turnover.groupby(key).sum()
if __name__ is '__main__':
    import_data()
    calculate_liquidity_measures()
    resample_data()
    spread.corrwith(illiquidity).mean()
    spread.corrwith(turnover).mean()  # spread与turnover正相关
    turnover.corrwith(illiquidity).mean()
    spread.T.corrwith(illiquidity.T).mean()
    spread.T.corrwith(turnover.T).mean()
    turnover.T.corrwith(illiquidity.T).mean()
    spread.corrwith(1 / close).mean()
    spread.T.corrwith(
        1 / close.T).mean()  # estimated-spread has a higher cross-sectional correlation with the reciprocal of close pric
