# estimate the (effective) bid-ask spread, and compare it with illiquidity measure(Amihud,2002), turnover ...
# paper: 【RFS.2017】A Simple Estimation of Bid-Ask Spreads from Daily Close, High, and Low Prices
#

import pandas as pd
from pandas import DataFrame
import numpy as np
from pandas.tseries.offsets import YearEnd
from datetime import datetime


def import_data():
    global data_path, amount, size_free, size_tot, adj_open, adj_close, low, high, close
    data_path = '/Users/harbes/data/xccdata'
    # data_path='F:/data/xccdata'
    data_set = 'PV_datetime'
    Vars = ['size_tot', 'size_free', 'adj_close', 'adj_open', 'amount', 'low', 'high', 'clsprc']
    data = pd.read_pickle(data_path + '/' + data_set)[Vars]
    amount = data['amount'].unstack() * 0.1  # 保持单位与size一致
    size_free = data['size_free'].unstack()
    size_tot = data['size_tot'].unstack()
    adj_open = data['adj_open'].unstack()
    adj_close = data['adj_close'].unstack()
    close = data['clsprc'].unstack()
    low = data['low'].unstack()
    high = data['high'].unstack()


if __name__ is '__main__':
    import_data()
    spread = 2 * np.sqrt(
        (np.log(close) - np.log(high * low) * 0.5) * (np.log(close) - np.log((high * low).shift(-1)) * 0.5))
    illiquidity = np.abs((adj_close - adj_open) / adj_open) * 10000 / amount
    turnover = amount / size_free
    spread.corrwith(illiquidity).mean()
    spread.corrwith(turnover).mean()  # spread与turnover正相关？？？
    turnover.corrwith(illiquidity).mean()
    spread.T.corrwith(illiquidity.T).mean()
    spread.T.corrwith(turnover.T).mean()
    turnover.T.corrwith(illiquidity.T).mean()
