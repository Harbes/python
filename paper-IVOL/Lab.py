# used to accomplish the IVOL in 《empirical asset pricing》
# 1. calculate the market beta
# 2. calculate the returns of two mimicking portfolios--- size & BM

import pandas as pd
import numpy as np
from time import time
from pandas.tseries.offsets import MonthEnd
from scipy.stats import skew,kurtosis
#import matplotlib.pyplot as plt

def import_data():
    global data_path,pv,open_price,close_price,index_ret,rtn
    data_path = 'E:/data/NewData/'  # '/Users/harbes/data/NewData/' #
    pv = pd.read_pickle(data_path + 'PV_datetime')[['adj_close', 'adj_open', 'size_tot']]
    close_price = pv['adj_close'].unstack()
    # close_price.index=pd.to_datetime(close_price.index.astype(int).astype(str),format='%Y%m%d')
    open_price = pv['adj_open'].unstack()
    # open_price.index=pd.to_datetime(open_price.index.astype(int).astype(str),format='%Y%m%d')
    index_ret = pd.read_pickle(data_path + 'index_ret').set_index(['index_code', 'trddt'])['pctchange'].loc['000016.SH']
    index_ret.index = pd.to_datetime(index_ret.index.astype(int).astype(str), format='%Y%m%d')
    rtn = (close_price - open_price) / open_price * 100
    rtn['index'] = index_ret  # ['pctchange']
    # 注意，rtn的最后一列是月份序号
    rtn['month'] = (rtn.index.year - rtn.index[0].year) * 12 + rtn.index.month
    index_ret = pd.DataFrame({'pctchange': index_ret,
                              'month': (index_ret.index.year - index_ret.index[0].year) * 12 + index_ret.index.month})


def cal_beta(periods,save_data=False):
    '''
    calculate the market beta
    :param periods: int{measure periods}
    :param save_data:
    :return:
    '''
    beta = pd.DataFrame(index=pd.date_range('20050101', '20180301', freq='M'), columns=rtn.columns)
    t0 = time()
    for i in range(periods, 159):
        beta.iloc[i - 1] = rtn.iloc[:, :-1][(i - periods+1 <= rtn['month']) & (rtn['month'] <= i)].cov(min_periods=20).iloc[-1,
                           :-1] / \
                           np.var(index_ret['pctchange'][(index_ret['month'] >= i - periods+1) & (index_ret['month'] <= i)])
    print(time() - t0)
    if save_data:
        beta.to_pickle(data_path + 'beta_daily_'+str(periods)+'M')
    else:
        return beta

def import_book():
    '''
    注意：book是由assets-liability得来的，数据中有负数存在
    monthly
    :return:
    '''
    BS = pd.read_pickle(data_path + '/BS')[['fin_year', 'stkcd', 'tot_assets', 'tot_liab']].set_index(
        ['fin_year', 'stkcd']).sort_index()
    book = BS['tot_assets'] - BS['tot_liab']
    book = book.drop(book.index[book.index.duplicated(keep='last')]).unstack() * 1e-8
    book.index = pd.to_datetime(book.index.astype(str), format='%Y%m%d')
    book.index.name = 'trddt'
    book = book.resample('M').first().ffill()
    return book[book>0]
def import_size():
    '''
    monthly
    :return:
    '''
    GroupBy = lambda x: x.year * 100 + x.month
    size = pv['size_tot'].unstack().groupby(GroupBy).last()*1e-4
    size.index = pd.to_datetime(size.index.astype(str), format='%Y%m') + MonthEnd()
    return size
def cal_size_BM(save_data=False):
    book=import_book()
    size=import_size()
    BM=book/size
    if save_data:
        size.to_pickle(data_path+'size_monthly')
        BM.to_pickle(data_path+'BM_monthly')
    else:
        return size,BM
def cal_mom(periods,save_data=False):
    GroupBy = lambda x: x.year * 100 + x.month
    cls=close_price.groupby(GroupBy).nth(-1)
    #cls.index=pd.to_datetime(cls.index.astype(str),format='%Y%m')
    opn=open_price.groupby(GroupBy).nth(0)
    mom=(cls.shift(1)-opn.shift(periods))/opn.shift(periods)*100
    mom.index = pd.to_datetime(mom.index.astype(str), format='%Y%m')+MonthEnd()
    if save_data:
        mom.to_pickle(data_path+'Mom_'+str(periods)+'M')
    else:
        return mom
def cal_rev(save_data=False):
    GroupBy = lambda x: x.year * 100 + x.month
    cls = close_price.groupby(GroupBy).nth(-1)
    opn = open_price.groupby(GroupBy).nth(0)
    rev = (cls - opn) / opn*100
    rev.index = pd.to_datetime(rev.index.astype(str), format='%Y%m') + MonthEnd()
    if save_data:
        rev.to_pickle(data_path+'Reversal_M')
    else:
        return rev