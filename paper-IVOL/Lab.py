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
    pv = pd.read_pickle(data_path + 'PV_datetime')#[['adj_close', 'adj_open', 'size_tot']]
    close_price = pv['adj_close'].unstack()
    # close_price.index=pd.to_datetime(close_price.index.astype(int).astype(str),format='%Y%m%d')
    open_price = pv['adj_open'].unstack()
    # open_price.index=pd.to_datetime(open_price.index.astype(int).astype(str),format='%Y%m%d')
    index_ret = pd.read_pickle(data_path + 'index_ret').set_index(['index_code', 'trddt'])['pctchange'].loc['000016.SH']
    index_ret.index = pd.to_datetime(index_ret.index.astype(int).astype(str), format='%Y%m%d')
    rtn = (close_price - open_price) / open_price * 100
    # 注意，rtn的最后两列是index return 和 月份序号
    rtn['index'] = index_ret  # ['pctchange']

    rtn['month'] = (rtn.index.year - rtn.index[0].year) * 12 + rtn.index.month
    index_ret = pd.DataFrame({'pctchange': index_ret,
                              'month': (index_ret.index.year - index_ret.index[0].year) * 12 + index_ret.index.month})
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
    # 下面这一小段代码似乎不需要，只是为了防止日期不是月末
    By=lambda x:x.year*100+x.month
    book=book.groupby(By).nth(0)
    book.index= pd.to_datetime(book.index.astype(str), format='%Y%m')+MonthEnd()
    return book[book>0]


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
        # 也可以通过去均值的矩阵相乘
        beta.iloc[i - 1] = rtn.iloc[:, :-1][(i - periods+1 <= rtn['month']) & (rtn['month'] <= i)].cov(min_periods=20).iloc[-1,
                           :-1] / \
                           np.var(index_ret['pctchange'][(index_ret['month'] >= i - periods+1) & (index_ret['month'] <= i)])
    print(time() - t0)
    if save_data:
        beta.to_pickle(data_path + 'beta_daily_'+str(periods)+'M')
    else:
        return beta
def cal_size(save_data=False):
    '''
    monthly
    :return:
    '''
    GroupBy = lambda x: x.year * 100 + x.month
    size = pv['size_tot'].unstack().groupby(GroupBy).last()*1e-4
    size.index = pd.to_datetime(size.index.astype(str), format='%Y%m') + MonthEnd()
    if save_data:
        size.to_pickle(data_path+'size_monthly')
    return size
def cal_BM(save_data=False):
    book=import_book()
    size=cal_size()
    BM=book/size
    if save_data:
        BM.to_pickle(data_path+'BM_monthly')
    else:
        return BM
def cal_mom(periods,save_data=False):
    GroupBy = lambda x: x.year * 100 + x.month
    cls=close_price.groupby(GroupBy).nth(-1)
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
def cal_illiq(periods,save_data=False):
    GroupBy = lambda x: x.year * 100 + x.month
    amount=pv['amount'].unstack()
    rtn_amount=1e5*np.abs(close_price-open_price)/open_price/amount
    rtn_amount_sum=rtn_amount.groupby(GroupBy).sum()
    rtn_amount_count=rtn_amount.groupby(GroupBy).count()
    _sum=0;_count=0;
    for i in range(periods):
        _sum +=rtn_amount_sum.shift(i)
        _count +=rtn_amount_count.shift(i)
    illiq=_sum/_count
    illiq.index=pd.to_datetime(illiq.index.astype(str),format='%Y%m')+MonthEnd()
    illiq.iloc[:periods-1]=np.nan
    if save_data:
        illiq.to_pickle(data_path+'illiq_'+str(periods)+'M')
    else:
        return illiq
def cal_skew(periods,save_data=False):
    t_skew=pd.DataFrame(index=pd.date_range('20050101','20180301',freq='M'),columns=rtn.columns[:-2])
    for i in range(periods,159):
        t_skew.iloc[i-1]=rtn[(rtn['month']>=i-periods+1)&(rtn['month']<=i)].iloc[:,:-2].skew()
    if save_data:
        t_skew.to_pickle(data_path+'skew_'+str(periods)+'M')
    else:
        return t_skew
def cal_coskew(periods,save_data=False):
    coskew=pd.DataFrame(index=pd.date_range('20050101','20180301',freq='M'),columns=rtn.columns[:-2])
    rtn.insert(len(rtn.columns)-1,'index^2',rtn['index']**2);rtn.iloc[:5,-5:]
    for i in range(periods,159):
        tmp=rtn[(rtn['month'] >= i - periods + 1) & (rtn['month'] <= i)].iloc[:, :-1].cov()
        coskew.iloc[i-1]=(np.linalg.pinv(tmp.iloc[-2:,-2:])@tmp.iloc[-2:,:-2])[1]
    rtn.drop('index^2', axis=1, inplace=True)
    if save_data:
        coskew.to_pickle(data_path+'coskew_'+str(periods)+'M')
    else:
        return coskew
def cal_iskew():
    #cal_SMB
    #cal_HML
    pass
def cal_vol(periods,save_data=False):
    t_vol=pd.DataFrame(index=pd.date_range('20050101','20180301',freq='M'),columns=rtn.columns[:-2])
    for i in range(periods,159):
        t_vol.iloc[i-1]=rtn[(rtn['month'] >= i - periods + 1) & (rtn['month'] <= i)].iloc[:, :-2].std()
    t_vol *= np.sqrt(12/periods)
    if save_data:
        t_vol.to_pickle(data_path+'TotVol_'+str(periods)+'M')
    else:
        return t_vol
def cal_vol_ss():
    vol_ss=pd.DataFrame(index=pd.date_range('20050101','20180301',freq='M'),columns=rtn.columns[:-2])
    rtn2=rtn**2
    for i in range(periods,159):
        vol_ss.iloc[i-1]=rtn2[(rtn2['month'] >= i - periods + 1) & (rtn2['month'] <= i)].iloc[:, :-2].mean()
    vol_ss *= np.sqrt(12 / periods)
    if save_data:
        vol_ss.to_pickle(data_path + 'Vol_SS' + str(periods) + 'M')
    else:
        return vol_ss
def cal_ivol():
    pass
def GroupBySingleFactor():
    pass
def cal_SMB():
    group_num=5
    percentile = np.linspace(0, 1, group_num + 1)
    label_ = [i + 1 for i in range(len(percentile)-1)]
    mark_ = DataFrame([pd.qcut(past_rtn.iloc[i],
                                       q=percentile, labels=label_momentum) for i in range(J + M, len(price0))],
                              index=price0.index[J + M:])
def cal_HML():
    pass

def cal_mimick_port(indi,rtn,weights):
    '''
    注意：输入indi时，一定要保证是从非NA数据开始的
    :param indi:
    :return:
    '''
    percentile = np.linspace(0, 1, group_num + 1)
    label_ = [i + 1 for i in range(len(percentile) - 1)]
    mark_ = pd.DataFrame([pd.qcut(indi.iloc[i],q=percentile, labels=label_) for i in range(len(indi)-1)],
                      index=indi.index[1:]) # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
    valid_=~(np.isnan(indi.shift(1)) | np.isnan(rtn))
    if weights is None:
        df = pd.DataFrame()
        df['rtn'] = rtn[valid_].stack()
        df['ref'] = mark_.stack()
        tmp=df.groupby(level=0).apply(lambda g: g.groupby('ref').mean()).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    else:
        df1=pd.DataFrame()
        df2=pd.DataFrame()
        df1['rtn_w'] = (rtn*weights)[valid_].stack()
        df1['ref'] = mark_.stack()
        df2['rtn_w'] = weights[valid_].stack()
        df2['ref'] = mark_.stack()
        tmp1=df1.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp2 = df2.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp=(tmp1/tmp2).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    return tmp
BM=cal_BM()
size=cal_size()
BM=BM[BM.columns&size.columns]
size=size[BM.columns&size.columns]
a=cal_mimick_port(BM['2005':'2017'],rtn['2005':'2017'],None)#,size['2005':'2017'])
tmp=a[1]-a[5]
tmp.mean()/tmp.std()*np.sqrt(len(tmp))













