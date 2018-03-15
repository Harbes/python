# used to accomplish the IVOL in 《empirical asset pricing》
# 1. calculate the market beta
# 2. calculate the returns of two mimicking portfolios--- size & BM

import pandas as pd
import numpy as np
from time import time
from pandas.tseries.offsets import MonthEnd
#import matplotlib.pyplot as plt
data_path = 'E:/data/NewData/'  # '/Users/harbes/data/NewData/'#
def import_pv_index():
    global pv,open_price,close_price,index_ret,rtn,stock_pool
    pv = pd.read_pickle(data_path + 'PV_datetime')#[['adj_close', 'adj_open', 'size_tot']]
    close_price = pv['adj_close'].unstack()
    filter_ = pd.read_pickle(data_path + 'filtered_data')
    close_price=close_price[filter_]
    # close_price.index=pd.to_datetime(close_price.index.astype(int).astype(str),format='%Y%m%d')
    open_price = pv['adj_open'].unstack()
    # open_price.index=pd.to_datetime(open_price.index.astype(int).astype(str),format='%Y%m%d')
    rtn = (close_price - open_price) / open_price * 100
    rtn[rtn==0]=np.nan
    stock_pool=rtn.columns.values

    index_ret = pd.read_pickle(data_path + 'index_ret').set_index(['index_code', 'trddt'])['pctchange'].loc['000016.SH']
    index_ret.index = pd.to_datetime(index_ret.index.astype(int).astype(str), format='%Y%m%d')
    index_ret = pd.DataFrame({'pctchange': index_ret,
                              'month': (index_ret.index.year - index_ret.index[0].year) * 12 + index_ret.index.month})
    index_ret.drop(index_ret.loc['20160907'].name,inplace=True)
    # 注意，rtn的最后两列是index return 和 月份序号
    rtn['index'] = index_ret['pctchange']#目的是通过cov计算beta，同时通过market_beta是否为1检验code是否正确
    rtn['month'] = (rtn.index.year - rtn.index[0].year) * 12 + rtn.index.month

def cal_book(freq='M'):
    '''
    注意：book是由assets-liability得来的，数据中有负数存在
    monthly
    :return:
    '''
    BS = pd.read_pickle(data_path + '/BS')[['fin_year','ann_dt', 'stkcd', 'tot_assets', 'tot_liab']]#
    BS=BS[~(BS['ann_dt'].isnull())].set_index(['ann_dt', 'stkcd']).sort_index()
    book = BS['tot_assets'] - BS['tot_liab']
    book = book.drop(book.index[book.index.duplicated(keep='last')]).unstack() * 1e-8
    book.index = pd.to_datetime(book.index.astype(int).astype(str), format='%Y%m%d')
    book.index.name = 'trddt'
    book.loc[pd.to_datetime('20180228',format='%Y%m%d')]=np.nan # 为了resample时，将最后时间定在2018-2-28
    if freq=='M':
        book = book.resample('M').first().ffill()
        # 下面这一小段代码似乎不需要，只是为了防止日期不是月末
        By = lambda x: x.year * 100 + x.month
        book = book.groupby(By).nth(0)
        book.index = pd.to_datetime(book.index.astype(str), format='%Y%m') + MonthEnd()
    else:
        book = book.resample('D').first().ffill()
        book = book.reindex(rtn.index)
    return book[book > 0]

def check_match_index_columns():
    # TODO
    pass
def cal_beta(periods,save_data=False):
    '''
    calculate the market beta
    :param periods: int{measure periods}
    :param save_data:
    :return:
    '''
    beta = pd.DataFrame(index=pd.date_range('20050101', '20180301', freq='M'), columns=rtn.columns[:-1])
    for i in range(periods, 159):
        # 也可以通过cov/var，但是速度较慢
        tmp1=index_ret['pctchange'][(index_ret['month'] >= i - periods+1) & (index_ret['month'] <= i)]
        tmp2=rtn.iloc[:, :-1][(i - periods+1 <= rtn['month']) & (rtn['month'] <= i)]
        beta.iloc[i - 1] = ((tmp1-tmp1.mean()).values[:,None]*(tmp2-tmp2.mean())).sum()/((tmp1-tmp1.mean())**2).sum()
    if save_data:
        beta[beta!=0].to_pickle(data_path + 'beta_daily_'+str(periods)+'M')
    else:
        return beta[beta!=0]
def cal_size(freq='M',save_data=False):
    '''
    :return:
    '''
    if freq=='M':
        GroupBy = lambda x: x.year * 100 + x.month
        size = pv['size_tot'].unstack().groupby(GroupBy).last() * 1e-4
        size.index = pd.to_datetime(size.index.astype(str), format='%Y%m') + MonthEnd()
    else:
        size = pv['size_tot'].unstack() * 1e-4
    if save_data:
        size[size > 0].to_pickle(data_path+'size_'+freq)
    return size[size>0]
def cal_BM(size,freq='M',save_data=False):
    book=cal_book(freq=freq)
    #size=cal_size(freq=freq)
    BM=book/size
    if save_data:
        BM[BM>0].to_pickle(data_path+'BM_'+freq)
    else:
        return BM[BM>0]
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
        rev[rev!=0].to_pickle(data_path+'Reversal_M')
    else:
        return rev[rev!=0]
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
    coskew = pd.DataFrame(index=pd.date_range('20050101', '20180301', freq='M'), columns=rtn.columns[:-1])
    index_ret.insert(2, 'index^2', index_ret['pctchange'] ** 2)
    for i in range(periods,159):
        tmp1=rtn[(rtn['month'] >= i - periods + 1) & (rtn['month'] <= i)].iloc[:, :-1]
        tmp2=index_ret[(rtn['month'] >= i - periods + 1) & (rtn['month'] <= i)][['pctchange','index^2']]
        tmp3=np.linalg.pinv((tmp2-tmp2.mean()).values.T@(tmp2-tmp2.mean()))[1]
        coskew.iloc[i-1]=tmp3[0]*((tmp2 - tmp2.mean()).values[:,0][:,None] * (tmp1 - tmp1.mean())).sum()\
                         +tmp3[1]*((tmp2 - tmp2.mean()).values[:,1][:,None] * (tmp1 - tmp1.mean())).sum()
    index_ret.drop('index^2', axis=1, inplace=True)
    if save_data:
        coskew[coskew!=0].to_pickle(data_path+'coskew_'+str(periods)+'M')
    else:
        return coskew[coskew!=0]
def cal_iskew():
    # TODO
    #cal_SMB
    #cal_HML
    pass
def cal_vol(periods,save_data=False):
    t_vol=pd.DataFrame(index=pd.date_range('20050101','20180301',freq='M'),columns=rtn.columns[:-2])
    for i in range(periods,159):
        t_vol.iloc[i-1]=rtn[(rtn['month'] >= i - periods + 1) & (rtn['month'] <= i)].iloc[:, :-2].std()
    t_vol *= np.sqrt(12/periods)
    if save_data:
        t_vol.astype(float).to_pickle(data_path+'TotVol_'+str(periods)+'M')
    else:
        return t_vol.astype(float)
def cal_vol_ss(periods,save_data=False):
    vol_ss=pd.DataFrame(index=pd.date_range('20050101','20180301',freq='M'),columns=rtn.columns[:-2])
    for i in range(periods,159):
        vol_ss.iloc[i-1]=(rtn[(rtn['month'] >= i - periods + 1) & (rtn['month'] <= i)].iloc[:, :-2]**2).mean()
    vol_ss = np.sqrt(12 / periods)*np.sqrt(vol_ss.astype(float))
    if save_data:
        vol_ss.to_pickle(data_path + 'Vol_SS' + str(periods) + 'M')
    else:
        return vol_ss
def cal_ivol(periods,method='CAPM'):
    ivol = pd.DataFrame(index=pd.date_range('20050101', '20180301', freq='M'), columns=rtn.columns[:-2])
    if method =='CAPM':
        for i in range(periods, 159):
            # 也可以通过cov/var，但是速度较慢
            tmp1 = index_ret['pctchange'][(index_ret['month'] >= i - periods + 1) & (index_ret['month'] <= i)]
            tmp2 = rtn.iloc[:, :-2][(i - periods + 1 <= rtn['month']) & (rtn['month'] <= i)]
            beta_tmp = ((tmp1 - tmp1.mean()).values[:, None] * (tmp2 - tmp2.mean())).sum() / (
                        (tmp1 - tmp1.mean()) ** 2).sum()
            ivol.iloc[i-1]=(tmp2-tmp1[:,None]@beta_tmp.values[None,:]).std(axis=0)
    elif method=='FF':
        fSMB=SMB.iloc[:,0]-SMB.iloc[:,-1]
        fHML=HML.iloc[:,0]-HML.iloc[:,-1]
        X=pd.DataFrame({'pctchange':index_ret['pctchange'],'SMB':fSMB,'HML':fHML,'month':index_ret['month']})
        for i in range(periods, 159):
            # 也可以通过cov/var，但是速度较慢
            tmp1 = X[(X['month'] >= i - periods + 1) & (X['month'] <= i)][['pctchange', 'SMB', 'HML']]
            tmp2 = rtn.iloc[:, :-2][(i - periods + 1 <= rtn['month']) & (rtn['month'] <= i)]
            XY=pd.DataFrame({'MKT':(tmp1['pctchange'].values[:,None]*tmp2).mean(),
                             'SMB':(tmp1['SMB'].values[:,None]*tmp2).mean(),
                             'HML':(tmp1['HML'].values[:, None] * tmp2).mean()}).T
            beta_tmp=np.linalg.pinv((tmp1-tmp1.mean()).cov())@(XY.loc[['MKT','SMB','HML']])
            ivol.iloc[i-1]=(tmp2-tmp1@beta_tmp).std(axis=0)
    ivol=ivol.astype(float)*np.sqrt(12/periods)
    return ivol

def cal_SMB_HML(freq='M',weight=False):
    size = cal_size(freq=freq)
    BM = cal_BM(size,freq=freq)
    stock_pool=BM.columns&size.columns
    size=size[stock_pool]
    BM=BM[stock_pool]
    if freq=='M':
        ret = cal_rev()
        ret = ret[ret!=0][stock_pool]
    else:
        ret = rtn[rtn!=0].iloc[:,:-2][stock_pool]
    if weight:
        return cal_mimick_port1(size['2005':'2018-02'], ret['2005':'2018-02'], size),cal_mimick_port1(BM['2005':'2018-02'],ret['2005':'2018-02'],size)
    else:
        return cal_mimick_port1(size['2005':'2018-02'], ret['2005':'2018-02'], None),cal_mimick_port1(BM['2005':'2018-02'],ret['2005':'2018-02'],None)
def cal_mimick_port1(indi,rtn,weights):
    '''
    用法：cal_mimick_port(BM['2005':'2017'],rtn['2005':'2017'],None)或者
         cal_mimick_port(BM['2005':'2017'],rtn['2005':'2017'],size['2005':'2017'])
    注意：输入indi时，一定要保证是从非NA数据开始的;而且数据一定要对齐！！！
    :param indi:
    :return:
    '''
    group_num=5
    percentile = np.linspace(0, 1, group_num + 1) #[0.0,0.3,0.7,1.0]# 也可以自定义percentile，例如
    label_ = [i + 1 for i in range(len(percentile) - 1)]
    mark_ = pd.DataFrame([pd.qcut(indi.iloc[i],q=percentile, labels=label_) for i in range(len(indi)-1)],
                      index=indi.index[1:]) # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
    valid_=~(pd.isnull(mark_) | pd.isnull(rtn.iloc[1:])) # valid的股票要满足：当期有前一个月的indicator信息；当期保证交易
    if weights is None:
        df = pd.DataFrame()
        df['rtn'] = rtn[valid_].stack()
        df['ref'] = mark_[valid_].stack()
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
def cal_mimick_port2(indi1,indi2,rtn,weights,independent=True):

    group_num = 5
    percentile = np.linspace(0, 1, group_num + 1)  # 也可以自定义percentile，例如 [0.0,0.3,0.7,1.0]#
    label_ = [i + 1 for i in range(len(percentile) - 1)]
    if independent:
        mark_1 = pd.DataFrame([pd.qcut(indi1.iloc[i], q=percentile, labels=label_) for i in range(len(indi1) - 1)],
                              index=indi1.index[1:])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
        mark_2 = pd.DataFrame([pd.qcut(indi2.iloc[i], q=percentile, labels=label_) for i in range(len(indi2) - 1)],
                              index=indi2.index[1:])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
    else:
        mark_1 = pd.DataFrame([pd.qcut(indi1.iloc[i], q=percentile, labels=label_) for i in range(len(indi1) - 1)],
                              index=indi1.index[1:])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
        mark_2=pd.DataFrame(index=mark_1.index,columns=mark_1.columns)
        for l_ in label_:
            tmp=pd.DataFrame([pd.qcut(indi2.iloc[i][mark_1.iloc[i]==l_],q=percentile,labels=label_) for i in range(len(indi2)-1)])
            mark_2 = mark_2.combine_first(tmp)
    valid_ = ~(pd.isnull(mark_1+mark_2) | pd.isnull(rtn[1:]))  # valid的股票要满足：当期有前一个月的indicator信息；当期保证交易
    if weights is None:
        df = pd.DataFrame()
        df['rtn'] = rtn[valid_].stack()
        df['ref1'] = mark_1.stack()
        df['ref2'] = mark_2.stack()
        tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).mean()).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    else:
        df1 = pd.DataFrame()
        df2 = pd.DataFrame()
        df1['rtn_w'] = (rtn * weights)[valid_].stack()
        df1['ref1'] = mark_1.stack()
        df1['ref2'] = mark_2.stack()
        df2['rtn_w'] = weights[valid_].stack()
        df2['ref1'] = mark_1.stack()
        df2['ref2'] = mark_2.stack()
        tmp1 = df1.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).sum())
        tmp2 = df2.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).sum())
        tmp = (tmp1 / tmp2).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    return tmp

def describe(df,stats):
    d=df.describe(percentiles=[0.05,0.25,0.5,0.75,0.95])
    return d.append(df.reindex(d.columns, axis=1).agg(stats))
def func_percentile(n):
    # TODO 速度比较慢，慎用
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = '%s' % n
    return percentile_
def summary_vol_ivol():
    measure_periods=(1,3,6,12)
    vol=pd.DataFrame(index=measure_periods,columns=['count','mean','std','min','5%','25%','50%','75%','95%','max','skew','kurt'])
    vol_ss=pd.DataFrame(index=measure_periods,columns=['count','mean','std','min','5%','25%','50%','75%','95%','max','skew','kurt'])
    ivol_CAPM=pd.DataFrame(index=measure_periods,columns=['count','mean','std','min','5%','25%','50%','75%','95%','max','skew','kurt'])
    ivol_FF=pd.DataFrame(index=measure_periods,columns=['count','mean','std','min','5%','25%','50%','75%','95%','max','skew','kurt'])
    for i in measure_periods:
        vol.loc[i]=describe(cal_vol(i).T,['skew','kurt']).mean(axis=1)
        vol_ss.loc[i] = describe(cal_vol_ss(i).T, ['skew', 'kurt']).mean(axis=1)
        ivol_CAPM.loc[i]=describe(cal_ivol(i).T, ['skew', 'kurt']).mean(axis=1)
        ivol_FF.loc[i] = describe(cal_ivol(i,method='FF').T, ['skew', 'kurt']).mean(axis=1)
    return vol,vol_ss,ivol_CAPM,ivol_FF
if __name__ == '__main__':
    import_pv_index()
    #import_book(freq='D')
    vol, vol_ss, ivol_CAPM, ivol_FF=summary_vol_ivol() # 不同测度期限相差较大（这点与美国不同）
    vol
    vol_ss
    ivol_CAPM
    ivol_FF
