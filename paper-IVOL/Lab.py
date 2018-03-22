# used to accomplish the IVOL in 《empirical asset pricing》
# 1. calculate the market beta
# 2. calculate the returns of two mimicking portfolios--- size & BM

import pandas as pd
import numpy as np
from time import time
from pandas.tseries.offsets import MonthEnd
from dateutil.parser import parse
#import matplotlib.pyplot as plt
data_path = '/Users/harbes/data/NewData/'#'E:/data/NewData/'  #
def import_pv_index(del_limit_move=True):
    global pv,open_price,close_price,index_ret,rtn,stock_pool,rf
    rf = pd.read_excel(data_path + 'risk_free.xlsx')[2:].set_index(['Clsdt'])
    rf.index = pd.to_datetime(rf.index, format='%Y-%m-%d')
    pv = pd.read_pickle(data_path + 'PV_datetime')[['adj_close', 'adj_open', 'size_tot','amount']]
    close_price = pv['adj_close'].unstack()
    if del_limit_move:
        filter_ = pd.read_pickle(data_path + 'filtered_data')
    else:
        filter_=pd.read_pickle(data_path+'filtered_data_with_limit_move')
    close_price=close_price[filter_]
    # close_price.index=pd.to_datetime(close_price.index.astype(int).astype(str),format='%Y%m%d')
    open_price = pv['adj_open'].unstack()
    # open_price.index=pd.to_datetime(open_price.index.astype(int).astype(str),format='%Y%m%d')
    rtn = (close_price - open_price) / open_price * 100
    rtn[rtn==0]=np.nan
    rtn=rtn.sub(rf['Nrrdaydt'][rtn.index].astype(float), axis=0)
    stock_pool=rtn.columns.values

    index_ret = pd.read_pickle(data_path + 'index_ret').set_index(['index_code', 'trddt'])['pctchange'].loc['000016.SH']
    index_ret.index = pd.to_datetime(index_ret.index.astype(int).astype(str), format='%Y%m%d')
    index_ret = index_ret.sub(rf['Nrrdaydt'][index_ret.index].astype(float), axis=0)
    index_ret = pd.DataFrame({'pctchange': index_ret,
                              'month': (index_ret.index.year - index_ret.index[0].year) * 12 + index_ret.index.month})
    index_ret.drop(index_ret.loc['20160907'].name,inplace=True)
    # 注意，rtn的最后两列是index return 和 月份序号
    rtn['index'] = index_ret['pctchange']#目的是通过cov计算beta，同时通过market_beta是否为1检验code是否正确
    rtn['month'] = (rtn.index.year - rtn.index[0].year) * 12 + rtn.index.month


def cal_index_ret(freq='M'):
    index_ret = pd.read_pickle(data_path + 'index_ret').set_index(['index_code', 'trddt'])[['opnprc','clsprc']].loc['000016.SH']
    By=lambda x:x//100.0
    index_cls=index_ret['clsprc'].groupby(By).last()
    index_opn=index_ret['opnprc'].groupby(By).first()
    index_ret = (index_cls - index_opn) / index_opn * 100.0
    index_ret.index = pd.to_datetime(index_ret.index.astype(int).astype(str), format='%Y%m')+MonthEnd()
    index_ret = index_ret.sub(rf['Nrrmtdt'][index_ret.index].astype(float), axis=0)
    return index_ret

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
        beta[beta!=0].astype(float).to_pickle(data_path + 'beta_daily_'+str(periods)+'M')
    else:
        return beta[beta!=0].astype(float).iloc[:,:-1] # 剔除 column：index
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
    mom=(cls.shift(1)-opn.shift(periods-1))/opn.shift(periods-1)*100
    mom.index = pd.to_datetime(mom.index.astype(str), format='%Y%m')+MonthEnd()
    if save_data:
        mom.to_pickle(data_path+'Mom_'+str(periods)+'M')
    else:
        return mom
def cal_rev(del_rf=False,first_day=True,save_data=False):
    GroupBy = lambda x: x.year * 100 + x.month
    if first_day:
        cls = close_price.groupby(GroupBy).nth(-1)
        opn = open_price.groupby(GroupBy).nth(0)
    else:
        cls = close_price.groupby(GroupBy).last()
        opn = open_price.groupby(GroupBy).first()
    rev = (cls - opn) / opn*100
    rev.index = pd.to_datetime(rev.index.astype(str), format='%Y%m') + MonthEnd()
    if del_rf:
        rev=rev.sub(rf['Nrrmtdt'][rev.index].astype(float), axis=0)
    if save_data:
        rev[rev!=0].to_pickle(data_path+'Reversal_M')
    else:
        return rev#[rev!=0]
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
        t_skew.astype(float).to_pickle(data_path+'skew_'+str(periods)+'M')
    else:
        return t_skew.astype(float)
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
        coskew[coskew!=0].astype(float).to_pickle(data_path+'coskew_'+str(periods)+'M')
    else:
        return coskew[coskew!=0].astype(float).iloc[:,:-1]
def cal_iskew(periods,method='CAPM'):
    iskew = pd.DataFrame(index=pd.date_range('20050101', '20180301', freq='M'), columns=rtn.columns[:-2])
    if method == 'CAPM':
        for i in range(periods, 159):
            # 也可以通过cov/var，但是速度较慢
            tmp1 = index_ret['pctchange'][(index_ret['month'] >= i - periods + 1) & (index_ret['month'] <= i)]
            tmp2 = rtn.iloc[:, :-2][(i - periods + 1 <= rtn['month']) & (rtn['month'] <= i)]
            beta_tmp = ((tmp1 - tmp1.mean()).values[:, None] * (tmp2 - tmp2.mean())).sum() / (
                    (tmp1 - tmp1.mean()) ** 2).sum()
            iskew.iloc[i - 1] = (tmp2 - tmp1[:, None] @ beta_tmp.values[None, :]).skew()
    else:
        '''
        FF方法依赖于存在SMB和HML多空组合的return
        '''
        fSMB = SMB.iloc[:, 0] - SMB.iloc[:, -1]
        fHML = HML.iloc[:, 0] - HML.iloc[:, -1]
        X = pd.DataFrame({'pctchange': index_ret['pctchange'], 'SMB': fSMB, 'HML': fHML, 'month': index_ret['month']})
        for i in range(periods, 159):
            # 也可以通过cov/var，但是速度较慢
            tmp1 = X[(X['month'] >= i - periods + 1) & (X['month'] <= i)][['pctchange', 'SMB', 'HML']]
            tmp2 = rtn.iloc[:, :-2][(i - periods + 1 <= rtn['month']) & (rtn['month'] <= i)]
            XY = pd.DataFrame(
                {'MKT': ((tmp1['pctchange'] - tmp1['pctchange'].mean()).values[:, None] * (tmp2 - tmp2.mean())).mean(),
                 'SMB': ((tmp1['SMB'] - tmp1['SMB'].mean()).values[:, None] * (tmp2 - tmp2.mean())).mean(),
                 'HML': ((tmp1['HML'] - tmp1['HML'].mean()).values[:, None] * (tmp2 - tmp2.mean())).mean()}).T
            beta_tmp = np.linalg.pinv((tmp1 - tmp1.mean()).cov()) @ (XY.loc[['MKT', 'SMB', 'HML']])
            iskew.iloc[i - 1] = (tmp2 - tmp1 @ beta_tmp).skew()
    return iskew.astype(float)
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
def cal_ivol(periods,SMB,HML,method='CAPM'):
    ivol = pd.DataFrame(index=pd.date_range('20050101', '20180301', freq='M'), columns=rtn.columns[:-2])
    if method =='CAPM':
        for i in range(periods, 159):
            # 也可以通过cov/var，但是速度较慢
            tmp1 = index_ret['pctchange'][(index_ret['month'] >= i - periods + 1) & (index_ret['month'] <= i)]
            tmp2 = rtn.iloc[:, :-2][(i - periods + 1 <= rtn['month']) & (rtn['month'] <= i)]
            beta_tmp = ((tmp1 - tmp1.mean()).values[:, None] * (tmp2 - tmp2.mean())).sum() / (
                        (tmp1 - tmp1.mean()) ** 2).sum()
            ivol.iloc[i-1]=(tmp2-tmp1[:,None]@beta_tmp.values[None,:]).std()
    else:
        '''
        FF方法依赖于存在SMB和HML多空组合的return
        '''
        fSMB=SMB.iloc[:,0]-SMB.iloc[:,-1]
        fHML=HML.iloc[:,0]-HML.iloc[:,-1]
        X=pd.DataFrame({'pctchange':index_ret['pctchange'],'SMB':fSMB,'HML':fHML,'month':index_ret['month']})
        for i in range(periods, 159):
            # 也可以通过cov/var，但是速度较慢
            tmp1 = X[(X['month'] >= i - periods + 1) & (X['month'] <= i)][['pctchange', 'SMB', 'HML']]
            tmp2 = rtn.iloc[:, :-2][(i - periods + 1 <= rtn['month']) & (rtn['month'] <= i)]
            XY=pd.DataFrame({'MKT':((tmp1['pctchange']-tmp1['pctchange'].mean()).values[:,None]*(tmp2-tmp2.mean())).mean(),
                             'SMB':((tmp1['SMB']-tmp1['SMB'].mean()).values[:,None]*(tmp2-tmp2.mean())).mean(),
                             'HML':((tmp1['HML']-tmp1['HML'].mean()).values[:, None] * (tmp2-tmp2.mean())).mean()}).T
            beta_tmp=np.linalg.pinv((tmp1-tmp1.mean()).cov())@(XY.loc[['MKT','SMB','HML']])
            ivol.iloc[i-1]=(tmp2-tmp1@beta_tmp).std()
    ivol=ivol.astype(float)*np.sqrt(12/periods)
    return ivol
def cal_vol_year(periods,ret):
    mini_periods={12:10,24:20,36:24}
    return ret.rolling(periods,min_periods=mini_periods[periods]).std()*np.sqrt(12/periods)
def cal_vol_ss_year(periods,ret):
    mini_periods={12:10,24:20,36:24}
    return np.sqrt((ret**2).rolling(periods,min_periods=mini_periods[periods]).mean()*12/periods)
def cal_ivol_year(periods,ret,index_ret,SMB,HML,method='CAPM'):
    ivol = pd.DataFrame(index=pd.date_range('20050101', '20180301', freq='M'), columns=ret.columns)
    if method == 'CAPM':
        for i in range(periods, 159):
            # 也可以通过cov/var，但是速度较慢
            tmp1 = index_ret.iloc[i-periods:i]
            tmp2 = ret.iloc[i-periods:i]
            beta_tmp = ((tmp1 - tmp1.mean()).values[:, None] * (tmp2 - tmp2.mean())).sum() / (
                    (tmp1 - tmp1.mean()) ** 2).sum()
            ivol.iloc[i - 1] = (tmp2 - tmp1[:, None] @ beta_tmp.values[None, :]).std()
    else:
        '''
        FF方法依赖于存在SMB和HML多空组合的return
        '''
        #SMB,HML=cal_SMB_HML()
        fSMB = SMB.iloc[:, 0] - SMB.iloc[:, -1]
        fHML = HML.iloc[:, 0] - HML.iloc[:, -1]
        X = pd.DataFrame({'pctchange': index_ret, 'SMB': fSMB, 'HML': fHML})
        for i in range(periods, 159):
            # 也可以通过cov/var，但是速度较慢
            tmp1 = X.iloc[i-periods:i]
            tmp2 = ret.iloc[i-periods:i]
            XY = pd.DataFrame(
                {'MKT': ((tmp1['pctchange'] - tmp1['pctchange'].mean()).values[:, None] * (tmp2 - tmp2.mean())).mean(),
                 'SMB': ((tmp1['SMB'] - tmp1['SMB'].mean()).values[:, None] * (tmp2 - tmp2.mean())).mean(),
                 'HML': ((tmp1['HML'] - tmp1['HML'].mean()).values[:, None] * (tmp2 - tmp2.mean())).mean()}).T
            beta_tmp = np.linalg.pinv((tmp1 - tmp1.mean()).cov()) @ (XY.loc[['MKT', 'SMB', 'HML']])
            ivol.iloc[i - 1] = (tmp2 - tmp1 @ beta_tmp).std()
    ivol = ivol.astype(float) * np.sqrt(12 / periods)
    return ivol
def cal_SMB_HML(freq='M',group_num=10,weight=False):
    size = cal_size(freq=freq)
    BM = cal_BM(size,freq=freq)
    stock_pool=BM.columns&size.columns
    size=size[stock_pool]
    BM=BM[stock_pool]
    if freq=='M':
        ret = cal_rev()
        ret = ret[ret!=0][stock_pool].sub(rf['Nrrmtdt'][ret.index], axis=0).astype(float)
    else:
        ret = rtn.iloc[:,:-2][stock_pool]
    if weight:
        return cal_mimick_port1(size['2005':'2018-02'], ret['2005':'2018-02'], size,group_num),cal_mimick_port1(BM['2005':'2018-02'],ret['2005':'2018-02'],size,group_num)
    else:
        return cal_mimick_port1(size['2005':'2018-02'], ret['2005':'2018-02'], None,group_num),cal_mimick_port1(BM['2005':'2018-02'],ret['2005':'2018-02'],None,group_num)
def cal_mimick_port1(indi,rtn,weights,group_num):
    '''
    用法：cal_mimick_port(BM['2005':'2017'],rtn['2005':'2017'],None)或者
         cal_mimick_port(BM['2005':'2017'],rtn['2005':'2017'],size['2005':'2017'])
    注意：输入indi时，一定要保证是从非NA数据开始的;而且数据一定要对齐！！！
    :param indi:
    :return:
    '''
    #indi=vol.loc['200512':'201802']
    #rtn=ret.loc['200512':'201802']
    #group_num=5
    percentile = np.linspace(0, 1, group_num + 1) #[0.0,0.3,0.7,1.0]# 也可以自定义percentile，例如
    label_ = [i + 1 for i in range(len(percentile) - 1)]
    mark_ = pd.DataFrame([pd.qcut(indi.iloc[i],q=percentile, labels=label_) for i in range(len(indi)-1)],index=indi.index[1:]) # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
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
            tmp=pd.DataFrame([pd.qcut(indi2.iloc[i][mark_1.iloc[i]==l_],q=percentile,labels=label_) for i in range(len(indi2)-1)],index=mark_1.index)
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

def describe(df,stats=['skew','kurt']):
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
def summary_vol_ivol_year():
    index_ret=cal_index_ret()
    ret=cal_rev(del_rf=True)
    measure_periods=(12, 24, 36)
    vol = pd.DataFrame(index=measure_periods,
                       columns=['count', 'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max', 'skew', 'kurt'])
    vol_ss = pd.DataFrame(index=measure_periods,
                          columns=['count', 'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max', 'skew',
                                   'kurt'])
    ivol_CAPM = pd.DataFrame(index=measure_periods,
                             columns=['count', 'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max', 'skew',
                                      'kurt'])
    ivol_FF = pd.DataFrame(index=measure_periods,
                           columns=['count', 'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max', 'skew',
                                    'kurt'])
    for i in measure_periods:
        vol.loc[i] = describe(cal_vol_year(i,ret).iloc[i-1:].T, ['skew', 'kurt']).mean(axis=1)
        vol_ss.loc[i] = describe(cal_vol_ss_year(i,ret).iloc[i-1:].T, ['skew', 'kurt']).mean(axis=1)
        ivol_CAPM.loc[i] = describe(cal_ivol_year(i,ret,index_ret).T, ['skew', 'kurt']).mean(axis=1)
        ivol_FF.loc[i] = describe(cal_ivol_year(i,ret,index_ret,method='FF').T, ['skew', 'kurt']).mean(axis=1)
    return vol,vol_ss,ivol_CAPM,ivol_FF
def summary_skew_coskew_iskew():
    measure_periods = (1, 3, 6, 12)
    skew = pd.DataFrame(index=measure_periods,
                       columns=['count', 'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max', 'skew', 'kurt'])
    coskew = pd.DataFrame(index=measure_periods,
                          columns=['count', 'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max', 'skew',
                                   'kurt'])
    iskew_CAPM = pd.DataFrame(index=measure_periods,
                             columns=['count', 'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max', 'skew',
                                      'kurt'])
    iskew_FF = pd.DataFrame(index=measure_periods,
                           columns=['count', 'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max', 'skew',
                                    'kurt'])
    for i in measure_periods:
        skew.loc[i]=describe(cal_skew(i).T,['skew','kurt']).mean(axis=1)
        coskew.loc[i] = describe(cal_coskew(i).T, ['skew', 'kurt']).mean(axis=1)
        iskew_CAPM.loc[i]=describe(cal_iskew(i).T, ['skew', 'kurt']).mean(axis=1)
        iskew_FF.loc[i] = describe(cal_iskew(i,method='FF').T, ['skew', 'kurt']).mean(axis=1)
    return skew,coskew,iskew_CAPM,iskew_FF
def return_vol_ivol():
    df1=pd.DataFrame()
    df2=pd.DataFrame()
    df3=pd.DataFrame()
    df4=pd.DataFrame()
    SMB,HML=cal_SMB_HML(freq='D')
    for i in (1,3,6,12):
        df1['vol_'+str(i)+'M']=cal_vol(i).stack()
        df2['vol_ss_'+str(i)+'M']=cal_vol_ss(i).stack()
        df3['ivol_CAPM_'+str(i)+'M']=cal_ivol(i,SMB,HML).stack()
        df4['ivol_FF_'+str(i)+'M']=cal_ivol(i,SMB,HML,method='FF').stack()
    return df1, df2, df3, df4
def return_vol_ivol_year():
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    df3 = pd.DataFrame()
    df4 = pd.DataFrame()
    SMB, HML = cal_SMB_HML(freq='M')
    ret=cal_rev(del_rf=True)
    index_ret=cal_index_ret()
    for j in (1,2,3):
        df1['vol_'+str(j)+'Y']=cal_vol_year(j*12,ret).stack()
        df2['vol_ss_'+str(j)+'Y']=cal_vol_ss_year(j*12,ret).stack()
        df3['ivol_CAPM_'+str(j)+'Y']=cal_ivol_year(j*12,ret,index_ret,SMB,HML).stack()
        df4['ivol_FF_'+str(j)+'Y']=cal_ivol_year(j*12,ret,index_ret,SMB,HML,method='FF').stack()
    return   df1,df2,df3,df4

def NWest_mean(d,L=None):
    '''
    df不要以Nan开头，会引起误差;
    :param df:
    :param L:
    :return:
    '''
    df=d.copy()
    df-=df.mean()
    T=len(df)
    if L is None:
        L=int(T**0.25)
    w=1.0-np.arange(1,L+1)/(L+1.0)
    return np.sqrt(2.0*pd.DataFrame((df*df.shift(i+1)*w[i]).sum() for i in range(L)).sum()/T+df.var())/np.sqrt(T)
def NWest(e,X,L=None):
    T = len(e)
    if L is None:
        L = int(T ** 0.25) # or : L = 0.75*T**(1.0/3.0)-1
    w = 1.0 - np.arange(1, L + 1) / (L+1.0)
    X.insert(0,'c',np.ones(T))
    S=0.0
    for l in range(1,L+1):
        for i in range(l,T):
            S+=w[l-1]*e[i]*e[i-l]*(X.iloc[i][:,None]*X.iloc[i-l].values+X.iloc[i-l][:,None]*X.iloc[i].values)
    for i in range(T):
        S+=e[i]*e[i]*X.iloc[i][:,None]*X.iloc[i].values
    XX_1=np.linalg.pinv(X.T@X.values)
    X.drop('c', axis=1, inplace=True)
    return np.sqrt((XX_1@S@XX_1)[0,0])
def cal_FF_alpha(arr,period,index_ret,fSMB,fHML):
    '''

    :param arr:
    :param period:
    :param index_ret:
    :param fSMB:
    :param fHML:
    :return: alpha, s.e. of estimated alpha
    '''
    start={1:'200502',
           3:'200504',
           6:'200507',
           12:'200601',
           24:'200701',
           36:'200801'} # 回归有效数据的起始月
    X=pd.DataFrame({'index':index_ret,'SMB':fSMB,'HML':fHML})[['index','SMB','HML']].loc[start[period]:'201802']
    arr=arr.loc[start[period]:'201802']
    tmp1=X-X.mean(axis=0)
    tmp2=arr-arr.mean(axis=0)
    tmp3=arr-X@np.linalg.pinv(tmp1.cov())@(tmp1.mul(tmp2,axis=0).mean()).values
    return tmp3.mean(),NWest(tmp3-tmp3.mean(),X)


if __name__ == '__main__':
    import_pv_index()
    #import_book(freq='D')
    # summary statistics------vol
    SMB, HML = cal_SMB_HML(freq='M')
    vol_y, vol_ss_y, ivol_CAPM_y, ivol_FF_y = summary_vol_ivol_year()
    SMB,HML=cal_SMB_HML(freq='D')
    vol, vol_ss, ivol_CAPM, ivol_FF=summary_vol_ivol() # 不同测度期限相差较大（这点与美国不同）;FF3-ivol存在更大的偏度，这是为什么，反映了什么
    vol.to_csv(data_path+'Vol_.csv')
    vol_ss.to_csv(data_path+'Vol_SS_.csv')
    ivol_CAPM.to_csv(data_path+'iVol_CAPM_.csv')
    ivol_FF.to_csv(data_path+'iVol_FF_.csv')
    vol_y.to_csv(data_path+'Vol_1')
    vol_ss_y.to_csv(data_path+'Vol_SS_1.csv')
    ivol_CAPM_y.to_csv(data_path+'iVol_CAPM_1.csv')
    ivol_FF_y.to_csv(data_path+'iVol_FF_1.csv')
    # correlation------ within vol and between vol
    start_ = parse('20071201')
    end_ = parse('20180228')
    vol, vol_ss, ivol_CAPM, ivol_FF = return_vol_ivol()
    vol_y, vol_ss_y, ivol_CAPM_y, ivol_FF_y = return_vol_ivol_year()
    tmp=pd.concat([vol,vol_y,vol_ss,vol_ss_y,ivol_CAPM,ivol_CAPM_y,ivol_FF,ivol_FF_y],axis=1)
    corr_=tmp.loc[start_:end_].groupby(level=0).corr('pearson').groupby(level=1).mean() # pearson # spearman
    corr_=corr_.loc[corr_.columns];corr_
    #corr_.to_csv(data_path+'pearson_correlation.csv')
    # correlations------ivol and other variables
    beta=cal_beta(12).stack()
    size=cal_size(freq='M')#.stack()
    BM=cal_BM(size).stack();size=np.log(size).stack() #size=size.stack() #
    MOM=cal_mom(12).stack()
    rev=cal_rev().stack()
    illiq=cal_illiq(12).stack()
    skew=cal_skew(12).stack()
    coskew=cal_coskew(12).stack()
    iskew=cal_iskew(12).stack()
    beta.name='beta'
    size.name='size'
    BM.name='BM'
    MOM.name='mom'
    rev.name='rev'
    illiq.name='illiq'
    skew.name='skew'
    coskew.name='coskew'
    iskew.name='iskew'
    tmp=pd.concat([ivol_FF,ivol_FF_y,beta,size,BM,MOM,rev,illiq,skew,coskew,iskew],axis=1);tmp.head()
    #tmp['beta']=beta;tmp['size']=size
    #tmp['BM']=BM
    #tmp['mom']=MOM
    #tmp['rev']=rev
    #tmp['illiq']=illiq
    #tmp['coskew']=coskew.iloc[:,:-1] # 剔除市场
    #tmp['iskew']=iskew
    corr_ = tmp.loc[start_:end_].groupby(level=0).corr('spearman').groupby(level=1).mean()  # pearson # spearman
    corr_ = corr_.loc[corr_.columns];corr_
    corr_.to_csv(data_path+'spearman_corr_between_ivol_other_variables.csv')

    s=tmp.groupby(level=0).apply(describe)
    s[s==0]=np.nan
    s.mean(level=1).to_csv(data_path+'SummaryStatisticsOfOtherVariables.csv')

    del beta, size, BM, MOM, rev, illiq, skew, coskew, iskew

    # persistence analysis
    df= pd.concat([ivol_FF, ivol_FF_y], axis=1) ;columns_=df.columns
    tau = 12
    df_tau= df.groupby(level=1).shift(tau)
    pd.concat([df, df_tau], axis=1, keys=['df', 'df_'+str(tau)]).groupby(level=0).corr(method='pearson').mean(level=(1,2)).loc['df', 'df_'+str(tau)][columns_]

    # univariate portfolio analysis
    # 一个疑问：ivol与size之间是负相关的，因此，直观感觉上，VW会强化ivol puzzle，然而我的结果为什么是弱化ivol puzzle?
    # 如果剔除size中最小的30%，结果会如何？
    import_pv_index()
    periods=[1,3,6,12]
    periods_y=[12,24,36]

    SMB,HML=cal_SMB_HML(freq='D')
    SMB_m, HML_m = cal_SMB_HML(freq='M')
    ret=cal_rev(del_rf=True)
    size=cal_size(freq='M')
    BM=cal_BM(size,freq='M')
    index_ret_m = cal_index_ret(freq='M')
    fSMB = SMB_m.iloc[:, 0] - SMB_m.iloc[:, -1]
    fHML = HML_m.iloc[:, 0] - HML_m.iloc[:, -1]

    res_mean = pd.DataFrame(index=periods, columns=range(1, 13))
    res_t = pd.DataFrame(index=res_mean.index,columns=res_mean.columns)
    res_y_mean = pd.DataFrame(index=periods_y, columns=range(1, 13))
    res_y_t = pd.DataFrame(index=res_y_mean.index, columns=res_y_mean.columns)
    weights= weights=size.shift(1) # None #
    W='EW' if weights is None else 'VW'
    for i in periods:
        ivol=cal_ivol(i,SMB,HML,method='FF')
        tmp=cal_mimick_port1(ivol.loc['2005-'+str(i):'201802'],ret.loc['2005-'+str(i):'201802'],weights,10)
        tmp[11]=tmp.iloc[:,-1]-tmp.iloc[:,0]
        #tmp.mean()/NWest_mean(tmp)
        res_mean.loc[i]=tmp.mean()
        res_t.loc[i]=res_mean.loc[i]/NWest_mean(tmp)
        a, b=cal_FF_alpha(tmp[11],i,index_ret_m,fSMB,fHML)
        res_mean.loc[i, 12]=a
        res_t.loc[i,12]=a/b
    for i in periods_y:
        ivol=cal_ivol_year(i,ret,index_ret_m,SMB_m,HML_m,method='FF')
        tmp=cal_mimick_port1(ivol.loc['200'+str(4+int(i/12))+'12':'201802'],ret.loc['200'+str(4+int(i/12))+'12':'201802'],weights,10)
        tmp[11] = tmp.iloc[:, -1] - tmp.iloc[:, 0]
        # tmp.mean()/NWest_mean(tmp)
        res_y_mean.loc[i] = tmp.mean()
        res_y_t.loc[i] = res_y_mean.loc[i] / NWest_mean(tmp)
        a, b = cal_FF_alpha(tmp[11], i, index_ret_m, fSMB, fHML)
        res_y_mean.loc[i, 12] = a
        res_y_t.loc[i, 12] = a / b

    pd.concat((res_mean,res_y_mean)).to_csv(data_path+'uni_portfolio_analysis_'+W+'_mean.csv')
    pd.concat((res_t, res_y_t)).to_csv(data_path+'uni_portfolio_analysis_'+W+'_t_value.csv')

    # portfolio characteristics
    # 从结果看，似乎是BM、rev、skew、iskew驱动着ivol puzzle
    import_pv_index()
    ivol=cal_ivol(1,SMB,HML,method='FF')
    beta=cal_beta(12)
    size=cal_size(freq='M') # 波动率最小的三组，size很大；ivol后面七组之间的差距不是很大
    BM=cal_BM(size,freq='M').loc[size.index,size.columns]
    mom=cal_mom(12)
    rev=cal_rev()
    illiq=cal_illiq(12)
    skew=cal_skew(12)
    coskew=cal_coskew(12)
    iskew=cal_iskew(1)
    cha=(ivol,beta,size,BM,mom,rev,illiq,skew,coskew,iskew)
    port_cha=pd.DataFrame(index=range(len(cha)),columns=range(1,len(cha)+1))
    for i in range(len(cha)):
        port_cha.iloc[i]=cal_mimick_port1(ivol.loc['200601':'201802'],cha[i].shift(1).loc['200601':'201802'],None,10).mean()
    port_cha.to_csv(data_path+'Portfolio_Characteristics.csv')

    # bivariate portfolio analysis
    import_pv_index()
    SMB, HML = cal_SMB_HML(freq='D')
    ivol = cal_ivol(1, SMB, HML, method='FF')
    beta = cal_beta(12)
    size = cal_size(freq='M')  # 波动率最小的三组，size很大；ivol后面七组之间的差距不是很大
    BM = cal_BM(size, freq='M').loc[size.index, size.columns]
    mom = cal_mom(12)
    rev = cal_rev()
    illiq = cal_illiq(12)
    skew = cal_skew(12)
    coskew = cal_coskew(12)
    iskew = cal_iskew(1)
    cha = (ivol, beta, size, BM, mom, rev, illiq, skew, coskew, iskew)
    Return_m=pd.DataFrame(index=range(len(cha)),columns=range(1,7))
    Return_t=pd.DataFrame(index=Return_m.index,columns=Return_m.columns)
    FF_alpha = pd.DataFrame(index=Return_m.index, columns=Return_m.columns)
    FF_t = pd.DataFrame(index=Return_m.index, columns=Return_m.columns)
    SMB_m, HML_m = cal_SMB_HML(freq='M')
    index_ret_m = cal_index_ret(freq='M')
    fSMB = SMB_m.iloc[:, 0] - SMB_m.iloc[:, -1]
    fHML = HML_m.iloc[:, 0] - HML_m.iloc[:, -1]
    ret=cal_rev(del_rf=True)
    weights = None  #weights = size.loc['200512':'201802']  #
    W='EW' if weights is None else 'VW'
    independent = True  #independent = False  #
    inde = 'inde' if independent else 'depe'
    for i in range(len(cha)):
        tmp=cal_mimick_port2(cha[i].loc['200512':'201802'],ivol.loc['200512':'201802'],ret.loc['200512':'201802'],weights,independent=independent)
        long_short = (tmp.iloc[:, -1] - tmp.iloc[:, 0]).unstack()
        long_short[6]=long_short.mean(axis=1)
        Return_m.iloc[i,:]=long_short.mean()
        Return_t.iloc[i,:]=long_short.mean()/NWest_mean(long_short)
        for j in range(6):
            FF_alpha.iloc[i,j],se=cal_FF_alpha(long_short.iloc[:,j],12,index_ret_m,fSMB,fHML)
            FF_t.iloc[i,j]=FF_alpha.iloc[i,j]/se
    Return_m.to_csv(data_path+'Return_mean_'+W+'_'+inde+'.csv')
    Return_t.to_csv(data_path + 'Return_tValue_'+W+'_'+inde+'.csv')
    FF_alpha.to_csv(data_path+'FF_alpha_'+W+'_'+inde+'.csv')
    FF_t.to_csv(data_path+'FF_tValue_'+W+'_'+inde+'.csv')



    # Fama-MacBeth Regression
    import statsmodels.api as sm
    import_pv_index()
    beta=cal_beta(12).shift(1)
    size=cal_size(freq='M')
    BM=cal_BM(size,freq='M').loc[size.index,size.columns].shift(1);size=(np.log(size)).shift(1)

    SMB, HML = cal_SMB_HML(freq='D')
    ivol=cal_ivol(1,SMB,HML,method='CAPM').shift(1)
    mom=cal_mom(12).shift(1)
    rev=cal_rev().shift(1)
    illiq=cal_illiq(12).shift(1)
    skew=cal_skew(12).shift(1)
    coskew=cal_coskew(12).shift(1)
    iskew=cal_iskew(12,method='FF').shift(1)
    ret=cal_rev(del_rf=True)
    periods=size['200601':'201802'].index
    columns_ = ['ivol', 'beta','size','BM','illiq']  #columns_=['beta','size','BM','ivol','mom','rev','illiq','skew','coskew','iskew'] #
    params=pd.DataFrame(index=periods,columns=['const']+columns_)
    rsqurared_adj=pd.DataFrame(index=params.index,columns=['R2'])
    for t in periods:
        X=pd.DataFrame({'beta':beta.loc[t],'size':size.loc[t],'BM':BM.loc[t],'ret':ret.loc[t],'ivol':ivol.loc[t],'mom':mom.loc[t],'rev':rev.loc[t],'illiq':illiq.loc[t],'skew':skew.loc[t],'coskew':coskew.loc[t],'iskew':iskew.loc[t]})[columns_+['ret']]
        X=sm.add_constant(X).dropna()
        model_fit=sm.OLS(X.iloc[:,-1],X.iloc[:,:-1]).fit()
        params.loc[t]=model_fit.params
        rsqurared_adj.loc[t]=model_fit.rsquared_adj
    params.mean()
    params.mean()/NWest_mean(params)
    rsqurared_adj.mean()


    # conditional volatility from EGARCH
    from arch import arch_model
    from time import time
    import_pv_index(del_limit_move=False)
    SMB_m,HML_m=cal_SMB_HML(freq='M')
    fSMB = SMB_m.iloc[:, 0] - SMB_m.iloc[:, -1]
    fHML = HML_m.iloc[:, 0] - HML_m.iloc[:, -1]
    ret=cal_rev(del_rf=True,first_day=False)
    index_ret_m=cal_index_ret(freq='M')

    #y=ret.iloc[1:60]
    X=pd.DataFrame({'index':index_ret_m,'SMB':fSMB,'HML':fHML})[['index','SMB','HML']]
    e_ivol=pd.DataFrame(index=ret.index,columns=ret.columns)
    %%timeit
    for i in ret.index[31:32]:
        tmp1 = X.loc['200502':i]
        tmp2 = ret.loc['200502':i]
        XY = pd.DataFrame(
            {'MKT': ((tmp1['index'] - tmp1['index'].mean()).values[:, None] * (tmp2 - tmp2.mean())).mean(),
             'SMB': ((tmp1['SMB'] - tmp1['SMB'].mean()).values[:, None] * (tmp2 - tmp2.mean())).mean(),
             'HML': ((tmp1['HML'] - tmp1['HML'].mean()).values[:, None] * (tmp2 - tmp2.mean())).mean()}).T
        beta_tmp = np.linalg.pinv((tmp1 - tmp1.mean()).cov()) @ (XY.loc[['MKT', 'SMB', 'HML']])
        e = tmp2 - tmp1 @ beta_tmp;
        e -= e.mean()
        for col in e.columns:
            ee=e[col][e[col].isnull()[::-1].idxmax():]
            if len(ee)-1 >= 30:
                e_ivol.loc[i, col] = arch_model(ee[1:], mean='zero', vol='egarch', p=1, o=1, q=1).fit(disp='off',
                                                                                                 show_warning=False).conditional_volatility[-1]
            else:
                e_ivol.loc[i, col] = np.nan



    arch_model(y,mean='zero',vol='egarch',p=1,o=1,q=1).fit(disp='off',show_warning=False)#.conditional_volatility[-10:]
    periods=12
    cvol = pd.DataFrame(index=pd.date_range('20050101', '20180301', freq='M'), columns=rtn.columns[:-2])
    t0=time()
    for i in range(periods, 158):
        tmp = rtn[(rtn['month'] >= i - periods + 1) & (rtn['month'] <= i)]
        for c in cvol.columns:
            if (~tmp[c][-100].isnull()).sum()>100.0:
                cvol.iloc[i - 1][c] =arch_model(tmp[c].fillna(0),p=1,o=1,q=1,vol='egarch').fit(disp='off',show_warning=False).conditional_volatility[-1]
            else:
                cvol.iloc[i - 1][c] =np.nan
    cvol.to_csv(data_path+'conditonal_volatility_12M_EGARCH.csv')
    time()-t0





    arch_model(index_ret['pctchange'], p=1, o=1, q=1, vol='egarch').fit().conditional_volatility[-10:]
    arch_model(rtn['000001.SZ'].combine_first(index_ret['pctchange'])[:100].fillna(0), p=1, o=1, q=1, vol='egarch').fit().conditional_volatility[-10:]
    arch_model(rtn['000001.SZ'][:100].fillna(0), p=1, o=1, q=1, vol='egarch').fit(disp='off').conditional_volatility[-10:]


















X = pd.DataFrame({'index': index_ret_m, 'SMB': fSMB, 'HML': fHML})[['index', 'SMB', 'HML']].loc[start[12]:'201802']
X.insert(0,'c',np.ones(len(X)))
%timeit sm.OLS(long_short[1],X).fit().params
a,b=cal_FF_alpha(long_short.iloc[:,0],12,index_ret_m,fSMB,fHML);a/b


    import_pv_index()
    ivol_FF=cal_ivol(12,method='FF')
    vol=cal_vol(1)
    vol_ss=cal_vol_ss(1)
    ret=cal_rev(del_rf=True)
    SMB,HML=cal_SMB_HML(freq='D')
    tmp=cal_ivol(12,SMB,HML,method='FF')
    size=cal_size(freq='M')
    ivol_mimick=cal_mimick_port1(illiq.loc['200512':'201802'],ret.loc['200512':'201802'],None)#size.loc['200502':'201802'].shift(1))
    tmp=ivol_mimick.iloc[:,-1]-ivol_mimick.iloc[:,0];tmp.mean()/tmp.std()*np.sqrt(len(tmp))


import numpy as np
from arch import arch_model
from numba import jit
from numpy.linalg import pinv
%timeit arch_model(rtn.iloc[:200,0],p=1,q=1).fit()
a=np.random.randn(100,1000)
b=np.random.randn(10)
arr=np.random.randn(300);arr-=arr.mean()
arch_model(arr,p=1,q=1,vol='egarch').fit(disp='off',show_warning=False)
%timeit res=arch_model(arr-arr.mean(),p=1,o=1,q=1,vol='egarch').fit(disp='off',show_warning=False).conditional_volatility[-1]

#@jit
def garch1(arr,max_iter,init=(0.1,0.4,0.5)):
    theta=np.array(init)
    T = len(arr)

    for i in range(1,max_iter):
        punish = np.array((0, -0.2, -0.2)) * np.maximum(theta[1] + theta[2] - 1.0, 0)
        G = 0.0;
        J = 0.0
        sigma2_lag = arr[0] * arr[0]
        for t in range(1, T):
            tmp = np.array((1.0, arr[t-1] * arr[t-1], sigma2_lag))
            sigma2 = theta[0] + theta[1] * arr[t - 1]* arr[t - 1] + theta[2] * sigma2_lag
            sigma2_lag=sigma2
            g = 0.5*(arr[t] * arr[t] / sigma2 / sigma2-1.0 / sigma2) * tmp#+punish
            G += g
            J +=0.5/sigma2/sigma2*tmp[:,None]*tmp
            #H -=(0.5 / sigma2 / sigma2 - arr[t] * arr[t] / sigma2/ sigma2/ sigma2) * tmp[:,None]*tmp #利用信息矩阵：0.5 / sigma2 / sigma2 * h
        G /=(T-1)
        J /=(T-1)
        theta +=pinv(J)@G
    return theta

from modules.uni_tsa import uni_tsa
arr=np.random.randn(300)
%timeit uni_tsa(arr).garch()