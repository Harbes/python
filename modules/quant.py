import pandas as pd
import numpy as np
from pandas.tseries.offsets import MonthEnd,MonthBegin,Day
from dateutil.parser import parse
import sys
import warnings
warnings.filterwarnings("ignore")
data_path ='E:/data/NewData/'  #'/Users/harbes/data/NewData/'#
sys.setrecursionlimit(100000)
threading.stack_size(200000000)
thread = threading.Thread(target=your_code)
thread.start()
def _data_path():
    sys_platform=sys.platform
    if sys_platform =='win32':
        return 'E:/data/NewData/'
    elif sys_platform=='mac':
        return '/Users/harbes/data/NewData/'
    elif sys_platform=='linux':
        return '/home/harbes/data/NewData/'
    else:
        raise ValueError('These is no such systerm in your work-station')
def _GroupBy(freq):
    if freq=='M':
        return lambda x:x.year*100+x.month
    elif freq=='W':
        return lambda x: x.year * 100 + x.week
def _GetEndDateList(data,freq):
    if freq=='M':
        date_list=data.index.year*100+data.index.month
        date_list=pd.to_datetime(date_list.astype(str),format='%Y%m')+MonthEnd()
    elif freq=='W':
        date_list = data.index.year * 100 + data.index.week
        date_list=pd.to_datetime(date_list.astype(str).str.pad(7,side='right',fillchar='6'),format='%Y%W%w')
    return sorted(set(date_list))
def _resample_h2l(data,to_freq='M',n_th=0,by_position=True):
    '''
    resample the given data(high freq to low freq), like daily to weekly or monthly
    :param data:
    :param freq:
        'M': month
        'W': week
    :param n_th:
        0,-1,...
    :param by_position:
    :return:
    '''
    if to_freq=='M':
        By=lambda x:x.year*100+x.month
    else:
        By=lambda x:x.year*100+x.week
    tmp=data.groupby(By)
    if by_position:
        data=tmp.nth(n_th)
    elif n_th:
        data=tmp.first()
    else:
        data =tmp.last()
    if to_freq=='M':
        data.index=pd.to_datetime(data.index.astype(str),format='%Y%m')+MonthEnd()
    else:
        # 转成每周的固定day，0对应周日，1-6分别对应星期一到星期六
        data.index=pd.to_datetime(data.index.astype(str).str.pad(7,side='right',fillchar='6'),format='%Y%W%w')
    return data

#pd.to_datetime(tmp.index.astype(str).str.pad(7,side='right',fillchar='0'),format='%Y%W%w')
def import_data(PV_vars=None, BS_vars=None,Rf_freq=None):
    '''

    :param pv_vars[pv数据已经指定index变量了]: list of strings
        'dretwd',
        'pre_close', 'opnprc', 'high', 'low', 'clsprc', 'amount',
        'adj_pre_close','adj_open', 'adj_high', 'adj_low', 'adj_close',
        'vwap',
        'adj_factor', 'adj_vwap',
        'size_tot', 'size_free'

    :param BS_vars[BS数据尚未指定index变量]: list of strings
        'fin_year','ann_dt', 'stkcd', 'tot_assets', 'tot_liab'

    :param Rf_freq:
        'Y','M','D'
        [Rf数据来源于csmar]

    :return:
    '''
    if PV_vars is None and BS_vars is None and Rf_freq is None:
        raise ValueError("You haven't name the varables to be imported")
    data_path=_data_path()
    if PV_vars is not None:
        PV=pd.read_pickle(data_path+'PV_datetime')[PV_vars]
    else:
        PV=None
    if BS_vars is not None:
        BS = pd.read_pickle(data_path + 'BS')[['ann_dt', 'stkcd']+BS_vars]  #
        BS = BS[~(BS['ann_dt'].isnull())]
        BS['ann_dt'] = pd.to_datetime(BS['ann_dt'].astype(int).astype(str), format='%Y%m%d')
        BS=BS.set_index(['ann_dt', 'stkcd']).sort_index()
    else:
        BS=None
    if Rf_freq is not None:
        if Rf_freq=='Y':
            Rf = pd.read_excel(data_path + 'risk_free.xlsx')[2:][['Clsdt','Nrrdata']].set_index(['Clsdt'])
        elif Rf_freq=='M':
            Rf = pd.read_excel(data_path + 'risk_free.xlsx')[2:][['Clsdt','Nrrmtdt']].set_index(['Clsdt'])
        elif Rf_freq=='W':
            Rf = (1.0+pd.read_excel(data_path + 'risk_free.xlsx')[2:][['Clsdt', 'Nrrdata']].set_index(['Clsdt'])*0.01)**(1.0/52.0)*100.0-100.0
        else:
            Rf = pd.read_excel(data_path + 'risk_free.xlsx')[2:][['Clsdt','Nrrdaydt']].set_index(['Clsdt'])
        Rf.index=pd.to_datetime(Rf.index,format='%Y-%m-%d')
    else:
        Rf=None
    return PV,BS,Rf

#np.log(1.0+Rf['Nrrdata'].astype(float)*0.01)/np.log(Rf['Nrrmtdt'].astype(float)*0.01+1.0)
def cal_index_ret(freq='M',index_code='399300.SZ',del_Rf=True):
    '''

    :param index_code:
        '000016.SH'：上证50
        '399005.SZ'：中小板指数
        '399006.SZ'：创业板指
        '399300.SZ'：沪深300指数
        '399905.SZ'：中证500指数
    :param del_Rf:
        delete risk-free or not
    :return:
    '''
    data_path=_data_path()
    index_ret = pd.read_pickle(data_path + 'index_ret')[['index_code', 'trddt', 'opnprc', 'clsprc']]
    index_ret = index_ret[index_ret['index_code'] == index_code][['trddt', 'opnprc', 'clsprc']].set_index('trddt')
    index_ret.index=pd.to_datetime(index_ret.index.astype(int).astype(str),format='%Y%m%d')
    price0=_resample_h2l(index_ret['opnprc'],to_freq=freq,n_th=0)
    price1=_resample_h2l(index_ret['clsprc'],to_freq=freq,n_th=-1)
    index_ret=(price1-price0)/price0*100.0
    index_ret.index.name = 'trddt'
    if del_Rf:
        Rf=import_data(Rf_freq=freq)[2]
        #index_ret = index_ret.sub(Rf.loc[index_ret.index], axis=0) # 报错，Rf是一个DataFrame，而不是Series
        #index_ret=-Rf.loc[index_ret.index].sub(index_ret,axis=0) # 虽然成功，但是返回的是一个DataFrame，而不是Series
        #index_ret = index_ret.sub(Rf.loc[index_ret.index].iloc[:,0], axis=0)
        index_ret=(index_ret + 100.0) / (Rf.loc[index_ret.index].iloc[:, 0] * 0.01 + 1.0) - 100.0

    return index_ret.astype(float)
def cal_ret(freq='M',periods=1,del_Rf=True):
    # TODO 尚未剔除异常交易数据
    '''

    calculate the past returns,including the current period
    :param periods:
    :param freq:
    :param del_Rf:
    :return:
    '''
    p0 = 'adj_open';p1 = 'adj_close' # price0也可以使用 'adj_pre_close'
    PV=import_data(PV_vars=[p0,p1])[0]
    price0=_resample_h2l(PV[p0].unstack(),to_freq=freq,n_th=0).shift(periods-1)
    price1=_resample_h2l(PV[p1].unstack(),to_freq=freq,n_th=-1)
    ret=(price1-price0)/price0*100.0
    if del_Rf:
        Rf = import_data(Rf_freq=freq)[2]
        ret = (ret + 100.0).div((Rf.loc[ret.index].iloc[:, 0] * 0.01 + 1.0) ** periods, axis=0) - 100.0
    return ret.iloc[periods-1:].astype(float)
def cal_variables(var_list,start='20050101',end='20180228'):
    pass
def cal_SMB_HML(freq='M',size=None,BM=None,percentile1=None,percentile2=None,independent=True):
    if size is None:
        size=cal_size(freq=freq)
        BM=cal_BM(size=size,freq=freq)
    ret=cal_ret(freq=freq)
    if percentile1 is None:
        percentile1 = [0.0, 0.5, 1.0]  # size
        #label_1 = ('S', 'B')
        percentile2 = [0.0, 0.3, 0.7, 1.0]  # value
        #label_2 = ('L', 'M', 'H')
    label_1=[i+1 for i in range(len(percentile1)-1)]
    label_2=[i+1 for i in range(len(percentile2)-1)]
    if independent:
        #mark_1 = pd.DataFrame([pd.qcut(size.iloc[i], q=percentile1, labels=label_1) for i in size.index[:-1]],
        #                      index=size.index[:-1]) # 报错
        mark_1 = pd.DataFrame([pd.qcut(size.iloc[i], q=percentile1, labels=label_1) for i in range(len(size)-1)],
                              index=size.index[1:])
        mark_2 = pd.DataFrame([pd.qcut(BM.iloc[i], q=percentile2, labels=label_2) for i in range(len(BM) - 1)],
                              index=BM.index[1:])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
    else:
        mark_1 = pd.DataFrame([pd.qcut(size.iloc[i], q=percentile1, labels=label_1) for i in range(len(size) - 1)],
                              index=size.index[1:])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
        mark_2=pd.DataFrame(index=mark_1.index,columns=mark_1.columns)
        for l_ in label_1:
            tmp=pd.DataFrame([pd.qcut(BM.iloc[i][mark_1.iloc[i]==l_],q=percentile2,labels=label_2) for i in range(len(BM)-1)],index=mark_1.index)
            mark_2 = mark_2.combine_first(tmp)
    valid_ = ~(pd.isnull(mark_1 + mark_2) | pd.isnull(ret[1:]))  # valid的股票要满足：当期有前一个月的indicator信息；当期保证交易
    ret[valid_].index
    df = pd.DataFrame()
    df['rtn'] = ret[valid_].stack()
    df['ref1'] = mark_1.stack()
    df['ref2'] = mark_2.stack()
    tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).mean()).unstack()
    tmp.columns = tmp.columns.get_level_values(1)
    tmp.index.names = ('date', 'ref1')
    # HML=tmp.loc[(slice(None), 'S'), :].reset_index(level=1,drop=True).sub(tmp.loc[(slice(None), 'B'), :].reset_index(level=1,drop=True))[['L','M','H']]
    HML = tmp.mean(axis=0, level=0)[['L', 'M', 'H']]
    SMB = tmp.mean(axis=1).unstack()
    return SMB, HML


def cal_beta(periods,freq='M',index_ret_d=None,ret_d=None):
    if index_ret_d is None:
        index_ret_d=cal_index_ret(freq='D')
        ret_d=cal_ret(freq='D')
    EndDate_list=_GetEndDateList(ret_d,freq)
    beta=pd.DataFrame(index=EndDate_list,columns=ret_d.columns)
    if freq=='M':
        for edt in EndDate_list[periods-1:]:
            tmp1=index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt]\
                 -index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
            tmp2=ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt]\
                 -index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
            beta.loc[edt]=tmp2.mul(tmp1,axis=0).mean()/tmp1.var()
    elif freq=='W':
        for edt in EndDate_list[periods-1:]:
            tmp1=index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt]\
                 -index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean()
            tmp2=ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt]\
                 -index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean()
            beta.loc[edt]=tmp2.mul(tmp1,axis=0).mean()/tmp1.var()
    return beta[beta!=0.0].iloc[periods-1:].astype(float)
def cal_size(freq='M'):
    size=import_data(PV_vars=['size_tot'])[0].unstack()
    size.columns=size.columns.get_level_values(1)
    size=_resample_h2l(size,to_freq=freq,n_th=-1,by_position=False)
    return size[size>0.0]
def cal_BM(size=None,freq='M'):
    if size is None:
        size=cal_size(freq=freq)
    tmp=import_data(BS_vars=['tot_assets', 'tot_liab'])[1]
    book=tmp['tot_assets']-tmp['tot_liab']
    book = book.drop(book.index[book.index.duplicated(keep='last')]).unstack()
    book = book.resample('D').first().ffill()
    book=book.loc[size.index,size.columns] * 1e-3
    BM=book/size
    return BM[BM>0.0]
def cal_mom(periods,freq='M',opn=None,cls=None):
    if opn is None:
        pv=import_data(PV_vars=['adj_open','adj_close'])[0]
        opn=pv['adj_open'].unstack()
        cls=pv['adj_close'].unstack()
    open_price=_resample_h2l(opn,to_freq=freq,n_th=0,by_position=False).shift(periods-1)
    close_price=_resample_h2l(cls,to_freq=freq,n_th=-1,by_position=False).shift(1)
    return ((close_price-open_price)/open_price*100.0).iloc[periods-1:]
def cal_rev(periods=1,freq='M',opn=None,cls=None):
    if opn is None:
        pv=import_data(PV_vars=['adj_open','adj_close'])[0]
        opn=pv['adj_open'].unstack()
        cls=pv['adj_close'].unstack()
    open_price=_resample_h2l(opn,to_freq=freq,n_th=0,by_position=False).shift(periods-1)
    close_price=_resample_h2l(cls,to_freq=freq,n_th=-1,by_position=False)
    return ((close_price-open_price)/open_price*100.0).iloc[periods-1:]
def cal_illiq(periods=12,freq='M'):
    pv=import_data(PV_vars=['adj_open','adj_close','amount'])[0]
    amount = pv['amount'].unstack()
    open_price=pv['adj_open'].unstack()
    close_price=pv['adj_close'].unstack()
    illiq_d=1e5*np.abs(close_price - open_price) / open_price / amount
    if freq=='M':
        By = lambda x: x.year * 100 + x.month
    elif freq=='W':
        By = lambda x: x.year * 100 + x.week
    illiq_sum = illiq_d.groupby(By).sum()
    illiq_count = illiq_d.groupby(By).count()
    illiq_sum=illiq_sum.rolling(periods).sum()
    illiq_count=illiq_count.rolling(periods).sum()
    illiq=illiq_sum/illiq_count
    if freq=='M':
        illiq.index=pd.to_datetime(illiq.index.astype(str),format='%Y%m')+MonthEnd()
    elif freq=='W':
        illiq.index = pd.to_datetime(illiq.index.astype(str).str.pad(7,side='right',fillchar='6'),format='%Y%W%w')
    return illiq.iloc[periods-1:]

    _sum = 0;
    _count = 0;
    for i in range(periods):
        _sum += rtn_amount_sum.shift(i)
        _count += rtn_amount_count.shift(i)
    illiq = _sum / _count
    illiq.index = pd.to_datetime(illiq.index.astype(str), format='%Y%m') + MonthEnd()
    illiq.iloc[:periods - 1] = np.nan
    if save_data:
        illiq.to_pickle(data_path + 'illiq_' + str(periods) + 'M')
    else:
        return illiq
def cal_turnover(periods,freq='M',amount=None,size_free=None):
    if amount is None:
        pv=import_data(PV_vars=['amount','size_free'])[0]
        amount_=pv['amount'].unstack()
        size_free=pv['size_free'].unstack()
    turnover=amount_/size_free*10.0
    #turnover=turnover[turnover>0.0001]
    if freq=='M':
        By=lambda x:x.year*100+x.month
    elif freq=='W':
        By = lambda x: x.year * 100 + x.week
    turnover_sum=turnover.groupby(By).sum()
    turnover_count=turnover.groupby(By).count();turnover_count[turnover_count==0.0]=np.nan
    turnover=turnover_sum.rolling(periods).sum()/turnover_count.rolling(periods).sum()
    if freq == 'M':
        turnover.index = pd.to_datetime(turnover.index.astype(str), format='%Y%m') + MonthEnd()
    elif freq == 'W':
        turnover.index = pd.to_datetime(turnover.index.astype(str).str.pad(7, side='right', fillchar='6'), format='%Y%W%w')
    return turnover.iloc[periods-1:]
def cal_MaxRet(periods=1,freq='M',ret_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D',del_Rf=False)
    By=_GroupBy(freq)
    max_ret = ret_d.groupby(By).max().rolling(periods).max()
    if freq == 'M':
        max_ret.index = pd.to_datetime(max_ret.index.astype(str), format='%Y%m') + MonthEnd()
    elif freq == 'W':
        max_ret.index = pd.to_datetime(max_ret.index.astype(str).str.pad(7, side='right', fillchar='6'),
                                        format='%Y%W%w')
    return max_ret.iloc[periods-1:]
def cal_skew(periods,freq='M',ret_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D')
    EndDate_list=_GetEndDateList(ret_d,freq)
    if freq=='M':
        return pd.DataFrame((ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].skew() for edt in EndDate_list)
                 , index=EndDate_list).iloc[periods-1]
    elif freq=='W':
        return pd.DataFrame((ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].skew() for edt in EndDate_list)
                 , index=EndDate_list).iloc[periods-1]
def cal_coskew(periods,freq='M',ret_d=None,index_ret_d=None):
    if index_ret_d is None:
        ret_d=cal_ret(freq='D')
        index_ret_d=cal_index_ret(freq='D')
    X=pd.DataFrame({'index':index_ret_d,'index^2':index_ret_d**2.0})
    EndDate_list=_GetEndDateList(ret_d,freq)
    coskew=pd.DataFrame(index=EndDate_list,columns=ret_d.columns)
    if freq=='M':
        for edt in EndDate_list[periods - 1:]:
            tmp1=X.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt]\
                 -X.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
            tmp2=ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt]\
                 -ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
            tmp3 = np.linalg.pinv(tmp1.cov())[1]
            coskew.loc[edt] = tmp3[0] * tmp2.mul(tmp1['index'], axis=0).mean() + tmp3[1] * tmp2.mul(tmp1['index^2'],
                                                                                                    axis=0).mean()
    elif freq=='W':
        for edt in EndDate_list[periods - 1:]:
            tmp1=X.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt]\
                 -X.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean()
            tmp2=ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt]\
                 -ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean()
            tmp3=np.linalg.pinv(tmp1.cov())[1]
            coskew.loc[edt]=tmp3[0]*tmp2.mul(tmp1['index'],axis=0).mean()+tmp3[1]*tmp2.mul(tmp1['index^2'],axis=0).mean()
    return coskew.iloc[periods-1:].astype(float)








def cal_weekday_ret():
    pass





tmp.rolling('1M').std()
(tmp2-tmp2.mean()).mul(tmp1-tmp1.mean(),axis=0).mean()
(tmp2-tmp2.mean()).mul(tmp1-tmp1.mean(),axis=0)
tmp1.shape
%timeit ((tmp1-tmp1.mean()).values[:,None]*(tmp2-tmp2.mean())).mean()/tmp1.var()
%timeit ((tmp1-tmp1.mean()).values[:,None]*(tmp2-tmp2.mean())).mean()/((tmp1-tmp1.mean())**2).mean()

    if save_data:
        beta[beta!=0].astype(float).to_pickle(data_path + 'beta_daily_'+str(periods)+'M')
    else:
        return beta[beta!=0].astype(float).iloc[:,:-1] # 剔除 column：index


import statsmodels.api as sm
tmp=tmp1.copy()
tmp=sm.add_constant(tmp)
tmp['y']=tmp2['000002.SZ']
tmp=tmp.dropna().astype(float)
sm.OLS(tmp['y'],tmp[['const','index']]).fit().summary()


data_dict={'pv':['']}

