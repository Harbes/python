import pandas as pd
import numpy as np
from pandas.tseries.offsets import MonthEnd,MonthBegin,Day
from dateutil.parser import parse
import sys
import warnings
warnings.filterwarnings("ignore")
data_path ='E:/data/NewData/'  #'/Users/harbes/data/NewData/'#

def _data_path():
    sys_platform=sys.platform
    if sys_platform =='windows':
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
    if freq=='M':
        By = lambda x: int(x * 0.01)
        index_cls = index_ret['clsprc'].groupby(By).last()
        index_opn = index_ret['opnprc'].groupby(By).first()
        index_ret = (index_cls - index_opn) / index_opn * 100.0
        index_ret.index = pd.to_datetime(index_ret.index.astype(str), format='%Y%m') + MonthEnd()
    else:
        index_ret=(index_ret['clsprc']-index_ret['opnprc'])/index_ret['opnprc']*100.0
        index_ret.index = pd.to_datetime(index_ret.index.astype(int).astype(str), format='%Y%m%d')
    index_ret.index.name = 'trddt'
    if del_Rf:
        Rf=import_data(Rf_freq=freq)[2]
        #index_ret = index_ret.sub(Rf.loc[index_ret.index], axis=0) # 报错，Rf是一个DataFrame，而不是Series
        #index_ret=-Rf.loc[index_ret.index].sub(index_ret,axis=0) # 虽然成功，但是返回的是一个DataFrame，而不是Series
        #index_ret = index_ret.sub(Rf.loc[index_ret.index].iloc[:,0], axis=0)
        index_ret=(index_ret + 100.0) / (Rf.loc[index_ret.index].iloc[:, 0] * 0.01 + 1.0) - 100.0

    return index_ret
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
    price0=PV[p0].unstack()
    price1=PV[p1].unstack()
    if freq=='M':
        By=lambda x:x.year*100+x.month
        price0=price0.groupby(By).nth(0).shift(periods-1)
        price1=price1.groupby(By).nth(-1)
        ret=(price1-price0)/price0*100.0
        ret.index=pd.to_datetime(ret.index.astype(str),format='%Y%m')+MonthEnd()
    else:
        price0=price0.shift(periods-1)
        ret=(price1-price0)/price0*100.0
    if del_Rf:
        Rf = import_data(Rf_freq=freq)[2]
        ret = (ret.astype(float) + 100.0).div((Rf.loc[ret.index].iloc[:, 0] * 0.01 + 1.0) ** periods, axis=0) - 100.0
    return ret

ret_d=cal_ret(freq='D')
index_ret_d=cal_index_ret(freq='D')

def cal_variables(var_list,start='20050101',end='20180228'):
    pass

def cal_beta(periods,index_ret_d=None,ret_d=None):
    '''
    calculate the market beta
    :param periods: int{measure periods;number of month}
    :param save_data:
    :return:
    '''
    start_='20050101';end_='20180301'
    if index_ret_d is None or ret_d is None:
        ret=cal_ret(freq='D')
        stock_pool=ret.columns
        beta = pd.DataFrame(index=pd.date_range(start_, end_, freq='M'), columns=stock_pool)
        ret['index']=cal_index_ret(freq='D')
    else:
        ret=ret_d.copy()
        stock_pool = ret.columns
        beta = pd.DataFrame(index=pd.date_range(start_, end_, freq='M'), columns=stock_pool)
        ret['index']=index_ret_d.copy()
    ret['month']=(ret.index.year-ret.index.year[0])*12+ret.index.month
    for i in range(periods, ret['month']):
        # 也可以通过cov/var，但是速度较慢
        tmp1=ret['index'][(ret['month'] >= i - periods+1) & (ret['month'] <= i)]
        tmp2=ret[stock_pool][(i - periods+1 <= ret['month']) & (ret['month'] <= i)]
        year_=str(tmp1.index[-1].year);month_=str(tmp1.index[-1].month)
        beta.loc[year_+'-'+month_] = ((tmp2-tmp2.mean()).mul(tmp1-tmp1.mean(),axis=0).mean()/((tmp1-tmp1.mean())**2).mean())[None,:]
    return beta[beta!=0]
def cal_size(freq='M'):
    size=import_data(PV_vars=['size_tot'])[0].unstack()
    size.columns=size.columns.get_level_values(1)
    size=_resample_h2l(size,to_freq=freq,n_th=-1)
    return size[size>0.0]
def cal_BM(size=None,freq='M'):
    if size is None:
        size=cal_size(freq=freq)
    tmp=import_data(BS_vars=['tot_assets', 'tot_liab'])[1]
    book=tmp['tot_assets']-tmp['tot_liab']
    book = book.drop(book.index[book.index.duplicated(keep='last')]).unstack()
    book = book.resample('D').first().ffill()
    book=book.loc[size.index,size.columns] * 1e-8
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
    turnover=amount_/size_free*0.1
    turnover=turnover[turnover>0.0001]
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
def cal_MaxRet(periods=1,freq='M',ret=None):
    if ret is None:
        ret=cal_ret(freq='D',del_Rf=False)
    By=_GroupBy(freq)
    max_ret = ret.groupby(By).max().rolling(periods).max()
    if freq == 'M':
        max_ret.index = pd.to_datetime(max_ret.index.astype(str), format='%Y%m') + MonthEnd()
    elif freq == 'W':
        max_ret.index = pd.to_datetime(max_ret.index.astype(str).str.pad(7, side='right', fillchar='6'),
                                        format='%Y%W%w')
    return max_ret.iloc[periods-1:]
def cal_skew(periods,freq='M',ret=None):
    if ret is None:
        ret=cal_ret(freq='D')
    EndDate_list=_GetEndDateList(ret,freq)
    if freq=='M':
        return pd.DataFrame((ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].skew() for edt in EndDate_list)
                 , index=EndDate_list)
    elif freq=='W':
        return pd.DataFrame((ret.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].skew() for edt in EndDate_list)
                 , index=EndDate_list)
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




def cal_weekday_ret():
    pass


    if freq=='M':
        ret['period']=(ret.index.year-ret.index[0].year)*12+ret.index.month
    elif freq=='W':
        ret['period'] = (ret.index.year - ret.index[0].year) * 100 + ret.index.week



    for i in range(periods,159):
        t_skew.iloc[i-1]=rtn[(rtn['month']>=i-periods+1)&(rtn['month']<=i)].iloc[:,:-2].skew()
    if save_data:
        t_skew.astype(float).to_pickle(data_path+'skew_'+str(periods)+'M')
    else:
        return t_skew.astype(float)





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

