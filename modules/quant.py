import pandas as pd
import numpy as np
from pandas.tseries.offsets import MonthEnd,YearEnd,Week,Day
#from dateutil.parser import parse
import sys
import warnings
warnings.filterwarnings("ignore")
#data_path ='E:/data/NewData/'  #'/Users/harbes/data/NewData/'#
#data_path=_data_path()
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
        date_list=data.index.where(data.index == ((data.index + MonthEnd()) - MonthEnd()),
                                      data.index + MonthEnd())
        #date_list=pd.to_datetime(date_list.astype(str),format='%Y%m')+MonthEnd()
    elif freq=='W':
        week_day = 5  # 0-6分别对应周一至周日
        date_list = data.index.where(data.index == ((data.index + Week(weekday=week_day)) - Week()),
                                      data.index + Week(weekday=week_day))

        #date_list=pd.to_datetime(date_list.astype(str).str.pad(7,side='right',fillchar='6'),format='%Y%W%w')
    elif freq=='Y':
        date_list = data.index.where(data.index == ((data.index + YearEnd()) - YearEnd()),
                                      data.index + YearEnd())
        #date_list = pd.to_datetime(date_list.astype(str), format='%Y') + YearEnd()
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
        data.index = data.index.where(data.index == ((data.index + MonthEnd()) - MonthEnd()),
                                      data.index + MonthEnd())
    elif to_freq=='W':
        # By=lambda x:x.year*100+x.week # 此种方法转化为周末日期时会出现错误
        week_day=5 #0-6分别对应周一至周日
        data.index=data.index.where(data.index==((data.index+Week(weekday=week_day))-Week()),data.index+Week(weekday=week_day))

    elif to_freq=='Y':
        data.index = data.index.where(data.index == ((data.index + YearEnd()) - YearEnd()),
                                      data.index + YearEnd())
    By = lambda x: x
    tmp=data.groupby(By)
    #data.resample('D',fill_method='pad')
    if by_position:
        data=tmp.nth(n_th)
    elif n_th:
        data=tmp.last()
    else:
        data =tmp.first()
    #if to_freq=='W':
    #    data.drop(data.index.duplicated(keep=last),inplace=True)
    #if to_freq=='M':
    #    data.index=pd.to_datetime(data.index.astype(str),format='%Y%m')+MonthEnd()
    #elif to_freq=='W':
        # 转成每周的固定day，0对应周日，1-6分别对应星期一到星期六
        # data.index=pd.to_datetime(data.index.astype(str).str.pad(7,side='right',fillchar='6'),format='%Y%W%w')
    #elif to_freq=='Y':
    #    data.index=pd.to_datetime(data.index.astype(str),format='%Y')+YearEnd()
    return data
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
        percentile2 = [0.0, 0.3, 0.7, 1.0]  # value
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
    valid_ = ~(pd.isnull(mark_1 + mark_2) | pd.isnull(ret.iloc[1:]))  # valid的股票要满足：当期有前一个月的indicator信息；当期保证交易
    df = pd.DataFrame()
    df['rtn'] = ret.iloc[1:][valid_].stack()
    df['ref1'] = mark_1[valid_].stack()
    df['ref2'] = mark_2[valid_].stack()
    tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).mean()).unstack()
    tmp.columns = tmp.columns.get_level_values(1)
    tmp.index.names = ('date', 'ref1')
    HML = tmp.mean(axis=0, level=0)
    SMB = tmp.mean(axis=1).unstack()
    return SMB.iloc[:,-1]-SMB.iloc[:,0], HML.iloc[:,-1]-HML.iloc[:,0]
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
                 -ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
            beta.loc[edt]=tmp2.mul(tmp1,axis=0).mean()/tmp1.var()
    elif freq=='W':
        for edt in EndDate_list[periods-1:]:
            tmp1=index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt]\
                 -index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean()
            tmp2=ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt]\
                 -ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean()
            beta.loc[edt]=tmp2.mul(tmp1,axis=0).mean()/tmp1.var()
    #elif freq=='Y':
    #    for edt in EndDate_list[periods-1:]:
    #        tmp1=index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(years=periods) + Day():edt]\
    #             -index_ret_d.loc[edt - pd.tseries.offsets.DateOffset(years=periods) + Day():edt].mean()
    #        tmp2=ret_d.loc[edt - pd.tseries.offsets.DateOffset(years=periods) + Day():edt]\
    #             -ret_d.loc[edt - pd.tseries.offsets.DateOffset(years=periods) + Day():edt].mean()
    #        beta.loc[edt]=tmp2.mul(tmp1,axis=0).mean()/tmp1.var()
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
def cal_iskew(periods,freq='M',method='',ret=None,index_ret=None,SMB=None,HML=None):
    if ret is None:
        ret=cal_ret(freq=freq)
        index_ret=cal_index_ret(freq=freq)
        if method=='FF':
            SMB,HML=cal_SMB_HML(freq=freq)
    EndDate_list=_GetEndDateList(ret,freq)
    iskew = pd.DataFrame(index=EndDate_list, columns=ret.columns)
    if method == 'CAPM':
        if freq == 'M':
            for edt in EndDate_list[periods - 1:]:
                tmp1 = index_ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt] \
                       - index_ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
                tmp2 = ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt] \
                       - ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
                iskew.loc[edt] =(tmp2-tmp1[:,None]@tmp2.mul(tmp1, axis=0).mean()[None,:]/tmp1.var()).skew()
        elif freq == 'W':
            for edt in EndDate_list[periods - 1:]:
                tmp1 = index_ret.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt] \
                       - index_ret.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean()
                tmp2 = ret.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt] \
                       - ret.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean()
                iskew.loc[edt] = (tmp2 - tmp1[:, None] @ tmp2.mul(tmp1, axis=0).mean()[None, :] / tmp1.var()).skew()
    elif method=='FF':
        XY=pd.concat((pd.DataFrame({'index':index_ret,'SMB':SMB,'HML':HML}),ret),axis=1)
        if freq == 'M':
            for edt in EndDate_list[periods - 1:]:
                tmp1 = XY.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt] \
                       - XY.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
                #tmp2 = ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt] \
                #       - ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
                iskew.loc[edt] =(tmp1[ret.columns]-tmp1[['index','SMB','HML']]@np.linalg.pinv(tmp1[['index','SMB','HML']].cov())@np.array([
                    tmp1[ret.columns].mul(tmp1['index'],axis=0).mean(),
                    tmp1[ret.columns].mul(tmp1['SMB'], axis=0).mean(),
                    tmp1[ret.columns].mul(tmp1['HML'], axis=0).mean()
                ])).skew()
        elif freq == 'W':
            for edt in EndDate_list[periods - 1:]:
                tmp1 = XY.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt] \
                       - XY.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean()
                #tmp2 = ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt] \
                #       - ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
                iskew.loc[edt] =(tmp1[ret.columns]-tmp1[['index','SMB','HML']]@np.linalg.pinv(tmp1[['index','SMB','HML']].cov())@np.array([
                    tmp1[ret.columns].mul(tmp1['index'],axis=0).mean(),
                    tmp1[ret.columns].mul(tmp1['SMB'], axis=0).mean(),
                    tmp1[ret.columns].mul(tmp1['HML'], axis=0).mean()
                ])).skew()
    return iskew.iloc[periods-1:].astype(float)
def cal_mimick_port1(indi,freq='M',ret=None,weights=None,percentile=None):
    if percentile is None:
        percentile=np.arange(0.0,1.01,0.1)
    label_=[i for i in range(1,len(percentile))]
    if ret is None:
        ret=cal_ret(freq=freq)
    if len(indi)!=len(ret):
        raise ValueError('the index of Indicator does not match the return data!')
    mark_ = pd.DataFrame([pd.qcut(indi.iloc[i],q=percentile, labels=label_) for i in range(len(indi)-1)],index=indi.index[1:])
    valid_=~(pd.isnull(mark_) | pd.isnull(ret.iloc[1:]))
    if weights is None:
        df = pd.DataFrame()
        df['rtn'] = ret.iloc[1:][valid_].stack()
        df['ref'] = mark_[valid_].stack()
        tmp=df.groupby(level=0).apply(lambda g: g.groupby('ref').mean()).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    else:
        df1=pd.DataFrame()
        df2=pd.DataFrame()
        df1['ret_w'] = (ret*weights).iloc[1:][valid_].stack()
        df1['ref'] = mark_[valid_].stack()
        df2['ret_w'] = weights.iloc[1:][valid_].stack()
        df2['ref'] = mark_[valid_].stack()
        tmp1=df1.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp2 = df2.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp=(tmp1/tmp2).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    return tmp
def cal_mimick_port2(indi1,indi2,freq='M',ret=None,weights=None,percentile1=None,percentile2=None,independent=True):
    if percentile1 is None:
        percentile1=np.arange(0.0,1.01,0.2)
        percentile2=np.arange(0.0,1.01,0.2)
    label_1 = [i for i in range(1, len(percentile1))]
    label_2 = [i for i in range(1, len(percentile2))]
    if len(indi1) != len(ret) or len(indi2) != len(ret):
        raise ValueError('the index of Indicators does not match the return data!')
    if independent:
        mark_1 = pd.DataFrame([pd.qcut(indi1.iloc[i], q=percentile1, labels=label_1) for i in range(len(indi1) - 1)],
                              index=indi1.index[1:])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
        mark_2 = pd.DataFrame([pd.qcut(indi2.iloc[i], q=percentile2, labels=label_2) for i in range(len(indi2) - 1)],
                              index=indi2.index[1:])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
    else:
        mark_1 = pd.DataFrame([pd.qcut(indi1.iloc[i], q=percentile1, labels=label_1) for i in range(len(indi1) - 1)],
                              index=indi1.index[1:])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
        mark_2=pd.DataFrame(index=mark_1.index,columns=mark_1.columns)
        for l_ in label_1:
            tmp=pd.DataFrame([pd.qcut(indi2.iloc[i][mark_1.iloc[i]==l_],q=percentile,labels=label_) for i in range(len(indi2)-1)],index=mark_1.index)
            mark_2 = mark_2.combine_first(tmp)
    valid_ = ~(pd.isnull(mark_1+mark_2) | pd.isnull(ret.iloc[1:]))  # valid的股票要满足：当期有前一个月的indicator信息；当期保证交易
    if weights is None:
        df = pd.DataFrame()
        df['rtn'] = ret.iloc[1:][valid_].stack()
        df['ref1'] = mark_1[valid_].stack()
        df['ref2'] = mark_2[valid_].stack()
        tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).mean()).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    else:
        df1 = pd.DataFrame()
        df2 = pd.DataFrame()
        df1['rtn_w'] = (ret * weights.shift(1)).iloc[1:][valid_].stack()
        df1['ref1'] = mark_1[valid_].stack()
        df1['ref2'] = mark_2[valid_].stack()
        df2['rtn_w'] = weights.shift(1).iloc[1:][valid_].stack()
        df2['ref1'] = mark_1[valid_].stack()
        df2['ref2'] = mark_2[valid_].stack()
        tmp1 = df1.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).sum())
        tmp2 = df2.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).sum())
        tmp = (tmp1 / tmp2).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    return tmp
def describe(df,stats=['skew','kurt']):
    d=df.describe(percentiles=[0.05,0.25,0.5,0.75,0.95])
    return d.append(df.reindex(d.columns, axis=1).agg(stats))
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
def cal_FF_alpha(arr,freq='M',index_ret=None, SMB=None, HML=None):
    if SMB is None:
        size=cal_size(freq=freq)
        BM=cal_BM(size,freq=freq)
        SMB,HML=cal_SMB_HML(freq=freq,size=size,BM=BM)
    #start={1:'200502',3:'200504',6:'200507',12:'200601',24:'200701',36:'200801'} # 回归有效数据的起始月
    X= pd.DataFrame({'index':index_ret,'SMB':SMB, 'HML':HML})[['index', 'SMB', 'HML']]#.loc[start[period]:'201802']
    #arr=arr.loc[start[period]:'201802']
    tmp1=X-X.mean(axis=0)
    tmp2=arr-arr.mean(axis=0)
    tmp3=arr-X@np.linalg.pinv(tmp1.cov())@(tmp1.mul(tmp2,axis=0).mean()).values
    return tmp3.mean(),NWest(tmp3-tmp3.mean(),X)


def cal_weekday_ret():
    pass





tmp.rolling('1M').std()
(tmp2-tmp2.mean()).mul(tmp1-tmp1.mean(),axis=0).mean()
(tmp2-tmp2.mean()).mul(tmp1-tmp1.mean(),axis=0)
tmp1.shape
%timeit ((tmp1-tmp1.mean()).values[:,None]*(tmp2-tmp2.mean())).mean()/tmp1.var()
%timeit ((tmp1-tmp1.mean()).values[:,None]*(tmp2-tmp2.mean())).mean()/((tmp1-tmp1.mean())**2).mean()




import statsmodels.api as sm
tmp=tmp1.copy()
tmp=sm.add_constant(tmp)
tmp['y']=tmp2['000002.SZ']
tmp=tmp.dropna().astype(float)
sm.OLS(tmp['y'],tmp[['const','index']]).fit().summary()


data_dict={'pv':['']}

