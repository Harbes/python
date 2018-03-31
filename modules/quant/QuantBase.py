import pandas as pd
import numpy as np
from pandas.tseries.offsets import MonthEnd,YearEnd,Week,Day,DateOffset
import statsmodels.api as sm
#from dateutil.parser import parse
import sys
import warnings
warnings.filterwarnings("ignore")
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
def _resample_index(data,to_freq):
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
    return data
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
        if False:
            -1: the last non-NAN data
            0: the first non-NAN data
    :return:
    '''
    #if to_freq=='M':
    #    data.index = data.index.where(data.index == ((data.index + MonthEnd()) - MonthEnd()),
    #                                  data.index + MonthEnd())
    #elif to_freq=='W':
    #    # By=lambda x:x.year*100+x.week # 此种方法转化为周末日期时会出现错误
    #    week_day=5 #0-6分别对应周一至周日
    #    data.index=data.index.where(data.index==((data.index+Week(weekday=week_day))-Week()),data.index+Week(weekday=week_day))
    #elif to_freq=='Y':
    #    data.index = data.index.where(data.index == ((data.index + YearEnd()) - YearEnd()),
    #                                  data.index + YearEnd())
    #By = lambda x: x
    data=_resample_index(data,to_freq)
    tmp=data.groupby(level=0)
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
def import_data(PV_vars=None, BS_vars=None,Rf_freq=None,filter_data=True):
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
        if filter_data is True:
            filter_ =pd.read_pickle(data_path + 'filtered_data').stack()
            PV=PV[filter_]
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
def cal_index_ret(freq='M',index_code='000016.SH',del_Rf=True):
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
def cal_SMB_HML(freq='M',ret=None,size=None,BM=None,percentile1=None,percentile2=None,independent=True):
    if ret is None:
        ret = cal_ret(freq=freq)
        size=cal_size(freq=freq)
        BM=cal_BM(size=size,freq=freq)
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
        delta_date=DateOffset(months=periods)
    elif freq=='W':
        delta_date = DateOffset(weeks=periods)
    for edt in EndDate_list[periods - 1:]:
        tmp1 = index_ret_d.loc[edt - delta_date + Day():edt] \
               - index_ret_d.loc[edt - delta_date + Day():edt].mean()
        tmp2 = ret_d.loc[edt - delta_date + Day():edt] \
               - ret_d.loc[edt - delta_date + Day():edt].mean()
        beta.loc[edt] = tmp2.mul(tmp1, axis=0).mean() / tmp1.var()
    #if freq=='D'or freq=='D-W':
    #    return beta[beta!=0.0].resample('D',fill_method='pad').loc[ret_d.loc[EndDate_list[periods-1]:].index].astype(float)
    #else:
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
    if periods <2:
        raise ValueError('Periods must not be lower than 2 !!!')
    if opn is None:
        p0 = 'adj_open'
        p1 = 'adj_close'
        pv = import_data(PV_vars=[p0, p1])[0]
        opn = pv[p0].unstack()
        cls = pv[p1].unstack()
    open_price=_resample_h2l(opn,to_freq=freq,n_th=0,by_position=False).shift(periods-1)
    close_price=_resample_h2l(cls,to_freq=freq,n_th=-1,by_position=False).shift(1)
    return ((close_price-open_price)/open_price*100.0).iloc[periods-1:]
def cal_rev(periods=1,freq='M',opn=None,cls=None):
    if opn is None:
        p0 = 'adj_open'
        p1 = 'adj_close'
        pv=import_data(PV_vars=[p0,p1])[0]
        opn=pv[p0].unstack()
        cls=pv[p1].unstack()
    open_price=_resample_h2l(opn,to_freq=freq,n_th=0,by_position=False).shift(periods-1)
    close_price=_resample_h2l(cls,to_freq=freq,n_th=-1,by_position=False)
    return ((close_price-open_price)/open_price*100.0).iloc[periods-1:]
def cal_illiq(periods=12,freq='M',ret_d=None,amount=None):
    if ret_d is None:
        p0='adj_open';p1='adj_close'
        pv=import_data(PV_vars=[p0,p1,'amount'])[0]
        amount = pv['amount'].unstack()
        opn=pv[p0].unstack()
        ret_d=(pv[p1].unstack()-opn)/opn*100.0
    illiq_d=1e3*np.abs(ret_d) / amount
    illiq_d=_resample_index(illiq_d,freq)
    illiq_sum = illiq_d.groupby(level=0).sum()
    illiq_count = illiq_d.groupby(level=0).count()
    illiq_sum=illiq_sum.rolling(periods).sum()
    illiq_count=illiq_count.rolling(periods).sum()
    illiq=illiq_sum/illiq_count
    return illiq.iloc[periods-1:]
def cal_turnover(periods,freq='M',amount=None,size_free=None):
    if amount is None:
        pv=import_data(PV_vars=['amount','size_free'])[0]
        amount_=pv['amount'].unstack()
        size_free=pv['size_free'].unstack()
    turnover=amount_/size_free*10.0
    #turnover=turnover[turnover>0.0001]
    turnover=_resample_index(turnover,freq)
    turnover_sum=turnover.groupby(level=0).sum()
    turnover_count=turnover.groupby(level=0).count();turnover_count[turnover_count==0.0]=np.nan
    turnover=turnover_sum.rolling(periods).sum()/turnover_count.rolling(periods).sum()
    return turnover.iloc[periods-1:]
def cal_MaxRet(periods=1,freq='M',ret_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D',del_Rf=False)
    ret_d=_resample_index(ret_d,freq)
    max_ret = ret_d.groupby(level=0).max().rolling(periods).max()
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
def cal_iskew(periods,freq='M',method='FF',ret_d=None,index_ret_d=None,SMB_d=None,HML_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D')
        index_ret_d=cal_index_ret(freq='D')
        if method=='FF':
            size=cal_size(freq='D')
            BM=cal_BM(size=size,freq='D')
            SMB_d,HML_d=cal_SMB_HML(freq='D',ret=ret_d,size=size,BM=BM)
    EndDate_list=_GetEndDateList(ret_d,freq)
    iskew = pd.DataFrame(index=EndDate_list, columns=ret_d.columns)
    if freq == 'M':
        delta_date = DateOffset(months=periods)
    elif freq == 'W':
        delta_date = DateOffset(weeks=periods)
    if method == 'CAPM':
        for edt in EndDate_list[periods - 1:]:
            tmp1 = index_ret_d.loc[edt - delta_date + Day():edt] \
                   - index_ret_d.loc[edt - delta_date + Day():edt].mean()
            tmp2 = ret_d.loc[edt - delta_date + Day():edt] \
                   - ret_d.loc[edt - delta_date + Day():edt].mean()
            iskew.loc[edt] = (tmp2 - tmp1[:, None] @ tmp2.mul(tmp1, axis=0).mean()[None, :] / tmp1.var()).skew()
    elif method=='FF':
        XY=pd.concat((pd.DataFrame({'index':index_ret_d,'SMB':SMB_d,'HML':HML_d}),ret_d),axis=1)
        for edt in EndDate_list[periods - 1:]:
            tmp1 = XY.loc[edt - delta_date + Day():edt] \
                   - XY.loc[edt - delta_date + Day():edt].mean()
            iskew.loc[edt] = (tmp1[ret_d.columns] - tmp1[['index', 'SMB', 'HML']] @ np.linalg.pinv(
                tmp1[['index', 'SMB', 'HML']].cov()) @ np.array([
                tmp1[ret_d.columns].mul(tmp1['index'], axis=0).mean(),
                tmp1[ret_d.columns].mul(tmp1['SMB'], axis=0).mean(),
                tmp1[ret_d.columns].mul(tmp1['HML'], axis=0).mean()
            ])).skew()
    return iskew.iloc[periods-1:].astype(float)
def cal_vol(periods,freq='M',ZeroMean=False,ret_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D')
    EndDate_list=_GetEndDateList(ret_d,freq)
    if freq=='M':
        if ZeroMean:
            ret_d=ret_d**2.0
            return pd.DataFrame(
                (ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean() for edt in
                 EndDate_list)
                , index=EndDate_list).iloc[periods - 1].sqrt()
        else:
            return pd.DataFrame((ret_d.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].std() for edt in EndDate_list)
                 , index=EndDate_list).iloc[periods-1]
    elif freq=='W':
        if ZeroMean:
            ret_d=ret_d**2.0
            return pd.DataFrame(
                (ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].mean() for edt in
                 EndDate_list)
                , index=EndDate_list).iloc[periods - 1].sqrt()
        else:
            return pd.DataFrame((ret_d.loc[edt - pd.tseries.offsets.DateOffset(weeks=periods) + Day():edt].std() for edt in EndDate_list)
                 , index=EndDate_list).iloc[periods-1]
def cal_ivol(periods,freq='M',method='FF',ret_d=None,index_ret_d=None,SMB_d=None,HML_d=None):
    if ret_d is None:
        if ret_d is None:
            ret_d = cal_ret(freq='D')
            index_ret_d = cal_index_ret(freq='D')
            if method == 'FF':
                size = cal_size(freq='D')
                BM = cal_BM(size=size, freq='D')
                SMB_d, HML_d = cal_SMB_HML(freq='D', ret=ret_d, size=size, BM=BM)
    EndDate_list=_GetEndDateList(ret_d,freq)
    ivol = pd.DataFrame(index=EndDate_list, columns=ret_d.columns)
    if freq == 'M':
        delta_date = DateOffset(months=periods)
    elif freq == 'W':
        delta_date = DateOffset(weeks=periods)
    if method == 'CAPM':
        for edt in EndDate_list[periods - 1:]:
            tmp1 = index_ret_d.loc[edt - delta_date + Day():edt] \
                   - index_ret_d.loc[edt - delta_date + Day():edt].mean()
            tmp2 = ret_d.loc[edt - delta_date + Day():edt] \
                   - ret_d.loc[edt - delta_date + Day():edt].mean()
            ivol.loc[edt] = (tmp2 - tmp1[:, None] @ tmp2.mul(tmp1, axis=0).mean()[None, :] / tmp1.var()).std()
    elif method=='FF':
        XY=pd.concat((pd.DataFrame({'index':index_ret_d,'SMB':SMB_d,'HML':HML_d}),ret_d),axis=1)
        for edt in EndDate_list[periods - 1:]:
            tmp1 = XY.loc[edt - delta_date + Day():edt] \
                   - XY.loc[edt - delta_date + Day():edt].mean()
            # tmp2 = ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt] \
            #       - ret.loc[edt - pd.tseries.offsets.DateOffset(months=periods) + Day():edt].mean()
            ivol.loc[edt] = (tmp1[ret_d.columns] - tmp1[['index', 'SMB', 'HML']] @ np.linalg.pinv(
                tmp1[['index', 'SMB', 'HML']].cov()) @ np.array([
                tmp1[ret_d.columns].mul(tmp1['index'], axis=0).mean(),
                tmp1[ret_d.columns].mul(tmp1['SMB'], axis=0).mean(),
                tmp1[ret_d.columns].mul(tmp1['HML'], axis=0).mean()
            ])).std()
    return ivol.iloc[periods-1:].astype(float)

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
        df1['ret_w'] = (ret*weights.shift(1)).iloc[1:][valid_].stack()
        df1['ref'] = mark_[valid_].stack()
        df2['ret_w'] = weightsweights.shift(1).iloc[1:][valid_].stack()
        df2['ref'] = mark_[valid_].stack()
        tmp1=df1.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp2 = df2.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp=(tmp1/tmp2).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    tmp.columns.name=None
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
    tmp.columns.name = None
    return tmp
def describe(df,stats=['skew','kurt']):
    d=df.describe(percentiles=[0.05,0.25,0.5,0.75,0.95])
    return d.append(df.reindex(d.columns, axis=1).agg(stats))
def NWest_mean(dataframe,L=None):
    '''
    df不要以Nan开头，会引起误差;
    :param df:
    :param L:
    :return:
    '''
    df=dataframe-dataframe.mean()
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
    if index_ret is None:
        index_ret=cal_index_ret(freq=freq)
        size=cal_size(freq=freq)
        BM=cal_BM(size,freq=freq)
        SMB,HML=cal_SMB_HML(freq=freq,size=size,BM=BM)
    independent_vars=['index', 'SMB', 'HML']#[independent_vars+['ret']]
    XY=pd.DataFrame({'index':index_ret,'SMB':SMB, 'HML':HML,'ret':arr}).dropna()
    tmp1=XY-XY.mean(axis=0)
    FF_alpha=XY['ret']-XY[independent_vars]@np.linalg.pinv(tmp1[independent_vars].cov())@tmp1[independent_vars].mul(tmp1['ret'],axis=0).mean().values
    return FF_alpha.mean(),NWest(FF_alpha-FF_alpha.mean(),XY[independent_vars])
def CumRet(port_ret,method='simple',labels=None):
    m,*n=port_ret.shape
    if labels is not None:
        port_ret.columns=labels
    cum_ret=(port_ret*0.01+1.0).cumprod()
    return cum_ret
def backtest():
    # TODO
    pass
def Fama_MacBeth(var_list,freq='M'):
    # ret_freq/ret: rev
    # ret_d: beta,illiq,max_ret,skew,coskew,iskew,vol,ivol
    # index_ret_freq/index_ret:
    # index_ret_d: beta,coskew,iskew,ivol
    # SMB_d,HML_d: iskew,ivol
    # amount_d: illiq,turnover
    # size_free_d: turnover
    #
    if set(var_list).intersection(['beta','illiq','max_ret','skew','coskew','iskew','vol','ivol']):
        ret_d=cal_ret(freq='D')
    if set(var_list).intersection(['beta','coskew','iskew','ivol']):
        index_ret_d=cal_index_ret(freq='D')
    if set(var_list).intersection(['iskew','ivol']):
        SMB_d,HML_d=cal_SMB_HML(freq='D')
    if set(var_list).intersection(['illiq','turnover']):
        amount_d=import_data(PV_vars=['amount'])[0].unstack()
    if set(var_list).intersection(['turnover']):
        size_free_d=import_data(PV_vars=['size_free'])[0].unstack()
    ### beta,mom
    freq_periods1={
        'M':12,
        'W':52
    }
    ### illiq,turnover,skew,coskew,iskew,vol,ivol
    freq_periods2={
        'M':1,
        'W':4
    }
    var_dict={}
    ret = cal_ret(freq=freq)
    var_dict['ret']=ret.stack()
    if 'beta' in var_list:
        var_dict['beta']=cal_beta(freq_periods1[freq],freq=freq,index_ret_d=index_ret_d,ret_d=ret_d).shift(1).stack()
    if 'size' in var_list or 'BM' in var_list:
        size=cal_size(freq=freq)
        var_dict['BM']=cal_BM(size=size,freq=freq).shift(1).stack()
        var_dict['size']=size.shift(1).stack()*1e-5
    if 'mom' in var_list:
        var_dict['mom']=cal_mom(freq_periods1[freq],freq=freq).shift(1).stack()
    if 'rev' in var_list:
        var_dict['rev']=ret.shift(1).stack()
    if 'illiq' in var_list:
        var_dict['illiq'] = cal_illiq(freq_periods2[freq], freq=freq, ret_d=ret_d, amount=amount_d).shift(1).stack()
    if 'max_ret' in var_list:
        var_dict['max_ret']=cal_MaxRet(freq_periods2,freq=freq,ret_d=ret_d).shift(1).stack()
    if 'turnover' in var_list:
        var_dict['turnover'] = cal_turnover(freq_periods2[freq], freq=freq, amount=amount_d,
                                            size_free=size_free_d).shift(1).stack()
    if 'skew' in var_list:
        var_dict['skew'] = cal_skew(freq_periods2[freq], freq=freq, ret_d=ret_d).shift(1).stack()
    if 'coskew' in var_list:
        var_dict['coskew'] = cal_coskew(freq_periods2[freq], freq=freq, ret_d=ret_d, index_ret_d=index_ret_d).shift(
            1).stack()
    if 'iskew' in var_list:
        var_dict['iskew'] = cal_iskew(freq_periods2[freq], freq=freq, ret_d=ret_d, index_ret_d=index_ret_d,
                                      SMB_d=SMB_d, HML_d=HML_d).shift(1).stack()

    if 'vol' in var_list:
        var_dict['vol'] = cal_vol(freq_periods2[freq], freq=freq, ret_d=ret_d).shift(1).stack()
    if 'ivol' in var_list:
        var_dict['ivol'] = cal_ivol(freq_periods2[freq], freq=freq, ret_d=ret_d, index_ret_d=index_ret_d,
                                    SMB_d=SMB_d, HML_d=HML_d).shift(1).stack()

    #global var_dict
    XY=pd.DataFrame(var_dict)[['ret']+var_list].dropna()
    date_list=sorted(set(XY.index.get_level_values(0)))
    params_list=['const']+var_list
    fit_results=pd.DataFrame(index=date_list,columns=params_list+['adj.R2'])
    for t in date_list:
        model_fit= sm.OLS(XY.loc[t]['ret'], sm.add_constant(XY.loc[t][var_list])).fit()
        fit_results.loc[t,params_list]=model_fit.params
        fit_results.loc[t,'adj.R2']=model_fit.rsquared_adj
    return fit_results.iloc[:-1].mean(),fit_results.iloc[:-1][params_list].mean() / NWest_mean(fit_results.iloc[:-1][params_list])

from time import time
t0=time()
var_list=['ivol','beta','size','BM','mom','rev','']
a,b=Fama_MacBeth(var_list)
print(time()-t0)
ivol=cal_ivol(1, freq=freq, ret_d=ret_d, index_ret_d=index_ret_d,
                                        SMB_d=SMB_d, HML_d=HML_d)
mom=cal_mom(12,freq='M')
tmp=cal_mimick_port1(mom,freq='M',ret=ret.loc[mom.index])
tmp[len(tmp.columns)+1.0]=tmp.iloc[:,-1]-tmp.iloc[:,0]
tmp.mean()/NWest_mean(tmp)
cal_FF_alpha(tmp)