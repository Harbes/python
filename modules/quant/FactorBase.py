import pandas as pd
import numpy as np
from pandas.tseries.offsets import MonthEnd,YearEnd,Week,Day,DateOffset
import statsmodels.api as sm
#from dateutil.parser import parse
import sys
import warnings
warnings.filterwarnings("ignore")
#data_path=_data_path()
def GetDataPath():
    sys_platform=sys.platform
    if sys_platform =='win32':
        return 'E:/data/NewData/'
    elif sys_platform=='mac':
        return '/Users/harbes/data/NewData/'
    elif sys_platform=='linux':
        return '/home/harbes/data/NewData/'
    else:
        raise ValueError('These is no such systerm in your work-station')
def resample_index(data, to_freq):
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
def GetEndDateList(data, freq, trim_end=False):
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
    if trim_end:
        return sorted(set(date_list))[:-1]
    else:
        return sorted(set(date_list))
def resample_h2l(data, to_freq='M', n_th=0, by_position=True):
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
    data=resample_index(data, to_freq)
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
def GetDeltaDate(periods,freq):
    if freq == 'M':
        delta_date = periods*MonthEnd()#DateOffset(months=periods)
    elif freq == 'W':
        delta_date = DateOffset(weeks=periods)
    return delta_date
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
    data_path=GetDataPath()
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
def cal_index_ret(freq='M',index_code='000016.SH',del_Rf=True,trim_end=True):
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
    data_path=GetDataPath()
    index_ret = pd.read_pickle(data_path + 'index_ret')[['index_code', 'trddt', 'opnprc', 'clsprc']]
    index_ret = index_ret[index_ret['index_code'] == index_code][['trddt', 'opnprc', 'clsprc']].set_index('trddt')
    index_ret.index=pd.to_datetime(index_ret.index.astype(int).astype(str),format='%Y%m%d')
    price0=resample_h2l(index_ret['opnprc'], to_freq=freq, n_th=0)
    price1=resample_h2l(index_ret['clsprc'], to_freq=freq, n_th=-1)
    index_ret=(price1-price0)/price0*100.0
    index_ret.index.name = 'trddt'
    if del_Rf:
        Rf=import_data(Rf_freq=freq)[2]
        #index_ret = index_ret.sub(Rf.loc[index_ret.index], axis=0) # 报错，Rf是一个DataFrame，而不是Series
        #index_ret=-Rf.loc[index_ret.index].sub(index_ret,axis=0) # 虽然成功，但是返回的是一个DataFrame，而不是Series
        #index_ret = index_ret.sub(Rf.loc[index_ret.index].iloc[:,0], axis=0)
        index_ret=(index_ret + 100.0) / (Rf.loc[index_ret.index].iloc[:, 0] * 0.01 + 1.0) - 100.0
    if trim_end:
        return index_ret.iloc[:-1].astype(float)
    else:
        return index_ret.astype(float)
def cal_ret(freq='M',periods=1,del_Rf=True,trim_end=True):
    '''

    calculate the past returns,including the current period
    :param periods:
    :param freq:
    :param del_Rf:
    :return:
    '''
    p0 = 'adj_open';p1 = 'adj_close' # price0也可以使用 'adj_pre_close'
    PV=import_data(PV_vars=[p0,p1])[0]
    price0=resample_h2l(PV[p0].unstack(), to_freq=freq, n_th=0).shift(periods - 1)
    price1=resample_h2l(PV[p1].unstack(), to_freq=freq, n_th=-1)
    ret=(price1-price0)/price0*100.0
    if del_Rf:
        Rf = import_data(Rf_freq=freq)[2]
        ret = (ret + 100.0).div((Rf.loc[ret.index].iloc[:, 0] * 0.01 + 1.0) ** periods, axis=0) - 100.0
    if trim_end:
        return ret.iloc[periods - 1:-1].astype(float)
    else:
        return ret.iloc[periods-1:].astype(float)
def cal_SMB_HML(freq='M',ret=None,size=None,BM=None,percentile1=None,percentile2=None,independent=True):
    if ret is None:
        ret = cal_ret(freq=freq)
        size=cal_size(freq=freq)
        BM=cal_BM(size=size,freq=freq)
    ret=ret.loc[size.index[0]:]
    size=size.loc[:ret.index[-1]]
    BM=BM.loc[:ret.index[-1]]
    valid_ = ~(pd.isnull(size.shift(1)) | pd.isnull(BM.shift(1)) | pd.isnull(ret))  # TypeError: bad operand type for unary ~: 'float'--->index或columns不匹配
    size = size[valid_.shift(-1).fillna(False)]
    BM = BM[valid_.shift(-1).fillna(False)]
    ret = ret[valid_]
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
    #valid_ = ~(pd.isnull(mark_1 + mark_2) | pd.isnull(ret.iloc[1:]))  # valid的股票要满足：当期有前一个月的indicator信息；当期保证交易
    df = pd.DataFrame()
    df['rtn'] = ret.iloc[1:].stack()
    df['ref1'] = mark_1.stack()
    df['ref2'] = mark_2.stack()
    tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).mean()).unstack()
    tmp.columns = tmp.columns.get_level_values(1)
    tmp.index.names = ('trddt', 'ref1')
    HML = tmp.mean(axis=0, level=0)
    SMB = tmp.mean(axis=1).unstack()
    return SMB.iloc[:,-1]-SMB.iloc[:,0], HML.iloc[:,-1]-HML.iloc[:,0]
def cal_beta(periods,freq='M',index_ret_d=None,ret_d=None):
    if index_ret_d is None:
        index_ret_d=cal_index_ret(freq='D')
        ret_d=cal_ret(freq='D')
    EndDate_list=GetEndDateList(ret_d, freq,trim_end=True)
    beta=pd.DataFrame(index=EndDate_list,columns=ret_d.columns)
    delta_date = GetDeltaDate(periods,freq)
    for edt in EndDate_list[periods - 1:]:
        tmp1 = index_ret_d.loc[edt - delta_date + Day():edt] \
               - index_ret_d.loc[edt - delta_date + Day():edt].mean()
        tmp2 = ret_d.loc[tmp1.index] \
               - ret_d.loc[tmp1.index].mean()
        beta.loc[edt] = tmp2.mul(tmp1, axis=0).mean() / tmp1.var()
    #if freq=='D'or freq=='D-W':
    #    return beta[beta!=0.0].resample('D',fill_method='pad').loc[ret_d.loc[EndDate_list[periods-1]:].index].astype(float)
    #else:
    return beta[beta!=0.0].iloc[periods-1:].astype(float)
def cal_size(freq='M'):
    size=import_data(PV_vars=['size_tot'])[0].unstack()
    size.columns=size.columns.get_level_values(1)
    size=resample_h2l(size, to_freq=freq, n_th=-1, by_position=False)
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
    open_price=resample_h2l(opn, to_freq=freq, n_th=0, by_position=False).shift(periods - 1)
    close_price=resample_h2l(cls, to_freq=freq, n_th=-1, by_position=False).shift(1)
    return ((close_price-open_price)/open_price*100.0).iloc[periods-1:]
def cal_rev(periods=1,freq='M',opn=None,cls=None):
    if opn is None:
        p0 = 'adj_open'
        p1 = 'adj_close'
        pv=import_data(PV_vars=[p0,p1])[0]
        opn=pv[p0].unstack()
        cls=pv[p1].unstack()
    open_price=resample_h2l(opn, to_freq=freq, n_th=0, by_position=False).shift(periods - 1)
    close_price=resample_h2l(cls, to_freq=freq, n_th=-1, by_position=False)
    return ((close_price-open_price)/open_price*100.0).iloc[periods-1:]
def cal_illiq(periods=12,freq='M',ret_d=None,amount=None):
    if ret_d is None:
        p0='adj_open';p1='adj_close'
        pv=import_data(PV_vars=[p0,p1,'amount'])[0]
        amount = pv['amount'].unstack()
        opn=pv[p0].unstack()
        ret_d=(pv[p1].unstack()-opn)/opn*100.0
    illiq_d=1e3*np.abs(ret_d) / amount
    illiq_d=resample_index(illiq_d, freq)
    illiq_sum = illiq_d.groupby(level=0).sum()
    illiq_count = illiq_d.groupby(level=0).count()
    illiq_sum=illiq_sum.rolling(periods).sum()
    illiq_count=illiq_count.rolling(periods).sum()
    illiq=illiq_sum/illiq_count
    return illiq.iloc[periods-1:]
def cal_turnover(periods,freq='M',amount=None,size_free=None):
    if amount is None:
        pv=import_data(PV_vars=['amount','size_free'])[0]
        amount=pv['amount'].unstack()
        size_free=pv['size_free'].unstack()
    turnover=amount/size_free*10.0
    #turnover=turnover[turnover>0.0001]
    turnover=resample_index(turnover, freq)
    turnover_sum=turnover.groupby(level=0).sum()
    turnover_count=turnover.groupby(level=0).count();turnover_count[turnover_count==0.0]=np.nan
    turnover=turnover_sum.rolling(periods).sum()/turnover_count.rolling(periods).sum()
    return turnover.iloc[periods-1:]
def cal_MaxRet(periods=1,freq='M',ret_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D',del_Rf=False)
    ret_d=resample_index(ret_d, freq)
    max_ret = ret_d.groupby(level=0).max().rolling(periods).max()
    return max_ret.iloc[periods-1:]
def cal_skew(periods,freq='M',ret_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D')
    EndDate_list=GetEndDateList(ret_d, freq,trim_end=True)
    delta_date = GetDeltaDate(periods,freq)
    return pd.DataFrame(
        (ret_d.loc[edt - delta_date + Day():edt].skew() for edt in EndDate_list)
        , index=EndDate_list).iloc[periods - 1:]
def cal_coskew(periods,freq='M',ret_d=None,index_ret_d=None):
    if index_ret_d is None:
        ret_d=cal_ret(freq='D')
        index_ret_d=cal_index_ret(freq='D')
    X=pd.DataFrame({'index':index_ret_d,'index^2':index_ret_d**2.0})[['index','index^2']]
    EndDate_list=GetEndDateList(ret_d, freq,trim_end=True)
    coskew=pd.DataFrame(index=EndDate_list,columns=ret_d.columns)
    delta_date = GetDeltaDate(periods, freq)
    for edt in EndDate_list[periods - 1:]:
        tmp1 = X.loc[edt - delta_date + Day():edt] \
               - X.loc[edt - delta_date + Day():edt].mean()
        tmp2 = ret_d.loc[tmp1.index] \
               - ret_d.loc[tmp1.index].mean()
        tmp3 = np.linalg.pinv(tmp1.cov())[1]
        coskew.loc[edt] = tmp3[0] * tmp2.mul(tmp1['index'], axis=0).mean() \
                          + tmp3[1] * tmp2.mul(tmp1['index^2'],axis=0).mean()
    return coskew.iloc[periods-1:].astype(float)
def cal_iskew(periods,freq='M',method='FF',ret_d=None,index_ret_d=None,SMB_d=None,HML_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D')
        index_ret_d=cal_index_ret(freq='D')
        if method=='FF':
            size_d=cal_size(freq='D')
            BM_d=cal_BM(size=size_d,freq='D')
            SMB_d,HML_d=cal_SMB_HML(freq='D',ret=ret_d,size=size_d,BM=BM_d)
    EndDate_list=GetEndDateList(ret_d, freq,trim_end=True)
    iskew = pd.DataFrame(index=EndDate_list, columns=ret_d.columns)
    delta_date = GetDeltaDate(periods,freq)
    if method == 'CAPM':
        for edt in EndDate_list[periods - 1:]:
            tmp1 = index_ret_d.loc[edt - delta_date + Day():edt] \
                   - index_ret_d.loc[edt - delta_date + Day():edt].mean()
            tmp2 = ret_d.loc[tmp1.index] \
                   - ret_d.loc[tmp1.index].mean()
            iskew.loc[edt] = (tmp2 - tmp1[:, None] @ tmp2.mul(tmp1, axis=0).mean()[None, :] / tmp1.var()).skew()
    elif method=='FF':
        X = pd.DataFrame({'index': index_ret_d, 'SMB': SMB_d, 'HML': HML_d})[['index', 'SMB', 'HML']].dropna()
        for edt in EndDate_list[periods - 1:]:
            tmp1 = X.loc[edt - delta_date + Day():edt] \
                   - X.loc[edt - delta_date + Day():edt].mean()
            tmp2 = ret_d.loc[tmp1.index] \
                   - ret_d.loc[tmp1.index].mean()
            iskew.loc[edt] = (tmp2 - tmp1 @ np.linalg.pinv(
                tmp1.cov()) @ np.array([
                tmp2.mul(tmp1['index'], axis=0).mean(),
                tmp2.mul(tmp1['SMB'], axis=0).mean(),
                tmp2.mul(tmp1['HML'], axis=0).mean()
            ])).skew()
    return iskew.iloc[periods-1:].astype(float)
def cal_vol(periods,freq='M',ZeroMean=False,ret_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D')
    EndDate_list=GetEndDateList(ret_d, freq,trim_end=True)
    delta_date = GetDeltaDate(periods, freq)
    if ZeroMean:
        ret_d = ret_d ** 2.0
        return pd.DataFrame(
            (ret_d.loc[edt - delta_date + Day():edt].mean() for edt in
             EndDate_list)
            , index=EndDate_list).iloc[periods - 1:].sqrt()
    else:
        return pd.DataFrame(
            (ret_d.loc[edt - delta_date + Day():edt].std() for edt in EndDate_list)
            , index=EndDate_list).iloc[periods - 1:]
def cal_ivol(periods,freq='M',method='FF',ret_d=None,index_ret_d=None,SMB_d=None,HML_d=None):
    if ret_d is None:
        ret_d = cal_ret(freq='D')
        index_ret_d = cal_index_ret(freq='D')
        if method == 'FF':
            size_d = cal_size(freq='D')
            BM_d = cal_BM(size=size_d, freq='D')
            SMB_d, HML_d = cal_SMB_HML(freq='D', ret=ret_d, size=size_d, BM=BM_d)
    EndDate_list=GetEndDateList(ret_d, freq,trim_end=True)
    ivol = pd.DataFrame(index=EndDate_list, columns=ret_d.columns)
    delta_date = GetDeltaDate(periods, freq)
    if method == 'CAPM':
        for edt in EndDate_list[periods - 1:]:
            tmp1 = index_ret_d.loc[edt - delta_date + Day():edt] \
                   - index_ret_d.loc[edt - delta_date + Day():edt].mean()
            tmp2 = ret_d.loc[tmp1.index] \
                   - ret_d.loc[tmp1.index].mean()
            ivol.loc[edt] = (tmp2 - tmp1[:, None] @ tmp2.mul(tmp1, axis=0).mean()[None, :] / tmp1.var()).std()
    elif method=='FF':
        # TODO 自从trim_end之后就存在问题，而且concat消耗资源
        X=pd.DataFrame({'index':index_ret_d,'SMB':SMB_d,'HML':HML_d})[['index', 'SMB', 'HML']].dropna()
        for edt in EndDate_list[periods - 1:]:
            tmp1 = X.loc[edt - delta_date + Day():edt] \
                   - X.loc[edt - delta_date + Day():edt].mean()
            tmp2 = ret_d.loc[tmp1.index] \
                   - ret_d.loc[tmp1.index].mean()
            ivol.loc[edt] = (tmp2 - tmp1 @ np.linalg.pinv(
                tmp1.cov()) @ np.array([
                tmp2.mul(tmp1['index'], axis=0).mean(),
                tmp2.mul(tmp1['SMB'], axis=0).mean(),
                tmp2.mul(tmp1['HML'], axis=0).mean()
            ])).std()
    return ivol.iloc[periods-1:].astype(float)

def cal_mimick_port1(indi,freq='M',ret=None,weights=None,percentile=None):
    if percentile is None:
        percentile=np.arange(0.0,1.01,0.1)
    label_=[i for i in range(1,len(percentile))]
    if ret is None:
        ret=cal_ret(freq=freq).loc[indi.index[0]:]
    else:
        ret=ret.loc[indi.index[0]:]
    indi=indi.loc[:ret.index[-1]]
    valid_=~(pd.isnull(indi.shift(1)) | pd.isnull(ret)) # TypeError: bad operand type for unary ~: 'float'--->index或columns不匹配
    indi=indi[valid_.shift(-1).fillna(False)]
    ret=ret[valid_]
    mark_ = pd.DataFrame([pd.qcut(indi.iloc[i],q=percentile, labels=label_) for i in range(len(indi)-1)],index=indi.index[1:])
    if weights is None:
        df = pd.DataFrame()
        df['rtn'] = ret.iloc[1:].stack() # ValueError: Must pass DataFrame with boolean values only--->index或columns不匹配
        df['ref'] = mark_.stack()
        tmp=df.groupby(level=0).apply(lambda g: g.groupby('ref').mean()).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    else:
        df1=pd.DataFrame()
        df2=pd.DataFrame()
        df1['ret_w'] = (ret*weights.shift(1)).iloc[1:].stack()
        df1['ref'] = mark_.stack()
        df2['ret_w'] = weights.shift(1).iloc[1:].stack()
        df2['ref'] = mark_.stack()
        tmp1=df1.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp2 = df2.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp=(tmp1/tmp2).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    tmp.columns.name=None
    return tmp
def cal_mimick_port2(indi1,indi2,freq='M',ret=None,weights=None,percentile1=None,percentile2=None,independent=False):
    if percentile1 is None:
        percentile1=np.arange(0.0,1.01,0.2)
        percentile2=np.arange(0.0,1.01,0.2)
    label_1 = [i for i in range(1, len(percentile1))]
    label_2 = [i for i in range(1, len(percentile2))]
    if ret is None:
        ret = cal_ret(freq=freq).loc[indi1.index[0]:]
    else:
        ret = ret.loc[indi1.index[0]:]
    start_=indi1.index[0] if indi1.index[0]>indi2.index[0] else indi2.index[0]
    indi1 = indi1.loc[start_:ret.index[-1]].copy()
    indi2 = indi2.loc[start_:ret.index[-1]].copy()
    valid_ = ~(pd.isnull((indi1+indi2).shift(1)) | pd.isnull(ret))  # TypeError: bad operand type for unary ~: 'float'--->index或columns不匹配
    indi1 = indi1[valid_.shift(-1).fillna(False)]
    indi2 = indi2[valid_.shift(-1).fillna(False)]
    ret = ret[valid_]
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
            tmp=pd.DataFrame([pd.qcut(indi2.iloc[i][mark_1.iloc[i]==l_],q=percentile2,labels=label_2) for i in range(len(indi2)-1)],index=mark_1.index)
            mark_2 = mark_2.combine_first(tmp)
    #valid_ = ~(pd.isnull(mark_1+mark_2) | pd.isnull(ret.iloc[1:]))  # valid的股票要满足：当期有前一个月的indicator信息；当期保证交易
    if weights is None:
        df = pd.DataFrame()
        df['rtn'] = ret.iloc[1:].stack()
        df['ref1'] = mark_1.stack()
        df['ref2'] = mark_2.stack()
        tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).mean()).unstack()
        tmp.columns = tmp.columns.get_level_values(1)
    else:
        df1 = pd.DataFrame()
        df2 = pd.DataFrame()
        df1['rtn_w'] = (ret * weights.shift(1)).loc[ret.index[1:]].stack()
        df1['ref1'] = mark_1.stack()
        df1['ref2'] = mark_2.stack()
        df2['rtn_w'] = weights.shift(1).loc[ret.index[1:]].stack()
        df2['ref1'] = mark_1.stack()
        df2['ref2'] = mark_2.stack()
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
        amount_d=import_data(PV_vars=['amount'])[0].unstack() # 如果只有一个变量时，unstack()的结果是columns是multiindex
        amount_d.columns=amount_d.columns.get_level_values(1)
    if set(var_list).intersection(['turnover']):
        size_free_d=import_data(PV_vars=['size_free'])[0].unstack()
        size_free_d.columns = size_free_d.columns.get_level_values(1)
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
        var_dict['max_ret']=cal_MaxRet(freq_periods2[freq],freq=freq,ret_d=ret_d).shift(1).stack()
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
def GetVarsFromList(var_list,freq='M'):
    if set(var_list).intersection(['beta','illiq','max_ret','skew','coskew','iskew','vol','ivol']):
        ret_d=cal_ret(freq='D')
    if set(var_list).intersection(['beta','coskew','iskew','ivol']):
        index_ret_d=cal_index_ret(freq='D')
    if set(var_list).intersection(['iskew','ivol']):
        SMB_d,HML_d=cal_SMB_HML(freq='D')
    if set(var_list).intersection(['illiq','turnover']):
        amount_d=import_data(PV_vars=['amount'])[0].unstack() # 如果只有一个变量时，unstack()的结果是columns是multiindex
        amount_d.columns=amount_d.columns.get_level_values(1)
    if set(var_list).intersection(['turnover']):
        size_free_d=import_data(PV_vars=['size_free'])[0].unstack()
        size_free_d.columns = size_free_d.columns.get_level_values(1)
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
    var_dict['ret']=ret
    if 'skew' in var_list:
        var_dict['skew'] = cal_skew(freq_periods2[freq], freq=freq, ret_d=ret_d)
    if 'coskew' in var_list:
        var_dict['coskew'] = cal_coskew(freq_periods2[freq], freq=freq, ret_d=ret_d, index_ret_d=index_ret_d)
    if 'iskew' in var_list:
        var_dict['iskew'] = cal_iskew(freq_periods2[freq], freq=freq, ret_d=ret_d, index_ret_d=index_ret_d,
                                      SMB_d=SMB_d, HML_d=HML_d)
    if 'vol' in var_list:
        var_dict['vol'] = cal_vol(freq_periods2[freq], freq=freq, ret_d=ret_d)
    if 'ivol' in var_list:
        var_dict['ivol'] = cal_ivol(freq_periods2[freq], freq=freq, ret_d=ret_d, index_ret_d=index_ret_d,
                                    SMB_d=SMB_d, HML_d=HML_d)
    if 'beta' in var_list:
        var_dict['beta']=cal_beta(freq_periods1[freq],freq=freq,index_ret_d=index_ret_d,ret_d=ret_d)
    if 'size' in var_list or 'BM' in var_list:
        size=cal_size(freq=freq)
        var_dict['BM']=cal_BM(size=size,freq=freq)
        var_dict['size']=size*1e-5
    if 'mom' in var_list:
        var_dict['mom']=cal_mom(freq_periods1[freq],freq=freq)
    if 'rev' in var_list:
        var_dict['rev']=ret
    if 'illiq' in var_list:
        var_dict['illiq'] = cal_illiq(freq_periods2[freq], freq=freq, ret_d=ret_d, amount=amount_d)
    if 'max_ret' in var_list:
        var_dict['max_ret']=cal_MaxRet(freq_periods2[freq],freq=freq,ret_d=ret_d)
    if 'turnover' in var_list:
        var_dict['turnover'] = cal_turnover(freq_periods2[freq], freq=freq, amount=amount_d,
                                            size_free=size_free_d)

    return var_dict
def SinglePortAnalysis(var_list,var_dict=None,freq='M',value_weighted=False):
    # TODO 待解决重复计算相同数据的非效率问题
    if var_dict is None:
        var_dict = GetVarsFromList(var_list, freq)
    index_ret = cal_index_ret(freq=freq)
    if 'size' in var_list or 'BM' in var_list:
        SMB, HML = cal_SMB_HML(freq=freq,ret=var_dict['ret'],size=var_dict['size'],BM=var_dict['BM'])
    else:
        SMB, HML = cal_SMB_HML(freq=freq)
    if value_weighted:
        if 'size' in var_list or 'BM' in var_list:
            weights=var_dict['size']
        else:
            weights=cal_size(freq=freq)
    else:
        weights=None
    percentile_=np.arange(0.0,1.01,0.1)
    portfolio_mean=pd.DataFrame(index=var_list,columns=[i for i in range(1,len(percentile_)+1)]+['alpha'])
    portfolio_t=pd.DataFrame(index=portfolio_mean.index,columns=portfolio_mean.columns)
    for var in var_list:
        tmp=cal_mimick_port1(indi=var_dict[var],freq=freq,ret=var_dict['ret'],weights=weights,percentile=percentile_)
        tmp[len(percentile_)]=tmp.iloc[:,-1]-tmp.iloc[:,0]
        portfolio_mean.loc[var]=tmp.mean()
        portfolio_t.loc[var]=portfolio_mean.loc[var]/NWest_mean(tmp)
        portfolio_mean.loc[var,'alpha'],se=cal_FF_alpha(tmp[len(percentile_)],freq=freq,index_ret=index_ret,SMB=SMB,HML=HML)
        portfolio_t.loc[var,'alpha']=portfolio_mean.loc[var,'alpha']/se
    return pd.DataFrame({'mean':portfolio_mean.stack(),'t':portfolio_t.stack()}).T
def DoublePortAnalysis(var_list,var2,var_dict=None,freq='M',value_weighted=False):
    # TODO 待解决重复计算相同数据的非效率问题
    if var_dict is None:
        var_dict = var_dict = GetVarsFromList(var_list+[var2], freq)
    index_ret = cal_index_ret(freq=freq)
    if 'size' in var_list or 'BM' in var_list:
        SMB, HML = cal_SMB_HML(freq=freq, ret=var_dict['ret'], size=var_dict['size'], BM=var_dict['BM'])
    else:
        SMB, HML = cal_SMB_HML(freq=freq)
    if value_weighted:
        if 'size' in var_list or 'BM' in var_list:
            weights = var_dict['size']
        else:
            weights = cal_size(freq=freq)
    else:
        weights=None
    percentile1 = np.arange(0.0, 1.01, 0.2)
    percentile2 = np.arange(0.0, 1.01, 0.2)
    portfolio_mean = pd.DataFrame(index=var_list, columns=range(1, len(percentile1) + 1))
    portfolio_t = pd.DataFrame(index=portfolio_mean.index, columns=portfolio_mean.columns)
    portfolio_alpha=pd.DataFrame(index=portfolio_mean.index, columns=portfolio_mean.columns)
    portfolio_alpha_t = pd.DataFrame(index=portfolio_mean.index, columns=portfolio_mean.columns)
    for var in var_list:
        tmp = cal_mimick_port2(indi1=var_dict[var],indi2=var_dict[var2], freq=freq, ret=var_dict['ret'], weights=weights,
                               percentile1=percentile1,percentile2=percentile2)
        long_short=(tmp.iloc[:,-1]-tmp.iloc[:,0]).unstack()
        long_short[len(percentile1)] = long_short.mean(axis=1)
        portfolio_mean.loc[var] = long_short.mean()
        portfolio_t.loc[var] = portfolio_mean.loc[var] / NWest_mean(long_short)
        for port in long_short.columns:
            portfolio_alpha.loc[var,port], se = cal_FF_alpha(long_short[port], freq=freq, index_ret=index_ret,
                                                            SMB=SMB, HML=HML)
            portfolio_alpha_t.loc[var, port] = portfolio_alpha.loc[var,port]/se
    return pd.DataFrame({'mean': portfolio_mean.stack(), 't': portfolio_t.stack(),'alpha':portfolio_alpha.stack(),'alpha_t':portfolio_alpha_t.stack()}).T





