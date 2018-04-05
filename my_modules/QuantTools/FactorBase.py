import pandas as pd
from pandas import DataFrame,Series,qcut
import numpy as np
from pandas.tseries.offsets import MonthEnd,YearEnd,Week,Day,DateOffset
import statsmodels.api as sm
from copy import deepcopy
from scipy.stats import mstats
from dateutil.parser import parse
import sys
import warnings
warnings.filterwarnings("ignore")
#data_path=_data_path()
#
def InputOptions():
    # TODO
    # 时间：投资时间点，调仓频率
    #      因子计算窗口期：日历窗口期or交易日窗口期
    # 收益率计算：method=arithmetic、log
    #      数据：close-close，close-open
    pass
def PreCalData():
    # TODO
    # 例如：ret_d,ret_index_d,SMB_freq,HML_freq
    pass
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
def InvestTime():
    start_=parse('20050101')
    end_=parse('20180301')
    return None


def resample_index(dat, to_freq):
    '''
    使用时一定要注意，此命令会更改数据的index；因此，凡是涉及输入的数据使用此命令时，一定要使用copy()，以防出错
    :param data:
    :param to_freq:
    :return:
    '''
    data=dat.copy()
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
    '''
    trim主要用于日度数据resample成低频数据
    :param data:
    :param freq:
    :param trim_end:
    :return:
    '''
    if freq=='M':
        date_list=data.index.where(data.index == ((data.index + MonthEnd()) - MonthEnd()),
                                      data.index + MonthEnd())
    elif freq=='W':
        week_day = 5  # 0-6分别对应周一至周日
        date_list = data.index.where(data.index == ((data.index + Week(weekday=week_day)) - Week()),
                                      data.index + Week(weekday=week_day))
    elif freq=='Y':
        date_list = data.index.where(data.index == ((data.index + YearEnd()) - YearEnd()),
                                      data.index + YearEnd())
    if trim_end:
        return sorted(set(date_list))[:-1]
    else:
        return sorted(set(date_list))
def GetValueWeightedIndexReturn(freq):
    p0 = 'adj_open'
    p1 = 'adj_close'  # price0也可以使用 'adj_pre_close'
    PV = import_data(PV_vars=['size_tot',p0, p1])[0]
    price0 = resample_index(PV[p0].unstack(),to_freq=freq).groupby(level=0).nth(0).stack()
    price1 = resample_index(PV[p1].unstack(),to_freq=freq).groupby(level=0).nth(-1).stack()
    size=resample_index(PV['size_tot'].unstack(),to_freq=freq).groupby(level=0).nth(-1).stack()
    #df=pd.concat(((price1-price0)/price0,size),axis=1).dropna()
    df=pd.DataFrame({'ret':(price1-price0)/price0,'size':size}).dropna()
    return (df['ret']*df['size']).sum(level=0)/df['size'].sum(level=0)* 100.0
def resample_h2l(data, to_freq='M', n_th=0, by_position=True):
    '''
    使用时一定要注意，此命令会更改数据的index
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
    if by_position:
        return tmp.nth(n_th)
    elif n_th:
        return tmp.last()
    else:
        return tmp.first()
def GetDeltaDate(periods,freq):
    if freq == 'M':
        return MonthEnd(periods)#DateOffset(months=periods)
    elif freq == 'W':
        return DateOffset(weeks=periods)
    elif freq=='D':
        return Day(periods)
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
        资产负债表
        'fin_year','ann_dt', 'stkcd', 'tot_assets', 'tot_liab'
        【book数据可以直接使用'tot_shrhldr_eqy_excl_min_int】
        利润表
        's_fa_eps_basic'


    :param Rf_freq:
        'Y','M','D'
        [Rf数据来源于csmar]

    :return:
    '''
    if PV_vars is None and BS_vars is None and Rf_freq is None:
        raise ValueError("You haven't name the varables to be imported")
    data_path=GetDataPath()
    if PV_vars is not None:
        if len(PV_vars)>1:
            PV=pd.read_pickle(data_path+'PV_datetime')[PV_vars]
        else:
            PV = pd.read_pickle(data_path + 'PV_datetime')[PV_vars[0]]
        if filter_data is True:
            filter_ =pd.read_pickle(data_path + 'non_ST_IPO_NT_datetime').stack() # non_ST_IPO_NT_SmallVolume_datetime
            PV=PV[filter_]
    else:
        PV=None
    if BS_vars is not None:
        BS = pd.read_pickle(data_path + 'BS')[['ann_dt', 'stkcd']+BS_vars]  #'BS'
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
def ClipQuantile(indi,percentile,label,output_indi=True):
    '''
    label的标准是：给欲保留的数据设置大于0的数
    :param indi:
    :param percentile:
    :param label:
    :return:
    '''
    if len(label)!=len(percentile)-1:
        raise ValueError('len(label)!=len(percentile)-1')
    mark_=DataFrame([qcut(indi.loc[t],q=percentile,labels=label) for t in indi.index])
    if output_indi:
        return indi[mark_>0.0]
    else:
        return mark_>0.0
def winsorize(indi,limits):
    '''
    for example:
        percentile=[0.0,0.3,0.7,1.0]
        label=[0.0,1.0,0.0] # 保留非0标记的分位区间
    :param indi:
    :param percentile:
    :param label:
    :param freq:
    :return:
    '''
    return DataFrame([mstats.winsorize(indi.loc[t],limits=limits,inplace=True) for t in indi.index],index=indi.index,columns=indi.columns)
def winsorize_vars(var_dict,limits,var_list=None):
    '''
    会修改输入数据var_dict
    :param var_dict:
    :param percentile:
    :param label:
    :param var_list:
    :return:
    '''
    var_d=deepcopy(var_dict)
    if var_list is None:
        var_l=var_d.keys()
    else:
        var_l=set(var_list).intersection(var_d.keys())
    for i in var_l:
        var_d[i]=winsorize(var_d[i],limits)
    return var_d
def cal_index_ret(freq='M',index_code=None,del_Rf=True,trim_end=True):
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
    data_path = GetDataPath()
    if index_code is None:
        try:
            index_ret=pd.read_pickle(data_path+'value_weighted_index_return_'+freq)
        except(NameError,FileNotFoundError):
            index_ret=GetValueWeightedIndexReturn(freq)
            index_ret.to_pickle(data_path+'value_weighted_index_return_'+freq)
    else:
        index_ret = pd.read_pickle(data_path + 'index_ret')[['index_code', 'trddt', 'opnprc', 'clsprc']]
        index_ret = index_ret[index_ret['index_code'] == index_code][['trddt', 'opnprc', 'clsprc']].set_index('trddt')
        index_ret.index = pd.to_datetime(index_ret.index.astype(int).astype(str), format='%Y%m%d')
        price0 = resample_h2l(index_ret['opnprc'], to_freq=freq, n_th=0)
        price1 = resample_h2l(index_ret['clsprc'], to_freq=freq, n_th=-1)
        index_ret = (price1 - price0) / price0 * 100.0
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
    return beta[beta!=0.0].shift(1).iloc[periods:].astype(float)
def cal_size(freq='M',shift_num=1,exclude=False):
    size=import_data(PV_vars=['size_tot'])[0].unstack()
    #size.columns=size.columns.get_level_values(1)
    size=resample_h2l(size, to_freq=freq, n_th=-1, by_position=False)
    if exclude:
        size=ClipQuantile(size,[0.0,0.3,1.0],[-1.0,1.0])
    return size[size>0.0].shift(shift_num).iloc[shift_num:]*1.0e-5
def cal_BM(size=None,freq='M',measure='BM'):
    if measure=='BM':
        if size is None:
            size = cal_size(freq=freq, shift_num=0)
        book=import_data(BS_vars=['tot_shrhldr_eqy_excl_min_int'])[1]
        book = book.drop(book.index[book.index.duplicated(keep='last')]).unstack()['tot_shrhldr_eqy_excl_min_int']
        BM=book.resample('D').first().ffill().loc[size.index,size.columns]/size*1.0e-8
        return BM[BM > 0.0].shift(1).iloc[1:]
    elif measure=='EP':
        # TODO 寻找相关的数据，例如
        '''
        不知道为什么，cnrds上的每股收益与wind上的eps_basic对不上？
        '''
        price, book = import_data(PV_vars=['clsprc'], BS_vars=['s_fa_eps_basic'])[:2]
        book = book.drop(book.index[book.index.duplicated(keep='last')]).unstack()['s_fa_eps_basic']
        price = resample_h2l(price.unstack(), to_freq=freq, n_th=-1)
        #BM = book.resample('D').first().ffill().loc[price.index, price.columns] / price * 100.0
        BM = book.fillna(method='ffill').loc[price.index, price.columns] / price * 100.0
        return BM.shift(1).iloc[1:]
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
    return ((close_price-open_price)/open_price*100.0).shift(1).iloc[periods:]
def cal_rev(periods=1,freq='M',opn=None,cls=None):
    if opn is None:
        p0 = 'adj_open'
        p1 = 'adj_close'
        pv=import_data(PV_vars=[p0,p1])[0]
        opn=pv[p0].unstack()
        cls=pv[p1].unstack()
    open_price=resample_h2l(opn, to_freq=freq, n_th=0, by_position=False).shift(periods - 1)
    close_price=resample_h2l(cls, to_freq=freq, n_th=-1, by_position=False)
    return ((close_price-open_price)/open_price*100.0).shift(1).iloc[periods:]
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
    return illiq.shift(1).iloc[periods:]
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
    return turnover.shift(1).iloc[periods:]
def cal_MaxRet(periods=1,freq='M',ret_d=None):
    '''
    也许15个交易日的效果更好
    :param periods:
    :param freq:
    :param ret_d:
    :return:
    '''
    if ret_d is None:
        ret_d=cal_ret(freq='D',del_Rf=False)
    ret=resample_index(ret_d, freq)
    max_ret = ret.groupby(level=0).max().rolling(periods).max()
    return max_ret.shift(1).iloc[periods:]
def cal_skew(periods,freq='M',ret_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D')
    EndDate_list=GetEndDateList(ret_d, freq,trim_end=True)
    delta_date = GetDeltaDate(periods,freq)
    return pd.DataFrame(
        (ret_d.loc[edt - delta_date + Day():edt].skew() for edt in EndDate_list)
        , index=EndDate_list).shift(1).iloc[periods:]
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
    return coskew.shift(1).iloc[periods:].astype(float)
def cal_iskew(periods,freq='M',method='FF',ret_d=None,index_ret_d=None,SMB_d=None,HML_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D')
        index_ret_d=cal_index_ret(freq='D')
        if method=='FF':
            size_d = cal_size(freq='D', shift_num=0)
            BM_d = cal_BM(size_d, freq='D')
            SMB_d, HML_d = cal_SMB_HML(ret_d, size_d.shift(1).iloc[1:], BM_d)
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
    return iskew.shift(1).iloc[periods:].astype(float)

def cal_vol(periods,freq='M',ZeroMean=False,ret_d=None):
    if ret_d is None:
        ret_d=cal_ret(freq='D')
    EndDate_list=GetEndDateList(ret_d, freq,trim_end=True)
    delta_date = GetDeltaDate(periods, freq)
    if ZeroMean:
        ret_d = ret_d ** 2.0
        vol=DataFrame(
            (ret_d.loc[edt - delta_date + Day():edt].mean() for edt in
             EndDate_list)
            , index=EndDate_list).shift(1).iloc[periods:].sqrt()
    else:
        vol=pd.DataFrame(
            (ret_d.loc[edt - delta_date + Day():edt].std() for edt in EndDate_list)
            , index=EndDate_list).shift(1).iloc[periods:]
    return vol[vol>0.0].shift(1).iloc[periods:]
def cal_ivol(periods,freq='M',method='FF',ret_d=None,index_ret_d=None,SMB_d=None,HML_d=None):
    if ret_d is None:
        ret_d = cal_ret(freq='D')
        index_ret_d = cal_index_ret(freq='D')
        if method == 'FF':
            size_d = cal_size(freq='D',shift_num=0)
            BM_d = cal_BM(size_d, freq='D')
            SMB_d, HML_d = cal_SMB_HML( ret_d,size_d.shift(1).iloc[1:],BM_d)
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
    return ivol[ivol>0.0].shift(1).iloc[periods:].astype(float)
def cal_ivol_from_EGARCH():
    return pd.read_pickle(GetDataPath()+'Expected_iVol_using_OptimalEgarch').shift(1).loc['2007-09']
#ivol[ivol>0.0]
#tmp2.loc['2015-12','603999.SH']
def IndexAlign(*vars):
    l=len(vars)
    index=vars[0].index.copy()
    for i in range(1,l):
        index &=vars[i].index
    return [vars[i].loc[index] for i in range(l)]
def cal_SMB_HML(ret,size,BM,percentile1=None,percentile2=None,independent=True,exclude_30_small_size=False):
    if exclude_30_small_size:
        size = ClipQuantile(size, [0.0, 0.3, 1.0], [-1.0, 1.0])
    ret,size,BM=IndexAlign(ret,size,BM)
    valid_ = ~pd.isnull(BM+ret+size)  # TypeError: bad operand type for unary ~: 'float'--->index或columns不匹配
    size = size[valid_]
    BM = BM[valid_]
    ret = ret[valid_]
    if percentile1 is None:
        percentile1 = [0.0, 0.5, 1.0]  # size
        percentile2 = [0.0, 0.3, 0.7, 1.0]  # value
    label_1=[i+1 for i in range(len(percentile1)-1)]
    label_2=[i+1 for i in range(len(percentile2)-1)]
    if independent:
        #mark_1 = pd.DataFrame([pd.qcut(size.iloc[i], q=percentile1, labels=label_1) for i in size.index[:-1]],
        #                      index=size.index[:-1]) # 报错
        mark_1 = DataFrame([qcut(size.loc[i], q=percentile1, labels=label_1) for i in size.index])
        mark_2 = DataFrame([qcut(BM.loc[i], q=percentile2, labels=label_2) for i in BM.index])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
    else:
        mark_1 = DataFrame([qcut(size.loc[i], q=percentile1, labels=label_1) for i in size.index])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
        mark_2=DataFrame(index=mark_1.index,columns=mark_1.columns)
        for l_ in label_1:
            tmp=DataFrame([qcut(BM.loc[i][mark_1.iloc[i]==l_],q=percentile2,labels=label_2) for i in BM.index])
            mark_2 = mark_2.combine_first(tmp)
    #valid_ = ~(pd.isnull(mark_1 + mark_2) | pd.isnull(ret.iloc[1:]))  # valid的股票要满足：当期有前一个月的indicator信息；当期保证交易
    df = DataFrame()
    df['rtn'] = ret.stack()
    df['ref1'] = mark_1.stack()
    df['ref2'] = mark_2.stack()
    tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).mean()).unstack()['rtn']
    #tmp.columns = tmp.columns.get_level_values(1)
    tmp.index.names = ('trddt', 'ref1')
    HML = tmp.mean(axis=0, level=0)
    SMB = tmp.mean(axis=1).unstack()
    return SMB.iloc[:,-1]-SMB.iloc[:,0], HML.iloc[:,-1]-HML.iloc[:,0]
def cal_mimick_port1(indi, freq='M', ret=None, weights=None, percentile=None):
    if percentile is None:
        percentile=np.arange(0.0,1.01,0.2)
    label_=[i for i in range(1,len(percentile))]
    if ret is None:
        ret= cal_ret(freq=freq).loc[indi.index[0]:]
    else:
        ret= ret.loc[indi.index[0]:].copy()
    indi= indi.loc[:ret.index[-1]].copy()
    valid_=~pd.isnull(indi+ret) # TypeError: bad operand type for unary ~: 'float'--->index或columns不匹配
    indi=indi[valid_]
    ret=ret[valid_]
    mark_ = pd.DataFrame([pd.qcut(indi.loc[i],q=percentile, labels=label_) for i in indi.index])
    if weights is None:
        df = pd.DataFrame()
        df['rtn'] = ret.stack() # ValueError: Must pass DataFrame with boolean values only--->index或columns不匹配
        df['ref'] = mark_.stack()
        tmp=df.groupby(level=0).apply(lambda g: g.groupby('ref').mean()).unstack()['rtn']
        #tmp.columns = tmp.columns.get_level_values(1)
    else:
        df1=pd.DataFrame()
        df2=pd.DataFrame()
        df1['rtn_w'] = (ret*weights.loc[ret.index]).stack()
        df1['ref'] = mark_.stack()
        df2['rtn_w'] = weights.loc[ret.index][valid_].stack()
        #df2['ref'] = mark_.stack()
        df2['ref']=df1['ref']
        tmp1=df1.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp2 = df2.groupby(level=0).apply(lambda g: g.groupby('ref').sum())
        tmp=(tmp1/tmp2).unstack()['rtn_w']
        #tmp.columns = tmp.columns.get_level_values(1)
    #tmp.columns.name=None
    return tmp
def cal_mimick_port2(indi1,indi2,freq='M',ret=None,weights=None,percentile1=None,percentile2=None,independent=False):
    if percentile1 is None:
        percentile1=np.arange(0.0,1.01,0.2)
        percentile2=np.arange(0.0,1.01,0.2)
    label_1 = [i for i in range(1, len(percentile1))]
    label_2 = [i for i in range(1, len(percentile2))]
    if ret is None:
        ret = cal_ret(freq=freq)
    indi1,indi2,ret=IndexAlign(indi1,indi2,ret)
    valid_ = ~pd.isnull(indi1+indi2+ret)  # TypeError: bad operand type for unary ~: 'float'--->index或columns不匹配
    indi1 = indi1[valid_]
    indi2 = indi2[valid_]
    ret = ret[valid_]
    if independent:
        mark_1 = pd.DataFrame([pd.qcut(indi1.loc[i], q=percentile1, labels=label_1) for i in indi1.index])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
        mark_2 = pd.DataFrame([pd.qcut(indi2.loc[i], q=percentile2, labels=label_2) for i in indi2.index])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
    else:
        mark_1 = pd.DataFrame([pd.qcut(indi1.loc[i], q=percentile1, labels=label_1) for i in indi1.index])  # indi已经shift(1)了，也就是其时间index与holding period of portfolio是一致的
        mark_2=pd.DataFrame(index=mark_1.index,columns=mark_1.columns)
        for l_ in label_1:
            tmp=pd.DataFrame([pd.qcut(indi2.loc[i][mark_1.loc[i]==l_],q=percentile2,labels=label_2) for i in indi2.index])
            mark_2 = mark_2.combine_first(tmp)
    #valid_ = ~(pd.isnull(mark_1+mark_2) | pd.isnull(ret.iloc[1:]))  # valid的股票要满足：当期有前一个月的indicator信息；当期保证交易
    if weights is None:
        df = pd.DataFrame()
        df['rtn'] = ret.stack()
        df['ref1'] = mark_1.stack()
        df['ref2'] = mark_2.stack()
        tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).mean()).unstack()['rtn']
        #tmp.columns = tmp.columns.get_level_values(1)
    else:
        df1 = pd.DataFrame()
        df2 = pd.DataFrame()
        df1['rtn_w'] = (ret * weights.loc[ret.index]).stack()
        df1['ref1'] = mark_1.stack()
        df1['ref2'] = mark_2.stack()
        df2['rtn_w'] = weights.loc[ret.index][valid_].stack()
        df2['ref1'] = mark_1.stack()
        df2['ref2'] = mark_2.stack()
        tmp1 = df1.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).sum())
        tmp2 = df2.groupby(level=0).apply(lambda g: g.groupby(['ref1','ref2']).sum())
        tmp = (tmp1 / tmp2).unstack()['rtn_w']
        #tmp.columns = tmp.columns.get_level_values(1)
    #tmp.columns.name = None
    return tmp
def describe(df,stats=('skew','kurt')):
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
    '''

    :param e:
        e是残差序列
    :param X:
        X是不包含常数项的解释变量
    :param L:
        L是滞后阶数
    :return:
    '''
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
        ret=cal_ret(freq=freq)
        index_ret=cal_index_ret(freq=freq)
        size=cal_size(freq=freq,shift_num=0)
        BM=cal_BM(size,freq=freq)
        SMB,HML=cal_SMB_HML(ret,size.shift(1).iloc[1:],BM)
        #del ret,size,BM
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
def GetVarsFromList(var_list,freq):
    if set(var_list).intersection(['beta','illiq','max_ret','skew','coskew','iskew','vol','ivol']):
        ret_d=cal_ret(freq='D')
        if set(var_list).intersection(['iskew', 'ivol']):
            size_d=cal_size(freq='D',shift_num=0)
            BM_d=cal_BM(size_d,freq='M')
            SMB_d, HML_d = cal_SMB_HML(ret_d, size_d.shift(1).iloc[1:], BM_d)
    if set(var_list).intersection(['beta','coskew','iskew','ivol']):
        index_ret_d=cal_index_ret(freq='D')
    if set(var_list).intersection(['illiq','turnover']):
        amount_d=import_data(PV_vars=['amount'])[0].unstack() # 如果只有一个变量时，unstack()的结果是columns是multiindex
        #amount_d.columns=amount_d.columns.get_level_values(1)
    if set(var_list).intersection(['turnover']):
        size_free_d=import_data(PV_vars=['size_free'])[0].unstack()
        #size_free_d.columns = size_free_d.columns.get_level_values(1)
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
    if 'beta' in var_list:
        var_dict['beta']=cal_beta(freq_periods1[freq],freq=freq,index_ret_d=index_ret_d,ret_d=ret_d)
    if 'size' in var_list or 'BM' in var_list:
        size=cal_size(freq=freq,shift_num=0)
        var_dict['BM']=cal_BM(size=size,freq=freq)
        var_dict['size']=size.shift(1).iloc[1:]
    if 'mom' in var_list:
        var_dict['mom']=cal_mom(freq_periods1[freq],freq=freq)
    if 'rev' in var_list:
        var_dict['rev']=ret.shift(1).iloc[1:]
    if 'illiq' in var_list:
        var_dict['illiq'] = cal_illiq(freq_periods2[freq], freq=freq, ret_d=ret_d, amount=amount_d)
    if 'max_ret' in var_list:
        var_dict['max_ret']=cal_MaxRet(freq_periods2[freq],freq=freq,ret_d=ret_d)
    if 'turnover' in var_list:
        var_dict['turnover'] = cal_turnover(freq_periods2[freq], freq=freq, amount=amount_d,
                                            size_free=size_free_d)
    if 'skew' in var_list:
        var_dict['skew'] = cal_skew(freq_periods2[freq], freq=freq, ret_d=ret_d)
    if 'coskew' in var_list:
        var_dict['coskew'] = cal_coskew(4, freq=freq, ret_d=ret_d, index_ret_d=index_ret_d)
    if 'iskew' in var_list:
        var_dict['iskew'] = cal_iskew(freq_periods2[freq], freq=freq, ret_d=ret_d, index_ret_d=index_ret_d,
                                      SMB_d=SMB_d, HML_d=HML_d)
    if 'vol' in var_list:
        var_dict['vol'] = cal_vol(freq_periods2[freq], freq=freq, ret_d=ret_d)
    if 'ivol' in var_list:
        var_dict['ivol'] = cal_ivol(freq_periods2[freq], freq=freq, ret_d=ret_d, index_ret_d=index_ret_d,
                                    SMB_d=SMB_d, HML_d=HML_d)

    return var_dict
def Fama_MacBeth(var_list,freq='M',var_dict=None):
    # ret_freq/ret: rev
    # ret_d: beta,illiq,max_ret,skew,coskew,iskew,vol,ivol
    # index_ret_freq/index_ret:
    # index_ret_d: beta,coskew,iskew,ivol
    # SMB_d,HML_d: iskew,ivol
    # amount_d: illiq,turnover
    # size_free_d: turnover
    #
    if var_dict is None:
        var_dict=GetVarsFromList(var_list,freq)
    var_d={}
    for k in var_list+['ret']:
        var_d[k]=var_dict[k].stack()
    if 'size' in var_d.keys():
        var_d['size']=np.log(var_d['size'])
    XY=pd.DataFrame(var_d)[['ret']+var_list].dropna()
    date_list=sorted(set(XY.index.get_level_values(0)))
    params_list=['const']+var_list
    fit_results=DataFrame(index=date_list,columns=params_list+['adj.R2'])
    for t in date_list:
        model_fit= sm.OLS(XY.loc[t]['ret'], sm.add_constant(XY.loc[t][var_list])).fit()
        fit_results.loc[t,params_list]=model_fit.params
        fit_results.loc[t,'adj.R2']=model_fit.rsquared_adj
    return pd.DataFrame({'mean':fit_results.iloc[:-1].mean(),'t':fit_results.iloc[:-1][params_list].mean() / NWest_mean(fit_results.iloc[:-1][params_list])}).T[params_list+['adj.R2']]
def Get_Index_SMB_HML(freq,var_dict):
    index_ret = cal_index_ret(freq=freq)
    if 'size' in var_dict.keys():
        SMB, HML = cal_SMB_HML(var_dict['ret'], var_dict['size'], var_dict['BM'])
    else:
        size = cal_size(freq=freq, shift_num=0)
        BM = cal_BM(size, freq='M')
        SMB, HML = cal_SMB_HML(var_dict['ret'], size.shift(1).iloc[1:], BM)
    return index_ret,SMB,HML
def SinglePortAnalysis(var_list,var_dict=None,index_ret=None,SMB=None,HML=None,freq='M',value_weighted=False):
    if var_dict is None:
        var_dict = GetVarsFromList(var_list, freq)
    if index_ret is None:
        index_ret = cal_index_ret(freq=freq)
        if 'size' in var_list or 'BM' in var_list:
            SMB, HML = cal_SMB_HML(var_dict['ret'], var_dict['size'], var_dict['BM'])
        else:
            size = cal_size(freq=freq, shift_num=0)
            BM = cal_BM(size, freq='M')
            SMB, HML = cal_SMB_HML(var_dict['ret'], size.shift(1).iloc[1:], BM)
    if value_weighted:
        if 'size' in var_list or 'BM' in var_list:
            weights=var_dict['size']
        elif index_ret is None:
            weights=size.shift(1).iloc[1:] # 需要上面的size
        else:
            weights=cal_size(freq=freq)
    else:
        weights=None
    percentile_=np.arange(0.0,1.01,0.2)
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
def DoublePortAnalysis(var_list,var2,var_dict=None,index_ret=None,SMB=None,HML=None,freq='M',value_weighted=False):
    if var_dict is None:
        var_dict = GetVarsFromList(var_list+[var2], freq)
    if index_ret is None:
        index_ret = cal_index_ret(freq=freq)
        if 'size' in var_list or 'BM' in var_list:
            SMB, HML = cal_SMB_HML(var_dict['ret'], var_dict['size'], var_dict['BM'])
        else:
            size = cal_size(freq=freq, shift_num=0)
            BM = cal_BM(size, freq='M')
            SMB, HML = cal_SMB_HML(var_dict['ret'], size.shift(1).iloc[1:], BM)
    if value_weighted:
        if 'size' in var_list or 'BM' in var_list:
            weights=var_dict['size']
        elif index_ret is None:
            weights=size.shift(1).iloc[1:]
        else:
            weights=cal_size(freq=freq)
    else:
        weights=None
    percentile1 = np.arange(0.0, 1.01, 0.2)
    percentile2 = np.arange(0.0, 1.01, 0.2)
    portfolio_mean = DataFrame(index=set(var_list).difference([var2]), columns=range(1, len(percentile1) + 1))
    portfolio_t = DataFrame(index=portfolio_mean.index, columns=portfolio_mean.columns)
    portfolio_alpha=DataFrame(index=portfolio_mean.index, columns=portfolio_mean.columns)
    portfolio_alpha_t = DataFrame(index=portfolio_mean.index, columns=portfolio_mean.columns)
    for var in portfolio_mean.index:
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


