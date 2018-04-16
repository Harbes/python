import pandas as pd
from pandas import DataFrame,Series,qcut
import numpy as np
from pandas.tseries.offsets import MonthEnd,YearEnd,Week,Day,DateOffset,MonthBegin
import statsmodels.api as sm
from copy import deepcopy
from scipy.stats import mstats
from dateutil.parser import parse
import sys
#import warnings
#warnings.filterwarnings("ignore")

# 下面三条命令添加基本设置
# options=InputOptions()
# EndDate=EndDateOfPerInvestment()
# measure_window=VarsMeasureWindow()

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
def InputOptions():
    options={
        'data_path':GetDataPath(),
        'data_start':parse('20050104'),
        'start':parse('20050101'), # 数据的起始点
        'end':parse('20180228'),

        'calendar': 'calendar',  # 或者 'trading'
        'calendar_periods':MonthEnd(),#DateOffset(years=0,months=1,weeks=0,days=0),#Week(1,weekday=5),# # #  调仓周期
        'freq': 'M', # 当'calendar':'calendar'，一定要与'calendar_periods'中的日期频率对应；当'calendar':'trading'，一定为'D'
        'num_periods': 1.0,  # 当'calendar':'calendar'，一定要与'calendar_periods'中的数字对应；当'calendar':'trading'，表征调仓周期(交易日个数)
        'num_reverse':0,#搭配MonthEnd使用，0表示倒数第一天

        'open_price':'adj_open',
        'close_price':'adj_close',
        'pv_name':'PV_datetime', # price&volume
        'bs_name':'BS',# balance of sheet
        'bs_date_name':'ann_dt', #'fin_year',#
        'filter_data_name':'non_ST_IPO_NT_datetime' ## non_ST_IPO_NT_SmallVolume_datetime
    }

    return options
def VarsMeasureWindow():# 变量测度窗口期
    # TODO与InputOptions合并
    '''measure_window=VarsMeasureWindow()'''
    # 如果'calendar':'calendar'——》例如MonthEnd(12)或DateOffset(years=0,months=1,weeks=0,days=0)
    # 如果'calendar':'trading'——》30(>=0) ，也就是说数字搭配trading设置
    options={
        'beta':(DateOffset(months=12),Day(0)), #(45,0),#DateOffset(years=0,months=1,weeks=0,days=0) #
        'mom':(DateOffset(months=12), DateOffset(months=1)),#MonthEnd(12)
        'rev':( DateOffset(months=1),Day(0)),
        'illiq':(DateOffset(months=12),Day(0)),
        'turnover':(DateOffset(months=2),Day(0)),
        'max_ret':(DateOffset(months=2),Day(0)),
        'skew':(DateOffset(months=12),Day(0)),
        'coskew':(DateOffset(months=12),Day(0)),
        'iskew':(DateOffset(months=2),Day(0)),
        'vol':(DateOffset(months=12),Day(0)),
        'ivol':(DateOffset(months=1),Day(0))
    }
    return options
def EndDateOfPerInvestment():
    # calendar='calendar'
    #       fix_date: num; fix_num: 推荐使用对应日频freq，例如“15D”
    # calendar='trading'
    #
    if options['calendar']== 'calendar':
        tmp=pd.date_range(start=options['start'], end=options['end'], freq=options['calendar_periods']) - Day(options['num_reverse'])
        return tmp[tmp>options['data_start']]
    elif options['calendar'] == 'trading':
        try:
            return pd.read_pickle(options['data_path'] + 'trading_calendar')[options['num_periods']::options['num_periods']]
        except FileNotFoundError:
            trddt=GetTradingCalendar()
            pd.to_pickle(trddt, options['data_path'] + 'trading_calendar')
            return trddt[options['num_periods']::options['num_periods']]
def GetTradingCalendar():
    trddt = pd.read_pickle(options['data_path'] + 'index_ret')[['index_code','trddt']].set_index(['index_code','trddt']).loc['000016.SH'].index
    trddt = pd.to_datetime(trddt.astype(int).astype(str), format='%Y%m%d')
    return trddt
def PerCal_Start2End(EndDate,var_name):
    '''
    用于计算每期变量测度窗口期的start和end
    :param EndDate:
    :param var_name:
    :return:
    '''
    tmp=DataFrame(index=EndDate,columns=['start','end'])
    if not isinstance(measure_window[var_name][0],int):
        tmp['end'] = EndDate - measure_window[var_name][1]
        tmp['start']=EndDate-measure_window[var_name][0]+Day()
    else:
        trading_calendar = GetTradingCalendar()
        for t in tmp.index[tmp.index>=trading_calendar[measure_window[var_name][0]-1]]:
            tmp.loc[t, 'end'] = trading_calendar[trading_calendar <=t][-measure_window[var_name][1]-1]
            tmp.loc[t,'start']=trading_calendar[trading_calendar<=t][-measure_window[var_name][0]]
    tmp=tmp[tmp['start'] >= options['start']] # 如果是nan，会自动删除
    tmp[tmp<options['data_start']]=options['data_start']
    return tmp
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
    data_path=options['data_path']
    if PV_vars is not None:
        if len(PV_vars)>1:
            PV=pd.read_pickle(data_path+options['pv_name'])[PV_vars]
        else:
            PV = pd.read_pickle(data_path+options['pv_name'])[PV_vars[0]]
        if filter_data is True:
            filter_ =pd.read_pickle(data_path + options['filter_data_name']).stack()
            PV=PV[filter_]
    else:
        PV=None
    if BS_vars is not None:
        BS = pd.read_pickle(data_path + options['bs_name'])[[options['bs_date_name'], 'stkcd']+BS_vars]  #'BS'
        BS = BS[~(BS[options['bs_date_name']].isnull())]
        BS[options['bs_date_name']] = pd.to_datetime(BS[options['bs_date_name']].astype(int).astype(str), format='%Y%m%d')
        if len(BS_vars)>1:
            BS=BS.set_index([options['bs_date_name'], 'stkcd']).sort_index()
        else:
            BS = BS.set_index([options['bs_date_name'], 'stkcd']).sort_index()[BS_vars[0]]
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
def PrepareData(var_list):
    # TODO 变量归类，以便输入数据
    # 全局变量：options,measure_window
    # EndDateIndex
    # ret_d: beta
    #       PV[adj_open,adj_close]
    # index_ret_d: beta
    #       index_ret/PV[adj_open,adj_close,size_tot]/
    # size_d: size,BM
    # book: BM
    #
    pass
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
def AlignData(*vars):
    l=len(vars)
    index=vars[0].index.copy()
    columns=vars[0].columns.copy()
    index_not_null=~vars[0].isnull().all(axis=1)
    for i in range(1,l):
        index &=vars[i].index
        columns &=vars[i].columns
        index_not_null &=(~vars[i].isnull().all(axis=1))
    return (vars[i].reindex(index=index[index_not_null],columns=columns) for i in range(l))
def GetValueWeightedIndexReturn(p0,p1,size,EndDate):
    EndDateIndex = EndDate
    weights=size.resample('D').ffill().loc[EndDateIndex].shift(1)#.stack()
    if options['open_price']=='adj_open':
        price0=p0.resample('D').bfill().shift(-1).loc[EndDateIndex].shift(1)#.stack()
    elif options['open_price']=='adj_pre_close':
        price0=p0.resample('D').ffill().loc[EndDateIndex].shift(1)
    price1=p1.resample('D').ffill().loc[EndDateIndex]#.stack()
    valid=~(price1+price0+size).isnull()
    return ((price1-price0)/price0*size).sum(axis=1)/size[valid].sum(axis=1)*100.0
def GetDailyValueWeightedIndexReturn(p0,p1,size):
    valid=~(size+p0+p1).isnull()
    return (size*(p1-p0)/p0).sum(axis=1)/size[valid].sum(axis=1)*100.0
def cal_ret(p0,p1,EndDate,del_Rf=True):
    # TODO: 剔除无风险收益率
    # try:
    #     p0.iloc[0,0]
    # except (NameError,AttributeError):
    #     pv=import_data(PV_vars=['adj_open','adj_close'])[0]
    #     p0=pv['adj_open'].unstack()
    #     p1=pv['adj_close'].unstack()
    price0=p0.resample('D').bfill().shift(-1).loc[EndDate].shift(1)
    price1=p1.resample('D').ffill().loc[EndDate]
    ret= (price1-price0)/price0*100.0
    ret = ret[ret != 0.0]
    if del_Rf:
        Rf = import_data(Rf_freq=options['freq'])[2]
        ret = (ret + 100.0).div((Rf[Rf.columns[0]].reindex(index=ret.index) * 0.01 + 1.0) ** options[
            'num_periods'],axis=0) - 100.0
    return ret.astype(float)
def cal_ret_d(p0,p1,del_Rf=True):
    '''
    :param p0:
        unstack 后的 open_price数据
    :param p1:
        unstack 后的 close_price数据
    :return:
    '''
    ret = (p1 - p0) / p0*100.0
    ret=ret[ret != 0.0]
    if del_Rf:
        Rf = import_data(Rf_freq='D')[2]
        ret = (ret + 100.0).div(Rf[Rf.columns[0]].reindex(index=ret.index) * 0.01 + 1.0,axis=0) - 100.0
    return ret.astype(float)
def cal_index_ret(index_code=None,del_Rf=True):
    # TODO 有问题
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
    data_path = options['data_path']
    if index_code is None:
        index_ret=GetValueWeightedIndexReturn(p0,p1,size,EndDate)
    else:
        index_ret = pd.read_pickle(data_path + 'index_ret')[['index_code', 'trddt', 'opnprc', 'clsprc']]
        index_ret = index_ret[index_ret['index_code'] == index_code][['trddt', 'opnprc', 'clsprc']].set_index('trddt')
        index_ret.index = pd.to_datetime(index_ret.index.astype(int).astype(str), format='%Y%m%d')
        EndDateIndex=EndDateOfPerInvestment()
        index_ret=cal_ret(index_ret['opnprc'],index_ret['clsprc'],EndDateIndex)

    if del_Rf:
        Rf = import_data(Rf_freq=options['freq'])[2]
        index_ret = (index_ret + 100.0) / (Rf[Rf.columns[0]].reindex(index=ret.index) * 0.01 + 1.0)**options['num_periods'] - 100.0
    return index_ret.astype(float)
def cal_index_ret_d(index_code=None,del_Rf=True):
    data_path = options['data_path']
    if index_code is None:
        try:
            index_ret=pd.read_pickle(data_path+'daily_value_weighted_index_return')
        except FileNotFoundError:
            index_ret = GetDailyValueWeightedIndexReturn(p0,p1,size)
            index_ret.to_pickle(data_path + 'daily_value_weighted_index_return')
    else:
        index_ret = pd.read_pickle(data_path + 'index_ret')[['index_code', 'trddt', 'opnprc', 'clsprc']]
        index_ret = index_ret[index_ret['index_code'] == index_code][['trddt', 'opnprc', 'clsprc']].set_index('trddt')
        index_ret.index = pd.to_datetime(index_ret.index.astype(int).astype(str), format='%Y%m%d')
        #EndDateIndex = EndDateOfPerInvestment()
        index_ret =(index_ret['clsprc']-index_ret['opnprc'])/index_ret['opnprc']*100.0
    if del_Rf:
        Rf = import_data(Rf_freq='D')[2]
        index_ret = (index_ret + 100.0) / (Rf[Rf.columns[0]].reindex(index=index_ret.index) * 0.01 + 1.0)- 100.0
    return index_ret.astype(float)
def cal_beta(index_ret_d,ret_d,EndDate):
    start2end=PerCal_Start2End(EndDate,'beta')
    beta = pd.DataFrame(index=start2end.index, columns=ret_d.columns)
    for t in beta.index:
        tmp1 = index_ret_d.loc[start2end.loc[t,'start']:start2end.loc[t,'end']] \
               - index_ret_d.loc[start2end.loc[t,'start']:start2end.loc[t,'end']].mean()
        tmp2 = ret_d.loc[tmp1.index] \
               - ret_d.loc[tmp1.index].mean()
        beta.loc[t] = tmp2.mul(tmp1, axis=0).mean() / tmp1.var()
    return beta[beta!=0.0].shift(1).iloc[1:].astype(float)
def cal_size_BM(size_d,EndDate,shift_num=1,only_size=False):
    # TODO 当BM频率较高时，相对于FF的分组方式，分组信息中增加了动量或者反转以及五因子中的investment因子
    size = size_d.resample('D').ffill().loc[EndDate]* 1.0e-5
    if only_size:
        return size[size > 0.0].shift(shift_num).iloc[shift_num:]
    else:
        book = import_data(BS_vars=['tot_shrhldr_eqy_excl_min_int'])[1]
        BM = book.drop(book.index[book.index.duplicated(keep='last')]).unstack().resample('D').first().ffill()\
            .reindex(index=EndDate, columns=size.columns)/ size * 1.0e-8
        return size[size > 0.0].shift(shift_num).iloc[shift_num:],BM[BM > 0.0].shift(shift_num).iloc[shift_num:]
def cal_mom(p0,p1,EndDate):
    start2end = PerCal_Start2End(EndDate, 'mom')
    price0=p0.resample('D').bfill().reindex(index=start2end['start'])
    price0.index=start2end.index
    price1=p1.resample('D').ffill().reindex(index=start2end['end'])
    price1.index = start2end.index
    return ((price1-price0)/price0).shift(1).iloc[1:]*100.0
def cal_rev(p0,p1,EndDate):
    start2end = PerCal_Start2End(EndDate, 'rev')
    price0 = p0.resample('D').bfill().reindex(index=start2end['start'])
    price0.index = start2end.index
    price1 = p1.resample('D').ffill().reindex(index=start2end['end'])
    price1.index = start2end.index
    return ((price1 - price0) / price0).shift(1).iloc[1:] * 100.0
def cal_illiq(ret_d,amount_d,EndDate):
    start2end = PerCal_Start2End(EndDate, 'illiq')
    illiq_d=1e3*np.abs(ret_d) / amount_d
    illiq=pd.DataFrame(index=start2end.index, columns=ret_d.columns)
    for t in illiq.index:
        illiq.loc[t]=illiq_d.loc[start2end.loc[t,'start']:start2end.loc[t,'end']].mean()
    return illiq.shift(1).iloc[1:].astype(float)
def cal_turnover(amount_d,size_free_d,EndDate):
    turnover_d=amount_d/size_free_d*10.0
    #turnover=turnover[turnover>0.0001]
    start2end = PerCal_Start2End(EndDate, 'turnover')
    turnover= pd.DataFrame(index=start2end.index, columns=ret_d.columns)
    for t in turnover.index:
        turnover.loc[t]=turnover_d.loc[start2end.loc[t,'start']:start2end.loc[t,'end']].mean()
    return turnover.shift(1).iloc[1:].astype(float)
def cal_max_ret(ret_d,EndDate):
    '''
    也许15个交易日的效果更好
    :param periods:
    :param freq:
    :param ret_d:
    :return:
    '''
    start2end = PerCal_Start2End(EndDate, 'max_ret')
    max_ret= pd.DataFrame(index=start2end.index, columns=ret_d.columns)
    for t in max_ret.index:
        max_ret.loc[t] = ret_d.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']].max()
    return max_ret.shift(1).iloc[1:].astype(float)
def cal_skew(ret_d,EndDate):
    start2end = PerCal_Start2End(EndDate, 'skew')
    skew=pd.DataFrame(index=start2end.index,columns=ret_d.columns)
    for t in skew.index:
        skew.loc[t]=ret_d.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']].skew()
    return skew.shift(1).iloc[1:].astype(float)
def cal_coskew(ret_d,index_ret_d,EndDate):
    X=pd.DataFrame({'index':index_ret_d,'index^2':index_ret_d**2.0})[['index','index^2']]#.astype(float)
    start2end = PerCal_Start2End(EndDate, 'coskew')
    coskew=pd.DataFrame(index=start2end.index,columns=ret_d.columns)
    for t in coskew.index:
        tmp1 = X.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']] \
               - X.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']].mean()
        tmp2 = ret_d.loc[tmp1.index] \
               - ret_d.loc[tmp1.index].mean()
        tmp3 = np.linalg.pinv(tmp1.cov())[1]
        coskew.loc[t] = tmp3[0] * tmp2.mul(tmp1['index'], axis=0).mean() \
                          + tmp3[1] * tmp2.mul(tmp1['index^2'],axis=0).mean()
    return coskew.shift(1).iloc[1:].astype(float)*100.0
def cal_SMB_HML(ret,size,BM,percentile1=None,percentile2=None,independent=True,weighted=True):
    ret,size,BM=AlignData(ret,size,BM)
    valid_ = ~pd.isnull(BM+ret+size)  # TypeError: bad operand type for unary ~: 'float'--->index或columns不匹配
    #valid_=valid_
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
    if weighted:
        df1 = DataFrame()
        df2 = DataFrame()
        df1['rtn_w'] = (ret*size).stack()
        df1['ref1'] = mark_1.stack()
        df1['ref2'] = mark_2.stack()
        df2['rtn_w']=size.stack()
        df2['ref1']=mark_1.stack()
        df2['ref2']=mark_2.stack()
        tmp1 = df1.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).sum())
        tmp2 = df2.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).sum())
        tmp = (tmp1 / tmp2).unstack()['rtn_w']
        # tmp.index.names = ('trddt', 'ref1')
        # HML = tmp.mean(axis=0, level=0)
        # SMB = tmp.mean(axis=1).unstack()
    else:
        df = DataFrame()
        df['rtn'] = ret.stack()
        df['ref1'] = mark_1.stack()
        df['ref2'] = mark_2.stack()
        tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).mean()).unstack()['rtn']
        # tmp.index.names = ('trddt', 'ref1')
    HML = tmp.mean(axis=0, level=0)
    SMB = tmp.mean(axis=1).unstack()
    return SMB.iloc[:,-1]-SMB.iloc[:,0], HML.iloc[:,-1]-HML.iloc[:,0]
def cal_SMB_HML_FF(ret,EndDate,size=None,book=None,weights=None):
    # TODO return要和EndDate的频率保持一致
    percentile1 = [0.0, 0.5, 1.0]  # size
    percentile2 = [0.0, 0.3, 0.7, 1.0]  # value
    label_1 = [i + 1 for i in range(len(percentile1) - 1)]
    label_2 = [i + 1 for i in range(len(percentile2) - 1)]
    size, book = import_data(PV_vars=['size_tot'], BS_vars=['tot_shrhldr_eqy_excl_min_int'])[:2]
    BE = book.drop(book.index[book.index.duplicated(keep='last')]).unstack()
    BE = BE[BE.index.month == 12]
    BE = BE[BE > 0]
    size = size.unstack()
    ME = size.copy()
    ME = ME.resample('M').last()
    ME6 = ME[ME.index.month == 6]
    ME12 = ME[ME.index.month == 12]
    ME12.loc[parse('20041231')] = size.loc['2005-01-04']
    ME12 = ME12.sort_index()
    BM = BE.reindex(index=ME12.index, columns=ME12.columns) / ME12
    mark_1 = DataFrame([qcut(ME6.loc[i], q=percentile1, labels=label_1) for i in ME6.index])
    mark_1.index = mark_1.index + Day()
    mark_1 = mark_1.resample('D').ffill().reindex(index=EndDate)
    mark_2 = DataFrame([qcut(BM.loc[i], q=percentile2, labels=label_2) for i in BM.index])
    mark_2.index = mark_2.index + MonthBegin(7)
    mark_2 = mark_2.resample('D').ffill().reindex(index=EndDate)

    if weights is None:
        df = DataFrame()
        df['ret'] = ret.stack()
        df['ref1'] = mark_1.stack()
        df['ref2'] = mark_2.stack()
        df = df.dropna()
        tmp = df.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).mean())['ret'].unstack()
    else:
        weights = size.resample('D').ffill().reindex(index=EndDate).shift(1)
        df = DataFrame()
        df['ret'] = (ret * weights).stack()
        df['ref1'] = mark_1.stack()
        df['ref2'] = mark_2.stack()
        df['w'] = weights.stack()
        df = df.dropna()
        tmp1 = df.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).sum())['ret']
        tmp2 = df.groupby(level=0).apply(lambda g: g.groupby(['ref1', 'ref2']).sum())['w']
        tmp = (tmp1 / tmp2).unstack()
    rHML = tmp.mean(axis=0, level=0)
    rSMB = tmp.mean(axis=1).unstack()
    return rSMB.iloc[:, -1] - rSMB.iloc[:, 0], rHML.iloc[:, -1] - rHML.iloc[:, 0]
def cal_ivol_from_EGARCH():
    return pd.read_pickle(GetDataPath()+'Expected_iVol_using_OptimalEgarch').shift(1).loc['2007-09']
def cal_iskew(ret_d,index_ret_d,EndDate,SMB_d=None,HML_d=None,method='FF'):
    # todo 待检验
    start2end = PerCal_Start2End(EndDate, 'iskew')
    iskew = pd.DataFrame(index=start2end.index, columns=ret_d.columns)
    if method == 'CAPM':
        for t in iskew.index:
            tmp1 = index_ret_d.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']] \
                   - index_ret_d.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']].mean()
            tmp2 = ret_d.loc[tmp1.index] \
                   - ret_d.loc[tmp1.index].mean()
            iskew.loc[edt] = (tmp2 - tmp1[:, None] @ tmp2.mul(tmp1, axis=0).mean()[None, :] / tmp1.var()).skew()
    elif method=='FF':
        X = pd.DataFrame({'index': index_ret_d, 'SMB': SMB_d, 'HML': HML_d})[['index', 'SMB', 'HML']].dropna()
        for t in iskew.index:
            tmp1 = X.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']] \
                   - X.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']].mean()
            tmp2 = ret_d.loc[tmp1.index] \
                   - ret_d.loc[tmp1.index].mean()
            iskew.loc[edt] = (tmp2 - tmp1 @ np.linalg.pinv(
                tmp1.cov()) @ np.array([
                tmp2.mul(tmp1['index'], axis=0).mean(),
                tmp2.mul(tmp1['SMB'], axis=0).mean(),
                tmp2.mul(tmp1['HML'], axis=0).mean()
            ])).skew()
    return iskew.shift(1).iloc[1:].astype(float)
def cal_vol(ret_d,EndDate,ZeroMean=False):
    # TODO 待检验
    start2end = PerCal_Start2End(EndDate, 'vol')
    #vol = pd.DataFrame(index=start2end.index, columns=ret_d.columns)
    if ZeroMean:
        ret_d = ret_d ** 2.0
        vol=np.sqrt(DataFrame(
            (ret_d.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']].mean() for t in
             start2end.index)
            , index=start2end.index).shift(1).iloc[1:])
    else:
        vol=pd.DataFrame(
            (ret_d.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']].std() for t in start2end.index)
            , index=start2end.index).shift(1).iloc[1:]
    return vol[vol>0.0]
def cal_ivol(ret_d,index_ret_d,EndDate,SMB_d=None,HML_d=None,method='FF'):
    # todo 待检验;对于月度数据呢？？？
    start2end = PerCal_Start2End(EndDate, 'ivol')
    ivol = pd.DataFrame(index=start2end.index, columns=ret_d.columns)
    if method == 'CAPM':
        for t in ivol.index:
            tmp1 = index_ret_d.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']] \
                   - index_ret_d.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']].mean()
            tmp2 = ret_d.loc[tmp1.index] \
                   - ret_d.loc[tmp1.index].mean()
            ivol.loc[t] = (tmp2 - tmp1[:, None] @ tmp2.mul(tmp1, axis=0).mean()[None, :] / tmp1.var()).std()
    elif method=='FF':
        X=pd.DataFrame({'index':index_ret_d,'SMB':SMB_d,'HML':HML_d})[['index', 'SMB', 'HML']].dropna()
        for t in ivol.index:
            tmp1 = X.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']] \
                   - X.loc[start2end.loc[t, 'start']:start2end.loc[t, 'end']].mean()
            tmp2 = ret_d.loc[tmp1.index] \
                   - ret_d.loc[tmp1.index].mean()
            ivol.loc[t] = (tmp2 - tmp1 @ np.linalg.pinv(
                tmp1.cov()) @ np.array([
                tmp2.mul(tmp1['index'], axis=0).mean(),
                tmp2.mul(tmp1['SMB'], axis=0).mean(),
                tmp2.mul(tmp1['HML'], axis=0).mean()
            ])).std()
    return ivol[ivol>0.0].shift(1).iloc[1:].astype(float)
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


options=InputOptions()
EndDate=EndDateOfPerInvestment()
measure_window=VarsMeasureWindow()
PV=import_data(PV_vars=['size_tot',options['open_price'],options['close_price']])[0]
p0=PV[options['open_price']].unstack()
p1=PV[options['close_price']].unstack()
size=PV['size_tot'].unstack()
ret_d=cal_ret_d(p0,p1)
ret=cal_ret(p0,p1,EndDate)
# def acovf(arr,lag):
#     '''
#     支持 1-D和2-D；unbiased
#     :param arr:
#     :param lag:
#     :return:
#     '''
#     return np.nanmean((arr[lag:len(arr)] - np.mean(arr))*(arr[0:len(arr)-lag] - np.mean(arr)),axis=0)
# rho=((ret_d.iloc[:-1]-ret_d.mean()).shift(1)*(ret_d.iloc[1:]-ret_d.mean())).mean()/ret_d.var()
# from my_modules.TSAtools import pacfs
# ret_d_dm=ret_d-ret_d.mean()
# rho_d=(ret_d_dm.shift(1)*ret_d_dm).mean()/ret_d_dm.var()
# ret_dm=ret-ret.mean()
# rho=(ret_dm.shift(1)*ret_dm).mean()/ret_dm.var()

# ret_d.std().mean()*np.sqrt(220)
# ret.std().mean()*np.sqrt(12)

mp_d=[1,3,6,12]
mp_m=[12,24,36]
index_ret_d=cal_index_ret_d()
index_ret=cal_index_ret(index_code='399300.SZ')
#SMB_d,HML_d=cal_SMB_HML_FF(ret_d,ret_d.index)
size_d,BM_d=cal_size_BM(size,size.index)
size_m,BM_m=cal_size_BM(size,EndDate)
SMB_d,HML_d=cal_SMB_HML(ret_d,size_d,BM_d)
SMB_m,HML_m=cal_SMB_HML(ret,size_m,BM_m)
ivol_method='FF'
vol_name='vol' # 'vol_ss'
var_list=[]
for i in mp_d:
    var_list.append(vol_name+'_'+str(i)+'M')
for i in mp_m:
    var_list.append(vol_name+'_'+str(i//12)+'Y')

sum_stat=pd.DataFrame(index=var_list,columns=['count', 'mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max',
       'skew', 'kurt'])
for i in mp_d:
    #options['calendar_periods']=MonthEnd(i)
    #options['num_periods']=float(i)
    #EndDate=EndDateOfPerInvestment()
    measure_window[vol_name]=(MonthEnd(i),Day(0))
    if vol_name.startswith('vol'):
        vol = cal_vol(ret_d, EndDate,ZeroMean=True) * np.sqrt(220)
    else:
        vol=cal_ivol(ret_d,index_ret_d,EndDate,SMB_d=SMB_d,HML_d=HML_d,method=ivol_method)* np.sqrt(220)
    sum_stat.loc[vol_name+'_'+str(i)+'M']=describe(vol.T).mean(axis=1)
for i in mp_m:
    measure_window[vol_name] = (MonthEnd(i), Day(0))
    if vol_name.startswith('vol'):
        vol = cal_vol(ret, EndDate,ZeroMean=True) * np.sqrt(12)
    else:
        vol=cal_ivol(ret,index_ret,EndDate,SMB_d=SMB_m,HML_d=HML_m,method=ivol_method)* np.sqrt(12)
    sum_stat.loc[vol_name+'_'+str(i//12)+'Y'] = describe(vol.T).mean(axis=1)
sum_stat.to_csv(options['data_path']+vol_name+'_'+ivol_method+'.csv')
sum_stat
