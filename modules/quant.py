import pandas as pd
import numpy as np
from pandas.tseries.offsets import MonthEnd
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
        *_,Rf=import_data(Rf_freq=freq)
        #index_ret = index_ret.sub(Rf.loc[index_ret.index], axis=0) # 报错，Rf是一个DataFrame，而不是Series
        #index_ret=-Rf.loc[index_ret.index].sub(index_ret,axis=0) # 虽然成功，但是返回的是一个DataFrame，而不是Series
        #index_ret = index_ret.sub(Rf.loc[index_ret.index].iloc[:,0], axis=0)
        index_ret=(index_ret + 100.0) / (Rf.loc[index_ret.index].iloc[:, 0] * 0.01 + 1.0) - 100.0

    return index_ret
def cal_ret(periods=1,freq='M',del_Rf=True):
    # TODO 尚未剔除异常交易数据
    '''

    calculate the past returns,including the current period
    :param periods:
    :param freq:
    :param del_Rf:
    :return:
    '''
    p0 = 'adj_open';p1 = 'adj_close' # price0也可以使用 'adj_pre_close'
    PV,*_=import_data(PV_vars=[p0,p1])
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
        *_, Rf = import_data(Rf_freq=freq)
        ret = (ret.astype(float) + 100.0).div((Rf.loc[ret.index].iloc[:, 0] * 0.01 + 1.0) ** periods, axis=0) - 100.0
    return ret




def cal_variables(var_list,start='20050101',end='20180228'):
    pass

def cal_beta(periods,index_ret,ret_d):
    '''
    calculate the market beta
    :param periods: int{measure periods;number of month}
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





data_dict={'pv':['']}

