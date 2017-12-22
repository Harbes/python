import pandas as pd
from pandas import DataFrame
import numpy as np
from pandas.tseries.offsets import YearEnd
from datetime import datetime


def import_data():
    '''

    :param data_path: str
    :param data_set: str
    :param Vars: str or list of str
    :return:
    '''
    data_path = '/Users/harbes/data/xccdata'
    # data_path='F:/data/xccdata'
    data_set = 'PV_datetime'
    Vars = ['size_tot', 'size_free', 'adj_close', 'adj_open', 'amount']
    return pd.read_pickle(data_path+'/'+data_set)[Vars]

def weighted_rtn(rtn,weights=None):
    if weights is None:
        return rtn.mean(axis=1)
    else:
        rtn_w=rtn*weights
        rtn_w=rtn_w[~np.isnan(weights)]
        return rtn_w.sum(axis=1)/(weights[~np.isnan(rtn_w)].sum(axis=1))

if __name__ is '__main__':
    data=import_data()
    adj_open=data['adj_open'].unstack()
    adj_close=data['adj_close'].unstack()
    size=data['size_tot'].unstack()
    rtn=(adj_close-adj_open)/adj_open
    auto_corrs = rtn.corrwith(rtn.shift(1))
    t_stat = auto_corrs.mean() / auto_corrs.std() * len(auto_corrs)
    print(t_stat)
    rtn_vw=weighted_rtn(rtn,weights=size)
    print(rtn_vw.corr(rtn_vw.shift(1)))
    key=lambda x:x.year*100+x.month
    adj_open_m=adj_open.groupby(key).first();
    adj_close_m=adj_open.groupby(key).last()
    size_m = size.groupby(key).first()
    rtn_m=(adj_close_m-adj_open_m)/adj_open_m
    rtn_m_vw=weighted_rtn(rtn_m,weights=size_m)
    print(rtn_m_vw.corr(rtn_m_vw.shift(1)))

#run quant_beta_and_risk_premium.py
