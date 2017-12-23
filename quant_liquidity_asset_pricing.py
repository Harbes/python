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
    data = pd.read_pickle(data_path + '/' + data_set)[Vars]
    global amount, size_free, size_tot, adj_open, adj_close
    amount = data['amount'].unstack()
    size_free = data['size_free'].unstack()
    size_tot = data['size_tot'].unstack()
    adj_open = data['adj_open'].unstack()
    adj_close = data['adj_close'].unstack()


def weighted_rtn(rtn, weights=None):
    if weights is None:
        return rtn.mean(axis=1)
    else:
        rtn_w = rtn * weights
        rtn_w = rtn_w[~np.isnan(weights)]
        return rtn_w.sum(axis=1) / (weights[~np.isnan(rtn_w)].sum(axis=1))


def resample_data():
    key = lambda x: x.year * 100 + x.month
    global adj_open_m, adj_close_m, size_tot_m, turn_over_m
    adj_open_m = adj_open.groupby(key).first()
    adj_close_m = adj_close.groupby(key).last()
    size_tot_m = size_tot.groupby(key).last().shift(1).iloc[1:]
    turn_over_m = (amount / size_free).groupby(key).sum()
    # adj_open_m,adj_close_m,size_tot_m,turn_over_m


def mark_group(indicator1, indicator2=None):
    global mark_1, labels1
    label_num1 = 5;
    labels1 = [i + 1 for i in range(label_num1)];
    percentiles1 = np.linspace(0, 1, label_num1 + 1)
    mark_1 = DataFrame([pd.qcut(indicator1.loc[i], q=percentiles1, labels=labels1) for i in indicator1.index])
    if indicator2 is not None:
        global mark_2, labels2
        label_num2 = 5;
        labels2 = [i + 1 for i in range(label_num2)];
        percentiles2 = np.linspace(0, 1, label_num2 + 1)
        mark_2 = DataFrame(np.nan, index=mark_1.index, columns=mark_1.columns)
        for l_ in labels2:
            tmp = DataFrame([pd.qcut(indicator2.loc[i][mark_1.loc[i] == l_],
                                     q=percentiles2, labels=labels2) for i in indicator2.index],
                            index=indicator2.index)
            mark_2 = mark_2.combine_first(tmp)


def var_by_groups(var, weights=None):
    return DataFrame([[var.loc[i][mark_1.loc[i] == l_].mean() for l_ in labels1] for i in var.index], index=var.index,
                     columns=labels1)


if __name__ is "__main__":
    import_data()
    size_tot.loc[pd.to_datetime('2004-12-30')] = size_tot.iloc[0];
    size_tot = size_tot.sort_index()
    resample_data()
    mark_group(size_tot_m)
    rtn_m = (adj_close_m - adj_open_m) / adj_open_m
    rtn_by_size = var_by_groups(rtn_m)
    TurnOver_by_size = var_by_groups(turn_over_m)
    tmp = TurnOver_by_size[1] - TurnOver_by_size[5];
    tmp.mean() / tmp.std() * np.sqrt(len(tmp))
