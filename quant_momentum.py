# 中国市场似乎只存在短期反转效应：Monthly: J/M=1,1或者J/M=(1,2,4,9,10不是很强烈,因为VW不显著；3,5,6,7,8,11,12),0
#                            daily:  J/M=5,1
# 总的结果是short-term reversal，但是似乎也存在momentum，例如 第10组比第9组更高(equally-weighted)——》rtn.quantile(0.95,axis=1).mean()；
import pandas as pd
import numpy as np
from pandas import DataFrame
import time


def import_data():
    global price0, price1, adj_open, adj_close, size
    # data_path='F:/data/xccdata'
    data_path = '/Users/harbes/data/xccdata'
    # data=pd.read_pickle(data_path+'/PV_datetime')[['adj_open','adj_close']]
    data = pd.read_pickle(data_path + '/PV_datetime')[['adj_open', 'adj_close', 'size_tot']]  # 'size_free'
    price0 = data['adj_open'].unstack()
    price1 = data['adj_close'].unstack()
    size = data['size_tot'].unstack()  # size=data['size_free'].unstack()


def rtn_monthly():
    global price0, price1, rtn, size
    key = lambda x: x.year * 100 + x.month
    price0 = price0.groupby(key).first()  # 买入价
    price1 = price1.groupby(key).last()  # 卖出价
    price0.index = pd.to_datetime(price0.index.values.astype(str), format=('%Y%m'))
    price1.index = pd.to_datetime(price1.index.values.astype(str), format=('%Y%m'))
    size = size.groupby(key).last()
    size.index = pd.to_datetime(size.index.values.astype(str), format=('%Y%m'))
    rtn = ((price1 - price0) / price0).iloc[J + M:]
    rtn[rtn == 0] = np.nan

def rtn_daily():
    global price0, price1, rtn
    price0 = adj_open
    price1 = adj_close
    rtn = ((price1 - price0) / price0).iloc[J + M:]
    rtn[rtn == 0] = np.nan


def group_by_mom():
    global percentile, label_momentum, mark_momentum
    percentile = np.linspace(0, 1, group_num + 1)
    label_momentum = [i + 1 for i in range(group_num)]  # 1表示过去收益较差的组合，10表示过去收益较好的组合
    past_rtn = np.log(price1.shift(M + 1) / price0.shift(M + J));
    past_rtn[past_rtn == 0] = np.nan
    mark_momentum = DataFrame([pd.qcut(past_rtn.iloc[i],
                                       q=percentile, labels=label_momentum) for i in range(J + M, len(price0))],
                              index=price0.index[J + M:])

def rtn_EW_by_mom_JegadeeshTitman1993():
    # percentile = np.linspace(0, 1, group_num + 1)
    # label_momentum = [i + 1 for i in range(group_num)]  # 1表示过去收益较差的组合，10表示过去收益较好的组合
    # mark_momentum = DataFrame([pd.qcut((price1.iloc[i - 1 - M] - price0.iloc[i - J - M]) / price0.iloc[i - J - M],
    #                                   q=percentile, labels=label_momentum) for i in range(J + M, len(price0))],
    #                          index=price0.index[J + M:])
    momentum = DataFrame(
        [[rtn.iloc[i][mark_momentum.iloc[i] == k].mean() for k in label_momentum] for i in range(K - 1, len(rtn))],
        index=rtn.index[K - 1:], columns=label_momentum)
    if K > 1:
        for lag in range(1, K):
            tmp = DataFrame([[rtn.iloc[i][mark_momentum.iloc[i - lag] == k].mean() for k in label_momentum] for i in
                             range(K - 1, len(rtn))], index=rtn.index[K - 1:], columns=label_momentum)
            momentum += tmp
    momentum /= K
    return momentum

def rtn_EW_by_mom():
    momentum = DataFrame(
        [[rtn.loc[i][mark_momentum.loc[i] == k].mean() for k in label_momentum] for i in rtn.index],
        index=rtn.index, columns=label_momentum)
    return momentum


def rtn_VW_by_mom():
    global size
    size = size.shift(1)
    rtn_size=rtn*size
    momentum = DataFrame(
        [[rtn_size.loc[i][mark_momentum.loc[i] == k].sum() / size.loc[i][mark_momentum.loc[i] == k].sum() for k in
          label_momentum] for i in rtn.index],
        index=rtn.index, columns=label_momentum)
    return momentum


if __name__ == "__main__":
    t0 = time.time()
    group_num = 10
    J = 3  # past J-period return
    M = 0  # time lapse
    K = 3  # periods of overlapping ( JEGADEESH and TITMAN,1993)
    import_data()
    rtn_monthly()  # rtn_daily()  #
    group_by_mom()
    # rtn_JT1993 = rtn_EW_by_mom_JegadeeshTitman1993();
    rtn_EW = rtn_EW_by_mom();
    rtn_VW = rtn_VW_by_mom();
    print('running time:', time.time() - t0)
    print('GroupNums:', group_num, '\n', 'J/M:{a}/{b}'.format(a=J, b=M),'\n')
    # tmp1 = rtn_JT1993[group_num] - rtn_JT1993[1];print('JT1993:', '\n', rtn_JT1993.mean(), tmp1.mean(), '({})'.format(tmp1.mean() / tmp1.std() * np.sqrt(len(tmp1))),'\n')
    tmp2 = rtn_EW[group_num] - rtn_EW[1];
    print('Equal_weighted:', '\n', rtn_EW.mean(), tmp2.mean(),
          '({})'.format(tmp2.mean() / tmp2.std() * np.sqrt(len(tmp2))), '\n')
    tmp3 = rtn_VW[group_num] - rtn_VW[1];
    print('Value_weighted:', '\n', rtn_VW.mean(), tmp3.mean(),
          '({})'.format(tmp3.mean() / tmp3.std() * np.sqrt(len(tmp3))), '\n')

    ## 输入市场回报数据
    # rtn_index=pd.read_pickle('/Users/harbes/data/xccdata/essay/index_hs300_monthly')[f_momen.index] #直接用市场的数据似乎并不能说明 cov(Rm_t,Rm_t_1)正比于 cov(f_t,f_t_1)
    # rtn_index = rtn.mean(axis=1)[f_momen.index]
    # np.cov(rtn_index[1:], rtn_index[:-1]) # ( JEGADEESH and TITMAN,1993)
    # np.cov(f_momen[1:], f_momen[:-1])
