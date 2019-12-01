import pandas as pd
import numpy as np
from pandas.tseries.offsets import DateOffset
import scipy
import matplotlib.pyplot as plt
DPath='/Users/harbes/PycharmProjects/data/CNRDS/' #'E:/data/CNRDS/' #
indicators_all=pd.read_pickle(DPath+'indicators_all').replace([np.inf, -np.inf], np.nan)
PV=pd.read_pickle(DPath+'PV');PV.tail(10)
tmp=PV[['Scode','Trddt','Adclsprc']].iloc[1:].set_index(['Trddt','Scode'])
tmp=tmp.astype(float)
adjprc=tmp.loc[~tmp.index.duplicated(keep='last')].unstack() # 提出duplicates
adjprc.index=pd.to_datetime(adjprc.index,format='%Y-%m-%d')# 调整index格式
adjprc.columns=adjprc.columns.get_level_values(1)# 调整columns格式
GroupBy1=lambda x:x.year*100.0+x.month
#p0=adjprc.groupby(by=GroupBy).first();p0
ret=adjprc.groupby(by=GroupBy1).last().pct_change();
ret.index=pd.to_datetime(ret.index.astype(int).astype(str),format='%Y%m')

# 标准化indicators
indicators_all.index.names=['Trddt','Scode']
del_col=['CP', 'CTA', 'EY', 'ROIC', 'SC'] # 出现大量的nan，影响估计
indicators_all=indicators_all[indicators_all.columns.difference(del_col)]
indi_standardized=indicators_all.groupby(by=['Trddt']).apply(lambda x:(x-x.mean())/x.std())
ret[ret==0]=np.nan
ret1=ret.stack().reindex(indi_standardized.index)
ret1.index.names=['Trddt','Scode']
time_set=set(indi_standardized.index.get_level_values(0)) # 日期集合
indi_set=indi_standardized.columns # 指标集合

M_hat=pd.DataFrame(0.0,index=pd.MultiIndex.from_product([time_set,indi_set]),columns=indi_set).sort_index() # 用于储存每个时间的M_hat矩阵
slice_labels_matrix=pd.DataFrame(np.nan,index=time_set,columns=ret.columns).sort_index() # 用于储存股票在不同时间的分组／切片
slice_num=20
qcut_q=np.arange(slice_num+1)/slice_num
qcut_labels=range(1,slice_num+1) # labels总比q的长度小1
X_j_bar=pd.DataFrame(np.nan,index=pd.MultiIndex.from_product([time_set,qcut_labels]),columns=indi_set).sort_index() # 用于记录每个切片的characteristics的均值
for t in slice_labels_matrix.index:
    #对target进行slice
    slice_labels_matrix.loc[t]=pd.qcut(ret.loc[t],q=qcut_q,labels=qcut_labels)
    N = slice_labels_matrix.loc[t].count()
    for i in qcut_labels:
        X_j_bar.loc[(t,i)]=indi_standardized.loc[(t,slice(None))].loc[slice_labels_matrix.loc[t]==i].mean()#.fillna(0.0) # 存在不少的nan，删除指标？
        M_hat.loc[(t, slice(None)), :] +=np.kron(X_j_bar.loc[(t,i)],X_j_bar.loc[(t,i)]).reshape(indi_set.size,indi_set.size)\
                                         *indi_standardized.loc[(t,slice(None))].loc[slice_labels_matrix.loc[t]==i].shape[0]/N

for t in slice_labels_matrix.index:
    indi_standardized.loc[t].fillna(0.0).T@indi_standardized.loc[t].fillna(0)
    a,b=np.linalg.eigh(indi_standardized.loc[t].cov())
    tmp=np.zeros((indi_set.size,indi_set.size))


    for i in indi_standardized.loc[t].index:
        tmp+=np.kron(indi_standardized.loc[(t,i)].fillna(0.0),indi_standardized.loc[(t,i)].fillna(0.0)).reshape(indi_set.size,indi_set.size)
    tmp/=indi_standardized.loc[t].index.size;tmp
    eig_vals, _ = np.linalg.eigh(tmp);eig_vals#.cumsum() / eig_vals.sum()
    np.linalg.det(tmp)
M_hat.loc[t]
X_j_bar.loc[t].cov()

eig_vals,_=np.linalg.eigh(X_j_bar.loc[t].cov());eig_vals.cumsum()/eig_vals.sum()
eig_vals,_=np.linalg.eigh(M_hat.loc[t]);eig_vals.cumsum()/eig_vals.sum()