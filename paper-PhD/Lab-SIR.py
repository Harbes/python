import pandas as pd
import numpy as np
from pandas.tseries.offsets import DateOffset
import scipy
import matplotlib.pyplot as plt
DPath='E:/data/CNRDS/' #'/Users/harbes/PycharmProjects/data/CNRDS/' #
indicators_all=pd.read_pickle(DPath+'indicators_all').replace([np.inf, -np.inf], np.nan)
PV=pd.read_pickle(DPath+'PV');PV.tail(10)
tmp=PV[['Scode','Trddt','Adclsprc']].iloc[1:].set_index(['Trddt','Scode'])
tmp=tmp.astype(float)
adjprc=tmp.loc[~tmp.index.duplicated(keep='last')].unstack() # 提出duplicates
adjprc.index=pd.to_datetime(adjprc.index,format='%Y-%m-%d')# 调整index格式
adjprc.columns=adjprc.columns.get_level_values(1)# 调整columns格式
GroupBy1=lambda x:x.year*100.0+x.month
#p0=adjprc.groupby(by=GroupBy).first();p0
ret=adjprc.groupby(by=GroupBy1).last().pct_change();ret
ret.index=pd.to_datetime(ret.index.astype(int).astype(str),format='%Y%m')

# 标准化indicators
indicators_all.index.names=['Trddt','Scode']
indi_standardized=indicators_all.groupby(by=['Trddt']).apply(lambda x:(x-x.mean())/x.std())
ret[ret==0]=np.nan
ret1=ret.stack().loc[indi_standardized.index]
ret1.index.names=['Trddt','Scode']
time_set=set(indi_standardized.index.get_level_values(0)) # 日期
indi_set=indi_standardized.columns

M_hat=pd.DataFrame(0.0,index=pd.MultiIndex.from_product([time_set,indi_set]),columns=indi_set).sort_index()
slice_labels_matrix=pd.DataFrame(np.nan,index=time_set,columns=ret.columns).sort_index()
slice_num=20
qcut_q=np.arange(slice_num+1)/slice_num
qcut_labels=range(1,slice_num+1) # labels总比q的长度小1
X_j_bar=pd.DataFrame(np.nan,index=pd.MultiIndex.from_product([time_set,qcut_labels]),columns=indi_set).sort_index()
for t in time_set:
    #对target进行slice
    slice_labels_matrix.loc[t]=pd.qcut(ret.loc[t],q=qcut_q,labels=qcut_labels)
    N = slice_labels_matrix.loc[t].count()
    for i in qcut_labels:
        X_j_bar.loc[(t,i)]=indi_standardized.loc[(t,slice(None))].loc[slice_labels_matrix.loc[t]==i].mean().fillna(0.0)
        M_hat.loc[(t, slice(None)), :] +=np.kron(X_j_bar.loc[(t,i)],X_j_bar.loc[(t,i)]).reshape(67,67)\
                                         *indi_standardized.loc[(t,slice(None))].loc[slice_labels_matrix.loc[t]==i].shape[0]/N

for t in time_set:

eig_vals,_=np.linalg.eig(M_hat.loc[(t, slice(None)), :]);eig_vals