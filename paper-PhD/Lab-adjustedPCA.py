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

# step 1: X=a+b*Z
ret1=ret.stack().loc[indi_standardized.index]
ret1.index.names=['Trddt','Scode']
#indi_standardized.corrwith(ret1)
i_set=set(indi_standardized.index.get_level_values(0)) # 日期
indi_set=indi_standardized.columns
lambda_a=pd.DataFrame(np.nan,index=i_set,columns=indi_set).sort_index()
for i in lambda_a.index:
    lambda_a.loc[i]=indi_standardized.loc[i].corrwith(ret1.loc[i])#/ret1.loc[i].std()
lambda_a=lambda_a.astype(float)
#ii=lambda_a.index[6]
#indi_standardized.loc[ii,'SC'] # 估计的lambda_a存在missing value
#indicators_all.loc[ii,'SC']
selected=range(2)
lambda_a.iloc[:,selected].plot();plt.show()
lambda_a.iloc[:,selected].rolling(window=12,min_periods=1).mean().plot();plt.show()
lambda_a.iloc[:,selected].expanding().mean().plot();plt.show()
rho_lambda_a=lambda_a.corr().astype(float) # 为什么会出现大量的nan？？？---> 剔除inf后，rho中不再有nan，但是可能子样本中仍然会出现
rho_lambda_a.isnull().sum().sum()
eig_vals,eig_vecs=np.linalg.eigh(rho_lambda_a.astype(float))
eig_vals.cumsum()/eig_vals.sum()
f=lambda_a.fillna(0)@eig_vecs
f.iloc[:,[0,-1]].plot();plt.show()
f.iloc[:,-1].cumsum().plot();plt.show()
f_eig_vals=f/np.sqrt(eig_vals)
np.corrcoef(f.iloc[:,:5].T)

np.linalg.eigvalsh(rho_lambda_a.astype(float))
np.abs(lambda_a).mean().max()


rho_f=f.iloc[:,[0,1,2,3,-3,-2,-1]].rolling(window=12,min_periods=10).cov()
rho_f.loc[(slice(None),0),:].dropna().plot();plt.show()















lambda_a.isnull().sum()>0
col_nan=lambda_a.loc[:,lambda_a.isnull().sum()>0].columns;col_nan # lambda有缺失的变量

col_nornumll=rho_lambda_a.columns.difference(col_nan)
eig_vals, eig_vecs=np.linalg.eig(rho_lambda_a.loc[col_nornumll,col_nornumll])
eig_vals
eig_vals.cumsum()/eig_vals.sum()
