import pandas as pd
import numpy as np
from pandas import DataFrame
from datetime import datetime
import matplotlib.pyplot as plt


# 整理数据
#data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_open','adj_close']]
data=pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')[['amount','opnprc','clsprc','adj_open','adj_close','size_tot','size_free']]
filter_=pd.read_pickle('/Users/harbes/data/xccdata/filter') # 4672606个有效数据点(原来有6140094个数据点)

key=lambda x:x.year*100+x.month

size=data['size_tot'].unstack()[filter_==1]
size_month=size.groupby(key).mean()
size_month.index = pd.to_datetime(size_month.index.values.astype(str),format=('%Y%m'))


num_by_=10
label_=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)

# 按同期size分组，或者按前一期size进行分组
mark_=DataFrame([pd.qcut(size_month.loc[m],q=percentile,labels=label_) for m in size_month.index],index=size_month.index,columns=size_month.columns)
mark_=mark_.asfreq(freq='D',method='ffill').loc[size.index]

# 不同size组合收益，Equally-weighted
opnprc=data['adj_open'].unstack()#[filter_==1]
clsprc=data['adj_close'].unstack()#[filter_==1]
rtn=(clsprc-opnprc)/opnprc
rtn_port=DataFrame([[rtn.loc[m][mark_.loc[m]==l_].mean() for l_ in label_] for m in mark_.index],index=mark_.index,columns=label_)

# market portfolio
#rtn_port['M']=(size*rtn).mean(axis=1)/(size[~np.isnan(size*rtn)]).mean(axis=1) # value-weighted
rtn_port['M']=rtn.mean(axis=1) # equally-weighted

# 设置第二index
rtn_port['index']=(rtn_port.index.year-2005)*12+rtn_port.index.month-1
max_month=rtn_port['index'][-1]
rtn_port['trddt']=rtn_port.index
rtn_port=rtn_port.set_index(['trddt','index'])

# 计算 beta (前12个月数据滚动计算，每月计算一次)
#rtn_port.loc[(slice(None),range(0,12)),slice(None)]#.cov().loc['M',label_]
beta_=DataFrame([rtn_port.loc[(slice(None),range(np.maximum(0,i-12),i)),slice(None)].cov().loc['M',label_]/rtn_port.loc[(slice(None),range(np.maximum(0,i-12),i)),slice(None)].cov().loc['M','M'] for i in range(1,len(size_month))],index=size_month.index[1:],columns=label_)






