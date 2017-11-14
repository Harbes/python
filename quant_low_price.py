import pandas as pd
import numpy as np
from pandas import DataFrame

# 整理数据
#data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_open','adj_close']]
data=pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')[['amount','opnprc','clsprc','adj_open','adj_close','size_tot','size_free']]
filter_=pd.read_pickle('/Users/harbes/data/xccdata/filter') # 4672606个有效数据点(原来有6140094个数据点)
opnprc=data['adj_open'].unstack()[filter_==1]
clsprc=data['adj_close'].unstack()[filter_==1]
key = lambda x: x.year * 100 + x.month
price0 = opnprc.groupby(key).first() # 买入价
price1 = clsprc.groupby(key).last() # 卖出价
price0.index = pd.to_datetime(price0.index.values.astype(str),format=('%Y%m'))
price1.index = pd.to_datetime(price1.index.values.astype(str),format=('%Y%m'))
rtn = (price1-price0)/price0 # 月度收益数据

price=data['clsprc'].unstack()[filter_==1]
indi=price.groupby(key).mean()
indi.index = pd.to_datetime(indi.index.values.astype(str),format=('%Y%m'))

group_num=10
percentile=np.linspace(0,1,group_num+1)
label_indi=[i+1 for i in range(group_num)] # 1表示流动性较好的组合，另一个极端组合是流动性较差的组合

# 分组，并计算不同组合收益
mark_indi=DataFrame([pd.qcut(indi.iloc[i],q=percentile,labels=label_indi) for i in range(len(indi))],index=indi.index)
rtn_group=DataFrame([[rtn.iloc[i][mark_indi.iloc[i-1]==k].mean() for k in label_indi] for i in range(1,len(rtn))],index=rtn.index[1:],columns=label_indi)

(rtn_group+1).cumprod().plot()
rtn_group.mean()
tmp=rtn_group[10]-rtn_group[1];tmp.mean()/tmp.std()*np.sqrt(len(tmp))




