import pandas as pd
import numpy as np
from pandas import DataFrame

# 整理数据
data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_open','adj_close']]
#open=data['adj_open'].unstack()
clsprc=data['adj_close'].unstack()
key = lambda x: x.year * 100 + x.month
price0 = clsprc.groupby(key).first() # 买入价
price1 = clsprc.groupby(key).last() # 卖出价
price0.index = pd.to_datetime(price0.index.values.astype(str),format=('%Y%m'))
price1.index = pd.to_datetime(price1.index.values.astype(str),format=('%Y%m'))

J=3
K=3
percentile=np.linspace(0,1,11)
label_momentum=[i+1 for i in range(10)] # 1表示过去收益较差的组合，10表示过去收益较好的组合
mark_momentum=DataFrame([pd.qcut((price1.iloc[i-1]-price0.iloc[i-J])/price0.iloc[i-J],q=percentile,labels=label_momentum) for i in range(J,len(price0))],index=price0.index[J:])
rtn=((price1-price0)/price0).iloc[J:]
momentum=DataFrame([[rtn.iloc[i][mark_momentum.iloc[i]==k].mean() for k in label_momentum] for i in range(K-1,len(rtn))],index=rtn.index[K-1:],columns=label_momentum)
if K>1:
    for lag in range(1,K):
        tmp=DataFrame([[rtn.iloc[i][mark_momentum.iloc[i-lag]==k].mean() for k in label_momentum] for i in range(K-1,len(rtn))],index=rtn.index[K-1:],columns=label_momentum)
        momentum += tmp
momentum /= K

f_momen=momentum[10]-momentum[1]
f_momen.mean()/f_momen.std()*np.sqrt(len(f_momen))


