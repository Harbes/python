###改进后的代码
import pandas as pd
from pandas import DataFrame as DF
from pandas import Series as SS
import numpy as np
#非必要
from time import time
now0=time()
#参数设置
group_number=5
group_label=[i+1 for i in range(group_number)]  
percentile=np.linspace(0,1,group_number+1) 

#数据导入
data=pd.read_pickle('F:/data/xccdata/PV')[['adj_close','size_free','size_tot']]
data=pd.read_pickle('/Users/harbes/data/xccdata/PV')#[['trddt','stkcd','adj_close','size_free','size_tot']]
#data=pd.read_pickle('F:/data/xccdata/PV')#[['stkcd','trddt','adj_close','size_free','size_tot']]


#数据处理与删除重复条目
data['trddt']=pd.to_datetime(data['trddt'].astype(int).astype(str),format='%Y%m%d')
data=data.set_index(['trddt','stkcd'])
data.drop_duplicates(subset=None, keep='last',inplace=True)
#del data['trddt']
#del data['stkcd']
data.sort_index().to_pickle('/Users/harbes/data/xccdata/PV_datetime')

#计算回报率
data_adj_close=data['adj_close'].unstack()
data_rtn=data_adj_close/data_adj_close.shift(1)-1

#分组与结果展示
data_size_free=data['size_free'].unstack()               
data_group=DF([pd.qcut(data_size_free.ix[i],q=percentile,labels=group_label) for i in range(len(data_size_free))],index=data_size_free.index,columns=data_size_free.columns)                    
data_rtn_group=[[data_rtn.ix[i+1].ix[data_group.ix[i]==j+1].sum()/(data_group.ix[i]==j+1).sum()for j in range(group_number)]for i in range(len(data_rtn)-1)]
data_rtn_group_sum=DF((np.array(data_rtn_group)+1).cumprod(axis=0),index=data_rtn.index[1:],columns=list('12345'))
data_rtn_group_sum.plot()

print(time()-now0)#20秒左右




tmp=data[:10]

