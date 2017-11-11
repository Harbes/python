import pandas as pd
import numpy as np
#from pandas import DataFrame
from datetime import datetime

# 整理数据
#data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_open','adj_close']]
data=pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')[['adj_open','adj_close']]
#open=data['adj_open'].unstack()
clsprc=data['adj_close'].unstack()
key = lambda x: x.year * 100 + x.month
price0 = clsprc.groupby(key).first() # 买入价
price1 = clsprc.groupby(key).last() # 卖出价
price0.index = pd.to_datetime(price0.index.values.astype(str),format=('%Y%m'))
price1.index = pd.to_datetime(price1.index.values.astype(str),format=('%Y%m'))

t_end=datetime.strptime('2016-12-01','%Y-%m-%d')
rtn_gross=(price1/price0).mean(axis=1).loc[:t_end]


key2=lambda x:x.month
key1=lambda x:x.year
rtn_by_month=rtn_gross.groupby([key1,key2]).last().unstack()
rtn_by_month.cumprod().plot()
(rtn_by_month-1).mean()/(rtn_by_month-1).std()*np.sqrt(len(rtn_by_month))
