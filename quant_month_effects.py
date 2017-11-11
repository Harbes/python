import pandas as pd
import numpy as np
#from pandas import DataFrame
from datetime import datetime

# 整理数据
#data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_open','adj_close']]
data=pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')[['adj_open','adj_close']]
open=data['adj_open'].unstack()
clsprc=data['adj_close'].unstack()
key = lambda x: x.year * 100 + x.month

index_s=pd.Series(open.index)
index_s.index=index_s.values
index_s.groupby(key).describe()['count'] # 不同月份交易日数
tmp = open.groupby(key).median()
price_open_1 = open.groupby(key).first() #  第一周
price_open_2 = open.groupby(key).nth(5) #  第二周
price_open_3 = open.groupby(key).nth(11) #  第三周
price_open_4 = open.groupby(key).nth(17) #  月末

price_close_1 = clsprc.groupby(key).nth(4) #  第一周
price_close_2 = clsprc.groupby(key).nth(10) #  第二周
price_close_3 = clsprc.groupby(key).nth(16) #  第三周
price_close_4 = clsprc.groupby(key).last() #  月末

price_open_1.index = pd.to_datetime(price_open_1.index.values.astype(str),format=('%Y%m'))
price_open_2.index = pd.to_datetime(price_open_2.index.values.astype(str),format=('%Y%m'))
price_open_3.index = pd.to_datetime(price_open_3.index.values.astype(str),format=('%Y%m'))
price_open_4.index = pd.to_datetime(price_open_4.index.values.astype(str),format=('%Y%m'))
price_close_1.index = pd.to_datetime(price_close_1.index.values.astype(str),format=('%Y%m'))
price_close_2.index = pd.to_datetime(price_close_2.index.values.astype(str),format=('%Y%m'))
price_close_3.index = pd.to_datetime(price_close_3.index.values.astype(str),format=('%Y%m'))
price_close_4.index = pd.to_datetime(price_close_4.index.values.astype(str),format=('%Y%m'))

t_end=datetime.strptime('2016-12-01','%Y-%m-%d')
rtn_gross=(price1/price0).mean(axis=1).loc[:t_end]


key2=lambda x:x.month
key1=lambda x:x.year
rtn_by_month=rtn_gross.groupby([key1,key2]).last().unstack()
rtn_by_month.cumprod().plot()
(rtn_by_month-1).mean()/(rtn_by_month-1).std()*np.sqrt(len(rtn_by_month))
