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
data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_close','size_free','size_tot']]
data=pd.read_pickle('/Users/harbes/data/xccdata/PV')#[['trddt','stkcd','adj_close','size_free','size_tot']]
#data=pd.read_pickle('F:/data/xccdata/PV')#[['stkcd','trddt','adj_close','size_free','size_tot']]


#数据处理与删除重复条目
#data['trddt']=pd.to_datetime(data['trddt'].astype(int).astype(str),format='%Y%m%d')
#data=data.set_index(['trddt','stkcd'])
#data.drop_duplicates(subset=None, keep='last',inplace=True)
#del data['trddt']
#del data['stkcd']
#data.sort_index().to_pickle('F:/data/xccdata/PV_datetime')

#计算回报率
data_adj_close=data['adj_close'].unstack()
data_rtn=(data_adj_close/data_adj_close.shift(1)-1)*100

#分组与结果展示
data_size_free=data['size_free'].unstack()               
data_group=DF([pd.qcut(data_size_free.iloc[i],q=percentile,labels=group_label) for i in range(len(data_size_free))],index=data_size_free.index,columns=data_size_free.columns)
data_rtn_group=[[data_rtn.iloc[i+1].loc[data_group.iloc[i]==j+1].sum()/(data_group.iloc[i]==j+1).sum()for j in range(group_number)]for i in range(len(data_rtn)-1)]
data_rtn_group_sum=DF((np.array(data_rtn_group)+1).cumprod(axis=0),index=data_rtn.index[1:],columns=list('12345'))
data_rtn_group_sum.plot()

print(time()-now0)#20秒左右




### 用于计算 半方差、方差等指标
key=lambda x:x.year*100+x.month
tmp=data_rtn[:100].groupby('key')
tmp.std()

def pos_semi_var(arr):
    return np.mean((arr-arr.mean())**2* (arr-arr.mean() > 0))

def neg_semi_var(arr):
    return np.mean((arr - arr.mean()) ** 2 * (arr-arr.mean()< 0))
def var(arr):
    return np.mean((arr-arr.mean())**2)

data_rtn[:100].groupby(key).apply(pos_semi_var).iloc[0,0]
data_rtn[:100].groupby(key).apply(neg_semi_var).iloc[0,0]
tmp=data_rtn[:500].groupby(key).apply(var)
tmp.index=pd.to_datetime(tmp.index.values.astype(str),format='%Y%m')




### B/M 中的book 数据

data=pd.read_pickle('F:/data/xccdata/BS')[['fin_year','stkcd','tot_assets','tot_liab']]
data['fin_year']=pd.to_datetime(data['fin_year'].astype(int).astype(str),format='%Y%m%d')
data['book']=data['tot_assets']-data['tot_liab']
from pandas.tseries.offsets import Day,MonthBegin
data['trddt']=data['fin_year']+MonthBegin()
data=data.set_index(['trddt','stkcd']).sort_index()
data=data['book'].unstack()
tmp=data.resample('M').first().ffill()
tmp.index=tmp.index-MonthBegin()
tmp.to_pickle('F:/data/xccdata/essay/book')


### market return数据
data=pd.read_pickle('F:/data/xccdata/index_ret')
data['trddt']=pd.to_datetime(data['trddt'].astype(int).astype(str),format='%Y%m%d')
data=data.set_index('trddt')
data=data.loc[data['index_code']=='399300.SZ']
data['pctchange'].to_pickle('F:/data/xccdata/essay/index_hs300_daily')

# 月初or月末
key=lambda x:x.year*100+x.month
tmp=data['clsprc'].groupby(key)
tmp=tmp.last()
tmp.index=pd.to_datetime(tmp.index.values.astype(str),format='%Y%m')
tmp.to_pickle('F:/data/xccdata/essay/index_hs300_monthend')



