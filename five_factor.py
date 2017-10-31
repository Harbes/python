import numpy as np
import pandas as pd
from datetime import datetime

################# 数据准备 #####################
## 营业利润 Profit
data=pd.read_pickle('/Users/harbes/data/xccdata/IS')[['stkcd','fin_year','oper_profit']].set_index(['fin_year','stkcd']).sort_index().unstack()
data.index=data.index.astype(str).to_datetime('%Y%m%d')
data.index.name='trddt'
#data.to_pickle('./five_factor/oper_profit')
from pandas.tseries.offsets import MonthBegin
data=data.resample('M').first().ffill()
data.index=data.index+MonthBegin()
#data.to_pickle('./five_factor/oper_profit_monthly')

## 权益的账面价值 book
data=pd.read_pickle('/Users/harbes/data/xccdata/BS')[['fin_year','stkcd','tot_assets','tot_liab']].set_index(['fin_year','stkcd']).sort_index()
data['book']=data['tot_assets']-data['tot_liab']
data=data['book'].unstack()
data.index=data.index.astype(str).to_datetime('%Y%m%d')
data.index.name='trddt'
#data.to_pickle('./five_factor/book')
data=data.resample('M').first().ffill()
data.index=data.index+MonthBegin()
#data.to_pickle('./five_factor/book_monthly')

## 市值 M 以及价格数据 adj_price
data=pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')[['size_free','size_tot']]
data=data['size_free'].unstack()
key=lambda x:x.year*100+x.month
data=data.groupby(key).last()
data.index=pd.to_datetime(data.index.values.astype(str),format='%Y%m')
#data.to_pickle('./five_factor/size_monthly')

## 投资 inv
data=pd.read_pickle('/Users/harbes/data/xccdata/BS')[['fin_year','stkcd','tot_assets','tot_liab']].set_index(['fin_year','stkcd']).sort_index()
data['book']=data['tot_assets']-data['tot_liab']
data=data['tot_assets'].unstack()
data.index=data.index.astype(str).to_datetime('%Y%m%d')
data.index.name='trddt'
data=data/data.shift(1)-1
data=data.resample('M').first().ffill()
data.index=data.index+MonthBegin()
data.to_pickle('./five_factor/inv_monthly')

## 价格数据 adj_price
data=pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')['adj_close'].unstack()
key=lambda x:x.year*100+x.month
data=data.groupby(key).last()
data.index=pd.to_datetime(data.index.values.astype(str),format='%Y%m')
#data.to_pickle('./five_factor/price_monthly')

################# 分组标记 #####################
num_by_size=2
label_size=[i+1 for i in range(num_by_size)] # 1表示小市值，2表示大市值
percentile=np.linspace(0,1,num_by_size+1)
size=pd.read_pickle('./five_factor/size_monthly')
mark_size=pd.DataFrame([pd.qcut(size.iloc[i],q=percentile,labels=label_size) for i in range(len(size))],index=size.index,columns=size.columns)
mark_size.to_pickle('./five_factor/mark_size')

label_BM=np.arange(3,0,-1) # 1表示小市值、高价值股票，3表示小市值、低价值成长股;4表示大市值、高价值股票，6表示大市值、低价值股票
label_BM=np.array([3,2,1])
percentile=np.array([0.0,0.3,0.7,1.0])
t0=datetime.strptime('2005-01-01','%Y-%m-%d')
book=pd.read_pickle('./five_factor/book_monthly').loc[t0:]
size=pd.read_pickle('./five_factor/size_monthly')
book=book[size.columns]
BM=book/size/10000.0
mark_BM=pd.DataFrame([pd.qcut(BM[mark_size==1].iloc[i],q=percentile,labels=label_BM) for i in range(len(BM))],index=BM.index)
mark_B_BM=pd.DataFrame([pd.qcut(BM[mark_size==2].iloc[i],q=percentile,labels=label_BM+3) for i in range(len(BM))],index=BM.index)
mark_BM=mark_BM.combine_first(mark_B_BM)
mark_BM.to_pickle('./five_factor/mark_BM')

