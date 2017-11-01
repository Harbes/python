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

## 市值 M
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

## 价格数据 adj_price => rtn
data=pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')['adj_close'].unstack()
key=lambda x:x.year*100+x.month
tmp=data.groupby(key).last()
tmp.index=pd.to_datetime(tmp.index.values.astype(str),format='%Y%m')
t_add=datetime.strptime('2004-12-01','%Y-%m-%d')
tmp.loc[t_add]=data.groupby(key).first().iloc[0]
tmp=tmp.sort_index()
rtn=tmp.pct_change()*100.0
#rtn.iloc[1:].to_pickle('./five_factor/rtn_monthly')



################# 分组标记 #####################
size=pd.read_pickle('./five_factor/size_monthly')
oper=pd.read_pickle('./five_factor/oper_profit_monthly')
stock_selected=oper.columns & size.columns
### 均匀分组 or 分开分组 ###
num_by_size=2
label_size=[i+1 for i in range(num_by_size)] # 1表示小市值，2表示大市值
percentile=np.linspace(0,1,num_by_size+1)
size=pd.read_pickle('./five_factor/size_monthly')[stock_selected]
mark_size=pd.DataFrame([pd.qcut(size.iloc[i],q=percentile,labels=label_size) for i in range(len(size))],index=size.index,columns=size.columns)
#mark_size.to_pickle('./five_factor/mark_size')

label_BM=np.arange(3,0,-1) # 1表示小市值、高价值股票，3表示小市值、低价值成长股;4表示大市值、高价值股票，6表示大市值、低价值股票
percentile=np.array([0.0,0.3,0.7,1.0])
t0=datetime.strptime('2005-01-01','%Y-%m-%d')
book=pd.read_pickle('./five_factor/book_monthly').loc[t0:][stock_selected]
size=pd.read_pickle('./five_factor/size_monthly')[stock_selected]
BM=book/size/10000.0
mark_BM=pd.DataFrame([pd.qcut(BM[mark_size==1].iloc[i],q=percentile,labels=label_BM) for i in range(len(BM))],index=BM.index)
mark_B_BM=pd.DataFrame([pd.qcut(BM[mark_size==2].iloc[i],q=percentile,labels=label_BM+3) for i in range(len(BM))],index=BM.index)
mark_BM=mark_BM.combine_first(mark_B_BM)
#mark_BM.to_pickle('./five_factor/mark_BM')

label_OP=np.arange(3,0,-1) # 1表示小市值、盈利稳健股票，3表示小市值、weak成长股;4表示大市值、robust股票，6表示大市值、weak股票
percentile=np.array([0.0,0.3,0.7,1.0])
t0=datetime.strptime('2005-01-01','%Y-%m-%d')
book=pd.read_pickle('./five_factor/book_monthly').loc[t0:][stock_selected]
oper=pd.read_pickle('./five_factor/oper_profit_monthly').loc[t0:][stock_selected]
OP=oper/book*10.0
mark_OP=pd.DataFrame([pd.qcut(OP[mark_size==1].iloc[i],q=percentile,labels=label_OP) for i in range(len(OP))],index=OP.index)
mark_B_OP=pd.DataFrame([pd.qcut(OP[mark_size==2].iloc[i],q=percentile,labels=label_OP+3) for i in range(len(OP))],index=OP.index)
mark_OP=mark_OP.combine_first(mark_B_OP)
#mark_OP.to_pickle('./five_factor/mark_OP')

label_INV=np.arange(1,4) # 1表示小市值、conservative股票，3表示小市值、aggressive成长股;4表示大市值、conservative股票，6表示大市值、aggressive股票
percentile=np.array([0.0,0.3,0.7,1.0])
t0=datetime.strptime('2005-01-01','%Y-%m-%d')
INV=pd.read_pickle('./five_factor/inv_monthly').loc[t0:][stock_selected]*100.0
mark_INV=pd.DataFrame([pd.qcut(INV[mark_size==1].iloc[i],q=percentile,labels=label_INV) for i in range(len(INV))],index=INV.index)
mark_B_INV=pd.DataFrame([pd.qcut(INV[mark_size==2].iloc[i],q=percentile,labels=label_INV+3) for i in range(len(INV))],index=INV.index)
mark_INV=mark_INV.combine_first(mark_B_INV)
#mark_INV.to_pickle('./five_factor/mark_INV')


### 独立分组 ###
num_by_size=2
label_size=[i+1 for i in range(num_by_size)] # 1表示小市值，2表示大市值
percentile=np.linspace(0,1,num_by_size+1)
size=pd.read_pickle('./five_factor/size_monthly')[stock_selected]
mark_size=pd.DataFrame([pd.qcut(size.iloc[i],q=percentile,labels=label_size) for i in range(len(size))],index=size.index,columns=size.columns)
#mark_size.to_pickle('./five_factor/mark_size')

label_BM=np.arange(3,0,-1) # 1表示高价值股票，3表示低价值成长股
percentile=np.array([0.0,0.3,0.7,1.0])
t0=datetime.strptime('2005-01-01','%Y-%m-%d')
book=pd.read_pickle('./five_factor/book_monthly').loc[t0:][stock_selected]
size=pd.read_pickle('./five_factor/size_monthly')[stock_selected]
BM=book/size/10000.0
mark_BM=pd.DataFrame([pd.qcut(BM.iloc[i],q=percentile,labels=label_BM) for i in range(len(BM))],index=BM.index)
#mark_BM.to_pickle('./five_factor/mark_BM')

label_OP=np.arange(3,0,-1) # 1表示盈利robust股票，3表示盈利weak成长股
percentile=np.array([0.0,0.3,0.7,1.0])
t0=datetime.strptime('2005-01-01','%Y-%m-%d')
book=pd.read_pickle('./five_factor/book_monthly').loc[t0:][stock_selected]
oper=pd.read_pickle('./five_factor/oper_profit_monthly').loc[t0:][stock_selected]
OP=oper/book*10.0
mark_OP=pd.DataFrame([pd.qcut(OP.iloc[i],q=percentile,labels=label_OP) for i in range(len(OP))],index=OP.index)
#mark_OP.to_pickle('./five_factor/mark_OP')

label_INV=np.arange(1,4) # 1表示conservative股票，3表示aggressive成长股
percentile=np.array([0.0,0.3,0.7,1.0])
t0=datetime.strptime('2005-01-01','%Y-%m-%d')
INV=pd.read_pickle('./five_factor/inv_monthly').loc[t0:][stock_selected]*100.0
mark_INV=pd.DataFrame([pd.qcut(INV.iloc[i],q=percentile,labels=label_INV) for i in range(len(INV))],index=INV.index)
#mark_INV.to_pickle('./five_factor/mark_INV')




######### 计算多空组合收益（市值加权）#########
### 对应"独立分组" ###
size=pd.read_pickle('./five_factor/size_monthly')[stock_selected]/100000.0
rtn=pd.read_pickle('./five_factor/rtn_monthly')[stock_selected]
#rtn=price.pct_change()*100.0

### size-B/M ###
S_BM=pd.DataFrame([[(rtn.iloc[i+1]*size.iloc[i])[mark_size.iloc[i]==1][mark_BM.iloc[i]==j+1].sum()
                    /(size.iloc[i]*(~np.isnan(rtn.iloc[i+1]*size.iloc[i])))[mark_size.iloc[i]==1][mark_BM.iloc[i]==j+1].sum() for j in range(3)]
                   for i in range(len(rtn)-1)],index=rtn.index[1:],columns=['H','N','L'])
B_BM=pd.DataFrame([[(rtn.iloc[i+1]*size.iloc[i])[mark_size.iloc[i]==2][mark_BM.iloc[i]==j+1].sum()
                    /(size.iloc[i]*(~np.isnan(rtn.iloc[i+1]*size.iloc[i])))[mark_size.iloc[i]==2][mark_BM.iloc[i]==j+1].sum() for j in range(3)]
                   for i in range(len(rtn)-1)],index=rtn.index[1:],columns=['H','N','L'])

S_OP=pd.DataFrame([[(rtn.iloc[i+1]*size.iloc[i])[mark_size.iloc[i]==1][mark_OP.iloc[i]==j+1].sum()
                    /(size.iloc[i]*(~np.isnan(rtn.iloc[i+1]*size.iloc[i])))[mark_size.iloc[i]==1][mark_OP.iloc[i]==j+1].sum() for j in range(3)]
                   for i in range(len(rtn)-1)],index=rtn.index[1:],columns=['R','N','W'])
B_OP=pd.DataFrame([[(rtn.iloc[i+1]*size.iloc[i])[mark_size.iloc[i]==2][mark_OP.iloc[i]==j+1].sum()
                    /(size.iloc[i]*(~np.isnan(rtn.iloc[i+1]*size.iloc[i])))[mark_size.iloc[i]==2][mark_OP.iloc[i]==j+1].sum() for j in range(3)]
                   for i in range(len(rtn)-1)],index=rtn.index[1:],columns=['R','N','W'])

S_INV=pd.DataFrame([[(rtn.iloc[i+1]*size.iloc[i])[mark_size.iloc[i]==1][mark_INV.iloc[i]==j+1].sum()
                    /(size.iloc[i]*(~np.isnan(rtn.iloc[i+1]*size.iloc[i])))[mark_size.iloc[i]==1][mark_INV.iloc[i]==j+1].sum() for j in range(3)]
                   for i in range(len(rtn)-1)],index=rtn.index[1:],columns=['C','N','A'])
B_INV=pd.DataFrame([[(rtn.iloc[i+1]*size.iloc[i])[mark_size.iloc[i]==2][mark_INV.iloc[i]==j+1].sum()
                    /(size.iloc[i]*(~np.isnan(rtn.iloc[i+1]*size.iloc[i])))[mark_size.iloc[i]==2][mark_INV.iloc[i]==j+1].sum() for j in range(3)]
                   for i in range(len(rtn)-1)],index=rtn.index[1:],columns=['C','N','A'])

tmp=(S_INV.mean(axis=1)-B_INV.mean(axis=1)+S_OP.mean(axis=1)-B_OP.mean(axis=1)+S_BM.mean(axis=1)-B_BM.mean(axis=1))/3
tmp.mean()/tmp.std()*np.sqrt(len(tmp))














