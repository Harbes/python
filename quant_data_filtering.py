import pandas as pd
import numpy as np
from pandas import DataFrame
data_path = 'E:/data/NewData/'  # '/Users/harbes/data/NewData/' #
data_path0 = 'E:/data/xccdata/'  # '/Users/harbes/data/xccdata/' #
# 涨跌停
data=pd.read_pickle(data_path+'PV_datetime')[['opnprc','high','low']]
opnprc=data['opnprc'].unstack()
high=data['high'].unstack()
low=data['low'].unstack()
up=(high-opnprc)/opnprc>=0.1
down=(low-opnprc)/opnprc<=-0.1
limit_move=np.logical_or(up,down)
limit_move.sum().sum() # 剔除113640个数据点
(~limit_move).to_pickle(data_path+'limit_move')


# 非ST股票，且非新股(约三个月)
ST=pd.read_pickle(data_path+'ST');ST.head() # N表示非ST，且剔除IPO三个月交易数据【财神已标记】
ST.index=pd.to_datetime(ST.index.astype(int).astype(str),format='%Y%m%d');ST.head()
(ST=='N').to_pickle(data_path+'ST_datetime')


# 停牌NT，以及停牌
NT=pd.read_pickle(data_path+'NT');NT.head() # 1表示停牌
NT.index=pd.to_datetime(NT.index.astype(int).astype(str),format='%Y%m%d');NT.head()
#(NT==0.0).to_pickle(data_path+'NT_datetime')

NT50=NT.rolling(50).mean()>0.99 # 标记超过停牌不小于50天的点
NT50 +=0;NT50.sum().sum()
tmp=NT50.copy();tmp.head()
for i in range(100):
    NT50[i+1:] +=tmp.shift(i+1)[i+1:] #停牌超过50天后的n个交易日都作废
#((NT+NT50)==0.0).head()
((NT+NT50)==0.0).to_pickle(data_path+'NT&NT50_datetime')


# 成交量较小还是成交金额？？？
data=pd.read_pickle(data_path+'PV_datetime')[['amount','clsprc']]
amount=data['amount'].unstack()
clsprc=data['clsprc'].unstack()
amount /=clsprc;amount.head()
amount[amount==0]=np.nan
amount_filter=DataFrame([pd.qcut(amount.iloc[i],q=[0.0,0.01,1.0],labels=[0,1]) for i in range(len(amount))])
amount_filter.to_pickle(data_path+'FilterSmallVolume')


# 归总筛选
filter_=amount_filter[amount_filter>0][ST=='N'][~limit_move][NT==0.0];filter_.head()
(filter_==1.0).to_pickle(data_path+'filtered_data') # 4672606个有效数据点(原来有6140094个数据点)
(filter_==1.0).sum(axis=1).sum()




