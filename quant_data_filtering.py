import pandas as pd
import numpy as np
from pandas import DataFrame
data_path=''

# 涨跌停
data=pd.read_pickle(data_path+'PV_datetime')[['opnprc','high','low']]
opnprc=data['opnprc'].unstack()
high=data['high'].unstack()
low=data['low'].unstack()
up=(high-opnprc)/opnprc>=0.1
down=(low-opnprc)/opnprc<=-0.1
limit_move=np.logical_or(up,down)
limit_move.sum().sum() # 剔除113640个数据点
limit_move.to_pickle(data_path+'limit_move')


# 非ST股票，且非新股(约三个月)
ST=pd.read_pickle(data_path+'ST') # N表示非ST，且剔除IPO三个月交易数据
ST.index=ST.index.astype(int).astype(str).to_datetime('%Y%m%d')
ST.to_pickle(data_path+'ST_datetime')


# 停牌NT，以及停牌
NT=pd.read_pickle(data_path+'NT') # 1表示停牌
NT.index=NT.index.astype(int).astype(str).to_datetime('%Y%m%d')
NT.to_pickle(data_path+'NT_datetime')

def filter_NT_extension(NT,up_limit1):
    T,N = NT.shape
    NT_filter=np.zeros((T,N))
    for i in range(N):
        for j in range(up_limit1,T):
            if NT.values[j,i]==0 and NT.values[j-up_limit1:j,i].mean()>0.5:
                NT_filter[j,i]=1.0
    return DataFrame(NT_filter,index=NT.index,columns=NT.columns)
NT_filter=filter_NT_extension(NT,120) # 剔除了119464个数据点
NT_filter.sum().sum()
NT_filter.to_pickle(data_path+'NT_filter_120')


# 成交量较小还是成交金额？？？
data=pd.read_pickle(data_path+'PV_datetime')['amount'].unstack()
data[data==0]=np.nan
amount_filter=DataFrame([pd.qcut(data.iloc[i],q=[0.0,0.01,1.0],labels=[0,1]) for i in range(len(data))])
amount_filter.to_pickle(data_path+'amount_filter')


# 归总筛选
filter=amount_filter[amount_filter==1][ST=='N'][~limit_move][NT==0][NT_filter==0]
filter.to_pickle(data_path+'filter') # 4672606个有效数据点(原来有6140094个数据点)



