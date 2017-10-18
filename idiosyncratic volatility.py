###改进后的代码
import pandas as pd
from pandas import DataFrame as DF
from pandas import Series as SS
import numpy as np



from datetime import datetime
t0=datetime.strptime('2005-02-01','%Y-%m-%d')
t1=datetime.strptime('2017-09-29','%Y-%m-%d')
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

# 计算月度回报
price0=pd.read_pickle('F:/data/xccdata/essay/index_hs300_monthstart')
price=pd.read_pickle('F:/data/xccdata/essay/index_hs300_monthend')
rtn=(price/price.shift(1)-1)*100
rtn[0]=(price[0]/price0[0]-1)*100
rtn.to_pickle('F:/data/xccdata/essay/index_hs300_monthly')



### 计算mimicking portfolio

data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_close','size_free','size_tot']]
#data=pd.read_pickle('/Users/harbes/data/xccdata/PV')[['trddt','stkcd','adj_close','size_free','size_tot']]


#计算股票月度回报率和size月度数据
data_adj_close=data['adj_close'].unstack()
key=lambda x:x.year*100+x.month
tmp=data_adj_close.groupby(key)
price1=tmp.last();price1.index=pd.to_datetime(price1.index.values.astype(str),format='%Y%m')
price0=tmp.first();price0.index=pd.to_datetime(price0.index.values.astype(str),format='%Y%m')
price0.to_pickle('F:/data/xccdata/essay/stocks_clsprc_monthstart')
price1.to_pickle('F:/data/xccdata/essay/stocks_clsprc_monthend')
rtn=(price1/price1.shift(1)-1)*100
rtn.iloc[0]=(price1.iloc[0]/price0.iloc[0]-1)*100
rtn.to_pickle('F:/data/xccdata/essay/stocks_rtn_monthly')



size=data['size_tot'].unstack()
key=lambda x:x.year*100+x.month
tmp=size.groupby(key)
size_=tmp.last()
size_.index=pd.to_datetime(size_.index.values.astype(str),format='%Y%m')
size_.to_pickle('F:/data/xccdata/essay/stocks_size_tot_monthend')





# 计算 mimicking portfolio —— SMB
rtn=pd.read_pickle('F:/data/xccdata/essay/stocks_rtn_monthly')
size=pd.read_pickle('F:/data/xccdata/essay/stocks_size_tot_monthend')

group_number=5
group_label=[i+1 for i in range(group_number)]
percentile=np.linspace(0,1,group_number+1)

data_group=DF([pd.qcut(size.iloc[i],q=percentile,labels=group_label) for i in range(len(size))],index=size.index,columns=size.columns)
data_rtn_group=[[rtn.iloc[i+1].loc[data_group.iloc[i]==j+1].sum()/(data_group.iloc[i]==j+1).sum()for j in range(group_number)]for i in range(len(rtn)-1)]
data_rtn_group=DF(data_rtn_group,index=rtn.index[1:],columns=list('12345'))
mimick=data_rtn_group['1']-data_rtn_group['5']
mimick.to_pickle('F:/data/xccdata/essay/SMB_tot')
#data_rtn_group_sum=DF((np.array(data_rtn_group)+1).cumprod(axis=0),index=rtn.index[1:],columns=list('12345'))

# 计算 mimicking portfolio —— HML
rtn=pd.read_pickle('F:/data/xccdata/essay/stocks_rtn_monthly')
size=pd.read_pickle('F:/data/xccdata/essay/stocks_size_free_monthend') # 按总市值做出来的价值因子是编辑显著的
book=pd.read_pickle('F:/data/xccdata/essay/book')/10000
BM=book/size
t0=rtn.index[0]
BM=BM.loc[t0:]

data_group=DF([pd.qcut(BM.iloc[i],q=percentile,labels=group_label) for i in range(len(BM))],index=BM.index,columns=BM.columns)
data_rtn_group=[[rtn.iloc[i+1].loc[data_group.iloc[i]==j+1].sum()/(data_group.iloc[i]==j+1).sum()for j in range(group_number)]for i in range(len(rtn)-1)]
data_rtn_group=DF(data_rtn_group,index=rtn.index[1:],columns=list('12345'))
mimick=data_rtn_group['5']-data_rtn_group['1']
mimick.mean()/mimick.std()*np.sqrt(len(mimick))
mimick.to_pickle('F:/data/xccdata/essay/HML_free')






# 计算 SMB 和 HML 的日度数据
data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_close','size_free','size_tot']]
data_adj_close=data['adj_close'].unstack()
rtn=(data_adj_close/data_adj_close.shift(1)-1)*100
size=data['size_free'].unstack()

group_number=5
group_label=[i+1 for i in range(group_number)]
percentile=np.linspace(0,1,group_number+1)

data_group=DF([pd.qcut(size.iloc[i],q=percentile,labels=group_label) for i in range(len(size))],index=size.index,columns=size.columns)
data_rtn_group=[[rtn.iloc[i+1].loc[data_group.iloc[i]==j+1].sum()/(data_group.iloc[i]==j+1).sum()for j in range(group_number)]for i in range(len(rtn)-1)]
data_rtn_group=DF(data_rtn_group,index=rtn.index[1:],columns=list('12345'))
mimick=data_rtn_group['1']-data_rtn_group['5']
mimick.mean()/mimick.std()*np.sqrt(len(mimick))
mimick.to_pickle('F:/data/xccdata/essay/SMB_free_daily')


from datetime import datetime
t0=datetime.strptime('2005-02-01','%Y-%m-%d')
t1=datetime.strptime('2017-09-29','%Y-%m-%d')
book=pd.read_pickle('F:/data/xccdata/essay/book')/10000
book=book.resample('D').first().ffill()
book=pd.DataFrame(book,index=data_adj_close.index,columns=data_adj_close.columns)
book.to_pickle('F:/data/xccdata/essay/book_daily')
size=data['size_tot'].unstack()
BM=book/size
BM=BM.loc[t0:t1]
rtn=rtn.loc[t0:t1]
data_group=DF([pd.qcut(BM.iloc[i],q=percentile,labels=group_label) for i in range(len(BM))],index=BM.index,columns=BM.columns)
data_rtn_group=[[rtn.iloc[i+1].loc[data_group.iloc[i]==j+1].sum()/(data_group.iloc[i]==j+1).sum()for j in range(group_number)]for i in range(len(rtn)-1)]
data_rtn_group=DF(data_rtn_group,index=rtn.index[1:],columns=list('12345'))
mimick=data_rtn_group['5']-data_rtn_group['1']
mimick.mean()/mimick.std()*np.sqrt(len(mimick))
mimick.to_pickle('F:/data/xccdata/essay/HML_tot_daily')























# 计算 三因子 alpha

import pandas as pd
from pandas import DataFrame as DF
from pandas import Series as SS
import numpy as np
from numpy.linalg import pinv

from datetime import datetime
t0=datetime.strptime('2005-02-01','%Y-%m-%d')
t1=datetime.strptime('2017-09-29','%Y-%m-%d')

data_adj_close=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_close']].unstack()
rtn=(data_adj_close/data_adj_close.shift(1)-1)*100
SMB=pd.read_pickle('F:/data/xccdata/essay/SMB_tot_daily') #153
SMB.name='SMB'
HML=pd.read_pickle('F:/data/xccdata/essay/HML_tot_daily')  #153
HML.name='HML'
market=pd.read_pickle('F:/data/xccdata/essay/index_hs300_daily')  #154
market.name='market'


c=pd.Series(1,index=market.index)
X=pd.concat([c,market,SMB,HML],axis=1)
X=X.loc[t0:t1]
rtn=rtn.loc[t0:t1]
rtn.index=(rtn.index.year).astype(str)+'-'+(rtn.index.month).astype(str).str.zfill(2)
rtn.index.name='date'
rtn=rtn.stack()
#rtn['date']=(rtn.index.get_level_values(0).year).astype(str)+'-'+(rtn.index.get_level_values(0).month).astype(str).str.zfill(2)
#rtn=rtn.set_index(['date',rtn.index.get_level_values(1)])

def func(arr):
    z= arr.values!=np.nan
    if z.sum()<=10:
        return 0.0
    else:
        x=(X[arr.index.get_level_values(0)[0]]).values[z,:]
        return pinv(x.T@x)@x.T@(arr.values[z])[0]
alpha=rtn[:100].groupby(['date','stkcd']).apply(func)


rtn.index.get_level_values(0).year








