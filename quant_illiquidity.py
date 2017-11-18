import pandas as pd
import numpy as np
from pandas import DataFrame

# 整理数据
#data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_open','adj_close']]
data=pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')[['amount','opnprc','clsprc','adj_open','adj_close','size_tot','size_free']]

amount=data['amount'].unstack()
#amount[amount==0]=np.nan # 将交易额为0的数据替换成NaN，并检验 np.invert(amount==0).all().all()
amount_filter=pd.read_pickle('/Users/harbes/data/xccdata/amount_filter') # 1 表示成交量满足要求的股票数据
amount=amount[amount_filter==1] # 剔除成交量较小的股票数据；这部分的剔除大约可以将最优组合的累计收益下降60%
amount /=1e5

opnprc=data['adj_open'].unstack()
clsprc=data['adj_close'].unstack()
limit_move=pd.read_pickle('/Users/harbes/data/xccdata/limit_move') # True表示涨跌停；剔除涨跌停的影响程度与剔除小额交易相近
ST=pd.read_pickle('/Users/harbes/data/xccdata/ST_datetime') # N 表示非ST，非新股 ；剔除ST以及新股后，基本与剔除所有因素一致
NT=pd.read_pickle('/Users/harbes/data/xccdata/NT_datetime') # 1 表示停牌；剔除NT的效果并不很显著
NT_filter=pd.read_pickle('/Users/harbes/data/xccdata/NT_filter_120') # 1 表示停牌超过一定时间的数据点;
clsprc=clsprc[NT==0][NT_filter==0][~limit_move][ST=='N'] # NT和ST大概能剔除最优组合的2-3倍涨幅
key = lambda x: x.year * 100 + x.month
price0 = opnprc.groupby(key).first() # 买入价
price1 = clsprc.groupby(key).last() # 卖出价
price0.index = pd.to_datetime(price0.index.values.astype(str),format=('%Y%m'))
price1.index = pd.to_datetime(price1.index.values.astype(str),format=('%Y%m'))
rtn = (price1-price0)/price0 # 月度收益数据

#opnprc=data['opnprc'].unstack();clsprc=data['clsprc'].unstack() # 是否复权，对分组结果影响不大
illiq=np.abs((clsprc-opnprc)/opnprc/amount) # 事实上也可以使用turnover(其倒数)，如下：
#illiq=data['size_tot'].unstack()/amount
illiq_m=illiq.groupby(key).mean()
illiq_m.index = pd.to_datetime(illiq_m.index.values.astype(str),format=('%Y%m'))




group_num=10
percentile=np.linspace(0,1,group_num+1)
label_illiq=[i+1 for i in range(group_num)] # 1表示流动性较好的组合，另一个极端组合是流动性较差的组合

# 分组，并计算不同组合收益
mark_illiq=DataFrame([pd.qcut(illiq_m.iloc[i],q=percentile,labels=label_illiq) for i in range(len(illiq_m))],index=illiq_m.index)
illiquidity=DataFrame([[rtn.iloc[i][mark_illiq.iloc[i-1]==k].mean() for k in label_illiq] for i in range(1,len(rtn))],index=rtn.index[1:],columns=label_illiq)


(illiquidity+1)[np.arange(1,11)].cumprod().plot()
zero_cost=illiquidity[10]-illiquidity[1]
zero_cost.mean()/zero_cost.std()*np.sqrt(len(zero_cost))

# 检查size，结果发现size确实对结果有影响，所以需要double-sort
size=data['size_tot'].unstack().groupby(key).mean()
size.index = pd.to_datetime(size.index.values.astype(str),format=('%Y%m'))
size=DataFrame([[size.iloc[i][mark_illiq.iloc[i]==k].mean() for k in label_illiq] for i in range(len(size))],index=size.index,columns=label_illiq)
size.mean()/1e5




# double-sort,first by size,then by illiq
group_num=5
percentile=np.linspace(0,1,group_num+1)
label_size=[i+1 for i in range(group_num)] # 1表示流动性较好的组合，另一个极端组合是流动性较差的组合
label_illiq=[i+1 for i in range(group_num)] # 1表示流动性较好的组合，另一个极端组合是流动性较差的组合

filter_=pd.read_pickle('/Users/harbes/data/xccdata/filter') # 4672606个有效数据点(原来有6140094个数据点)
size=data['size_tot'].unstack()[filter==1].groupby(key).mean()
size.index = pd.to_datetime(size.index.values.astype(str),format=('%Y%m'))

mark_size=DataFrame([pd.qcut(size.iloc[i],q=percentile,labels=label_size) for i in range(len(size))],index=size.index)
mark_illiq=DataFrame(np.nan,index=illiq_m.index,columns=illiq_m.columns)
for l_i in label_size:
    tmp=DataFrame([pd.qcut(illiq_m.iloc[i][mark_size.iloc[i]==l_i],q=percentile,labels=label_illiq) for i in range(len(illiq_m))],index=illiq_m.index)
    mark_illiq=mark_illiq.combine_first(tmp)

size_illiquidity=DataFrame(np.zeros((25,len(size)-1)),index=pd.MultiIndex.from_product([label_size,label_illiq]))
s_i_count=DataFrame(np.zeros((25,len(size)-1)),index=pd.MultiIndex.from_product([label_size,label_illiq])) # 用于检查分组是否均匀
for s in label_size:
    for i in label_illiq:
        #size_illiquidity.loc[(s,i),:]= pd.Series([rtn.iloc[j][mark_size.iloc[j-1]==s][mark_illiq.iloc[j-1]==i].mean() for j in range(1,len(rtn))])
        size_illiquidity.loc[(s,i),:]= pd.Series([rtn.iloc[j][np.logical_and(mark_size.iloc[j-1]==s,mark_illiq.iloc[j-1]==i)].mean() for j in range(1,len(rtn))])
        s_i_count.loc[(s, i), :]=pd.Series([rtn.iloc[j][mark_size.iloc[j-1]==s][mark_illiq.iloc[j-1]==i].count() for j in range(1,len(rtn))])
size_illiquidity.columns=rtn.index[1:]
# 结果显示
(size_illiquidity+1).loc[(1,slice(None)),slice(None)].T.cumprod().plot() # 为什么使用axis=1不能得到想要的结果？？？
(size_illiquidity+1).loc[(1,[1,2,3,4]),slice(None)].T.cumprod().plot() # 为什么使用axis=1不能得到想要的结果？？？
s_i_count.max()-s_i_count.min() # 检查分组是否均匀





# double-sort,size & illiq sorted independently
group_num=5
percentile=np.linspace(0,1,group_num+1)
label_size=[i+1 for i in range(group_num)] # 1表示流动性较好的组合，另一个极端组合是流动性较差的组合
label_illiq=[i+1 for i in range(group_num)] # 1表示流动性较好的组合，另一个极端组合是流动性较差的组合

size=data['size_tot'].unstack().groupby(key).mean()
size.index = pd.to_datetime(size.index.values.astype(str),format=('%Y%m'))

mark_size=DataFrame([pd.qcut(size.iloc[i],q=percentile,labels=label_size) for i in range(len(size))],index=size.index)
mark_illiq=DataFrame([pd.qcut(illiq_m.iloc[i],q=percentile,labels=label_illiq) for i in range(len(illiq_m))],index=illiq_m.index)

size_illiquidity=DataFrame(np.zeros((25,len(size)-1)),index=pd.MultiIndex.from_product([label_size,label_illiq]))
s_i_count=DataFrame(np.zeros((25,len(size)-1)),index=pd.MultiIndex.from_product([label_size,label_illiq])) # 用于检查分组是否均匀
for s in label_size:
    for i in label_illiq:
        mark_illiq_tmp=DataFrame
        size_illiquidity.loc[(s,i),:]= pd.Series([rtn.iloc[j][mark_size.iloc[j-1]==s][mark_illiq.iloc[j-1]==i].mean() for j in range(1,len(rtn))])
        s_i_count.loc[(s, i), :]=pd.Series([rtn.iloc[j][mark_size.iloc[j-1]==s][mark_illiq.iloc[j-1]==i].count() for j in range(1,len(rtn))])

# 结果的图形显示
(size_illiquidity+1).loc[(1,[2,3,4,5]),slice(None)].T.cumprod().plot() # 为什么使用axis=1不能得到想要的结果？？？

s_i_count.max()-s_i_count.min()
s_i_count.min()
s_i_count.loc[(1,5)]
s_i_count.loc[(1,1)] # 一个非常奇怪的、不符合预期的组合

























