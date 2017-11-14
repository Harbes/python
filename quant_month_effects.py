import pandas as pd
import numpy as np
from pandas import DataFrame
from datetime import datetime
import matplotlib.pyplot as plt

# 整理数据
#data=pd.read_pickle('F:/data/xccdata/PV_datetime')[['adj_open','adj_close']]
data=pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')[['amount','opnprc','clsprc','adj_open','adj_close','size_tot','size_free']]
filter_=pd.read_pickle('/Users/harbes/data/xccdata/filter') # 4672606个有效数据点(原来有6140094个数据点)
opnprc=data['adj_open'].unstack()
clsprc=data['adj_close'].unstack()
rtn=((clsprc-opnprc)/opnprc)[filter_==1]
rtn.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除


after_fes_data=[
    (datetime(2005, 2, 16), '2005'),
    (datetime(2006, 2, 6), '2006'),
    (datetime(2007, 2, 26), '2007'),
    (datetime(2008, 2, 13), '2008'),
    (datetime(2009, 2, 2), '2009'),
    (datetime(2010, 2, 22), '2010'),
    (datetime(2011, 2, 9), '2011'),
    (datetime(2012, 1, 30), '2012'),
    (datetime(2013, 2, 18), '2013'),
    (datetime(2014, 2, 7), '2014'),
    (datetime(2015, 2, 25), '2015'),
    (datetime(2016, 2, 15), '2016'),
    (datetime(2017, 2, 3), '2017'),
]



# 调取二月数据
rtn_feb=pd.concat([rtn['{}-02'.format(i)] for i in range(2005,2018)])
plt.plot((rtn_feb.mean(axis=1)+1).values.cumprod())
#rtn_feb.groupby(rtn_feb.index.year).describe()
#rtn_feb.mean(axis=1).groupby(rtn_feb.index.year).describe()

# 春节后n个有效交易日
rtn_after_fes=pd.concat([rtn.iloc[rtn.index.get_loc(i[0])-3:rtn.index.get_loc(i[0])+8] for i in after_fes_data])
plt.plot((rtn_after_fes.mean(axis=1)+1).values.cumprod(axis=0))
plt.show()
# 春节前n个有效交易日
rtn_before_fes=pd.concat([rtn.iloc[rtn.index.get_loc(i[0])-18:rtn.index.get_loc(i[0])-3] for i in after_fes_data])
plt.plot((rtn_before_fes.mean(axis=1)+1).values.cumprod(axis=0))

tmp=rtn_after_fes.mean(axis=1);tmp.mean()/tmp.std()*np.sqrt(len(tmp))




# 按size进行分组(每天)
size=data['size_tot'].unstack()[filter_==1]#.groupby(key).mean()
size.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
rtn=rtn.loc[size.index]
num_by_size=10
label_size=[i+1 for i in range(num_by_size)] # 1表示小市值，2表示大市值
percentile=np.linspace(0,1,num_by_size+1)
mark_size=DataFrame([pd.qcut(size.iloc[i],q=percentile,labels=label_size) for i in range(len(size))],index=size.index,columns=size.columns)
rtn_by_size=DataFrame([[rtn.iloc[i].loc[mark_size.iloc[i-1]==k].mean() for k in label_size] for i in range(1,len(rtn))],index=rtn.index[1:],columns=label_size)

# 春节后
rtn_after_fes=pd.concat([rtn_by_size.iloc[rtn_by_size.index.get_loc(i[0])-5:rtn_by_size.index.get_loc(i[0])+10] for i in after_fes_data])
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
# 春节前
rtn_before_fes=pd.concat([rtn_by_size.iloc[rtn_by_size.index.get_loc(i[0])-10:rtn_by_size.index.get_loc(i[0])-5] for i in after_fes_data])
plt.plot((rtn_before_fes+1).values.cumprod(axis=0))




# 按illiq进行分组(每天)【结果是与illiq无关】
amount=data['amount'].unstack()[filter_==1]#.groupby(key).mean()
opnprc=data['adj_open'].unstack()[filter_==1]
clsprc=data['adj_close'].unstack()[filter_==1]

illiq=np.abs((clsprc-opnprc)/opnprc/amount)*1e5
illiq.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
illiq[illiq==0]=np.nan
rtn=rtn.loc[illiq.index]
num_by_=10
label_illiq=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
mark_illiq=DataFrame([pd.qcut(illiq.iloc[i],q=percentile,labels=label_illiq) for i in range(len(illiq))],index=illiq.index,columns=illiq.columns)
rtn_by_illiq=DataFrame([[rtn.iloc[i].loc[mark_illiq.iloc[i-1]==k].mean() for k in label_illiq] for i in range(1,len(rtn))],index=rtn.index[1:],columns=label_illiq)

# 春节后
rtn_after_fes=pd.concat([rtn_by_illiq.iloc[rtn_by_illiq.index.get_loc(i[0])-5:rtn_by_illiq.index.get_loc(i[0])+10] for i in after_fes_data])
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
# 春节前
rtn_before_fes=pd.concat([rtn_by_illiq.iloc[rtn_by_illiq.index.get_loc(i[0])-10:rtn_by_illiq.index.get_loc(i[0])-5] for i in after_fes_data])
plt.plot((rtn_before_fes+1).values.cumprod(axis=0))




# 按过去收益(momentum or contrarian)进行分组(每天)【结果与illiq无关】
rtn[rtn==0]=np.nan
num_by_=10
label_past=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
mark_past=DataFrame([pd.qcut(rtn.iloc[i],q=percentile,labels=label_past) for i in range(len(rtn))],index=rtn.index,columns=rtn.columns)
rtn_by_past=DataFrame([[rtn.iloc[i].loc[mark_past.iloc[i-1]==k].mean() for k in label_past] for i in range(1,len(rtn))],index=rtn.index[1:],columns=label_past)

# 春节后
rtn_after_fes=pd.concat([rtn_by_past.iloc[rtn_by_past.index.get_loc(i[0])-5:rtn_by_past.index.get_loc(i[0])+10] for i in after_fes_data])
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
# 春节前
rtn_before_fes=pd.concat([rtn_by_past.iloc[rtn_by_past.index.get_loc(i[0])-8:rtn_by_past.index.get_loc(i[0])-5] for i in after_fes_data])
plt.plot((rtn_before_fes+1).values.cumprod(axis=0))




# 按size分组(前n个交易日数据) ; y
num_by_=10
label_=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
rtn_after_fes=DataFrame(0,index=np.array(after_fes_data)[:,0],columns=label_)
indi=data['size_tot'].unstack()[filter_==1]
indi.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
opnprc=opnprc.loc[indi.index]
clsprc=clsprc.loc[indi.index]
mark_=DataFrame([pd.qcut(indi.iloc[np.maximum(i-30,0):i].mean(),q=percentile,labels=label_) for i in range(1,len(indi))],index=rtn.index[1:],columns=rtn.columns)
for k,i in enumerate(after_fes_data):
    rtn_tmp=(clsprc.iloc[clsprc.index.get_loc(i[0])+10]-opnprc.iloc[opnprc.index.get_loc(i[0])-5])/opnprc.iloc[opnprc.index.get_loc(i[0])-5]
    for j in range(num_by_):
        rtn_after_fes.iloc[k,j]=rtn_tmp[mark_.iloc[mark_.index.get_loc(i[0])]==j+1].mean()
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
tmp=rtn_after_fes[10]-rtn_after_fes[1];tmp.mean()/tmp.std()*np.sqrt(len(tmp))



# 按illiq分组(前n个交易日数据) ; y
num_by_=10
label_=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
rtn_after_fes=DataFrame(0,index=np.array(after_fes_data)[:,0],columns=label_)
amount=data['amount'].unstack()[filter_==1]#.groupby(key).mean()
opnprc=data['adj_open'].unstack()[filter_==1]
clsprc=data['adj_close'].unstack()[filter_==1]
indi=np.abs((clsprc-opnprc)/opnprc/amount)*1e5
indi.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
opnprc=opnprc.loc[indi.index]
clsprc=clsprc.loc[indi.index]
mark_=DataFrame([pd.qcut(indi.iloc[np.maximum(i-30,0):i].mean(),q=percentile,labels=label_) for i in range(1,len(indi))],index=rtn.index[1:],columns=rtn.columns)
for k,i in enumerate(after_fes_data):
    rtn_tmp=(clsprc.iloc[clsprc.index.get_loc(i[0])+10]-opnprc.iloc[opnprc.index.get_loc(i[0])-5])/opnprc.iloc[opnprc.index.get_loc(i[0])-5]
    for j in range(num_by_):
        rtn_after_fes.iloc[k,j]=rtn_tmp[mark_.iloc[mark_.index.get_loc(i[0])]==j+1].mean()
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
tmp=rtn_after_fes[10]-rtn_after_fes[1];tmp.mean()/tmp.std()*np.sqrt(len(tmp))





# 按past分组(前n个交易日数据) ; y  ; 资本流入 or 流出 ？
num_by_=10
label_=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
rtn_after_fes=DataFrame(0,index=np.array(after_fes_data)[:,0],columns=label_)
opnprc=data['adj_open'].unstack()[filter_==1]
clsprc=data['adj_close'].unstack()[filter_==1]
indi=(clsprc-opnprc)/opnprc
indi.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
opnprc=opnprc.loc[indi.index]
clsprc=clsprc.loc[indi.index]
mark_=DataFrame([pd.qcut(indi.iloc[np.maximum(i-30,0):i].mean(),q=percentile,labels=label_) for i in range(1,len(indi))],index=rtn.index[1:],columns=rtn.columns)
for k,i in enumerate(after_fes_data):
    rtn_tmp=(clsprc.iloc[clsprc.index.get_loc(i[0])+10]-opnprc.iloc[opnprc.index.get_loc(i[0])-5])/opnprc.iloc[opnprc.index.get_loc(i[0])-5]
    for j in range(num_by_):
        rtn_after_fes.iloc[k,j]=rtn_tmp[mark_.iloc[mark_.index.get_loc(i[0])]==j+1].mean()
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
tmp=rtn_after_fes[10]-rtn_after_fes[1];tmp.mean()/tmp.std()*np.sqrt(len(tmp))





# 按price分组(前n个交易日数据) ; n
num_by_=10
label_=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
rtn_after_fes=DataFrame(0,index=np.array(after_fes_data)[:,0],columns=label_)
indi=data['opnprc'].unstack()[filter_==1]
indi.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
opnprc=opnprc.loc[indi.index]
clsprc=clsprc.loc[indi.index]
mark_=DataFrame([pd.qcut(indi.iloc[np.maximum(i-30,0):i].mean(),q=percentile,labels=label_) for i in range(1,len(indi))],index=rtn.index[1:],columns=rtn.columns)
for k,i in enumerate(after_fes_data):
    rtn_tmp=(clsprc.iloc[clsprc.index.get_loc(i[0])+10]-opnprc.iloc[opnprc.index.get_loc(i[0])-5])/opnprc.iloc[opnprc.index.get_loc(i[0])-5]
    for j in range(num_by_):
        rtn_after_fes.iloc[k,j]=rtn_tmp[mark_.iloc[mark_.index.get_loc(i[0])]==j+1].mean()
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
tmp=rtn_after_fes[10]-rtn_after_fes[1];tmp.mean()/tmp.std()*np.sqrt(len(tmp))




# 按delta_price * size 分组(前n个交易日数据) ; 非线性关系：组合收益先增后减
num_by_=10
label_=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
rtn_after_fes=DataFrame(0,index=np.array(after_fes_data)[:,0],columns=label_)
opnprc=data['adj_open'].unstack()[filter_==1]
clsprc=data['adj_close'].unstack()[filter_==1]
size=data['size_free'].unstack()[filter_==1]
indi=(clsprc-opnprc)*size # opnprc-clsprc#
indi.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
opnprc=opnprc.loc[indi.index]
clsprc=clsprc.loc[indi.index]
mark_=DataFrame([pd.qcut(indi.iloc[np.maximum(i-30,0):i].mean(),q=percentile,labels=label_) for i in range(1,len(indi))],index=indi.index[1:],columns=indi.columns)
for k,i in enumerate(after_fes_data):
    rtn_tmp=(clsprc.iloc[clsprc.index.get_loc(i[0])+10]-opnprc.iloc[opnprc.index.get_loc(i[0])-5])/opnprc.iloc[opnprc.index.get_loc(i[0])-5]
    for j in range(num_by_):
        rtn_after_fes.iloc[k,j]=rtn_tmp[mark_.iloc[mark_.index.get_loc(i[0])]==j+1].mean()
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
tmp=rtn_after_fes[5]-rtn_after_fes[1];tmp.mean()/tmp.std()*np.sqrt(len(tmp))

rtn_after_fes.mean(axis=0)




# 按delta_price * 资金净流入 分组(前n个交易日数据) ;
num_by_=10
label_=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
rtn_after_fes=DataFrame(0,index=np.array(after_fes_data)[:,0],columns=label_)
opnprc=data['adj_open'].unstack()[filter_==1]
clsprc=data['adj_close'].unstack()[filter_==1]
amount=data['amount'].unstack()[filter_==1]
indi=(clsprc-opnprc)*amount/(clsprc+opnprc) # opnprc-clsprc#
indi.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
opnprc=opnprc.loc[indi.index]
clsprc=clsprc.loc[indi.index]
mark_=DataFrame([pd.qcut(indi.iloc[np.maximum(i-30,0):i].mean(),q=percentile,labels=label_) for i in range(1,len(indi))],index=indi.index[1:],columns=indi.columns)
for k,i in enumerate(after_fes_data):
    rtn_tmp=(clsprc.iloc[clsprc.index.get_loc(i[0])+10]-opnprc.iloc[opnprc.index.get_loc(i[0])-5])/opnprc.iloc[opnprc.index.get_loc(i[0])-5]
    for j in range(num_by_):
        rtn_after_fes.iloc[k,j]=rtn_tmp[mark_.iloc[mark_.index.get_loc(i[0])]==j+1].mean()
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
tmp=rtn_after_fes[10]-rtn_after_fes[1];tmp.mean()/tmp.std()*np.sqrt(len(tmp))

rtn_after_fes.mean(axis=0)





# 按 signed-volume 分组(前n个交易日数据) ; y (前10个交易日资金净流入)
num_by_=10
label_=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
rtn_after_fes=DataFrame(0,index=np.array(after_fes_data)[:,0],columns=label_)
opnprc=data['adj_open'].unstack()[filter_==1]
clsprc=data['adj_close'].unstack()[filter_==1]
amount=data['amount'].unstack()[filter_==1]
indi=np.sign(clsprc-opnprc)*amount/(clsprc+opnprc) # opnprc-clsprc#
indi.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
opnprc=opnprc.loc[indi.index]
clsprc=clsprc.loc[indi.index]
mark_=DataFrame([pd.qcut(indi.iloc[np.maximum(i-15,0):i].mean(),q=percentile,labels=label_) for i in range(1,len(indi))],index=indi.index[1:],columns=indi.columns)
for k,i in enumerate(after_fes_data):
    rtn_tmp=(clsprc.iloc[clsprc.index.get_loc(i[0])+10]-opnprc.iloc[opnprc.index.get_loc(i[0])-5])/opnprc.iloc[opnprc.index.get_loc(i[0])-5]
    for j in range(num_by_):
        rtn_after_fes.iloc[k,j]=rtn_tmp[mark_.iloc[mark_.index.get_loc(i[0])]==j+1].mean()
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
tmp=rtn_after_fes[10]-rtn_after_fes[1];tmp.mean()/tmp.std()*np.sqrt(len(tmp))
rtn_after_fes.mean(axis=0)





# 按 size_free 分组(前n个交易日数据) ; y (前10个交易日资金净流入)
num_by_=10
label_=[i+1 for i in range(num_by_)] #
percentile=np.linspace(0,1,num_by_+1)
rtn_after_fes=DataFrame(0,index=np.array(after_fes_data)[:,0],columns=label_)
opnprc=data['adj_open'].unstack()[filter_==1]
clsprc=data['adj_close'].unstack()[filter_==1]
size=data['size_free'].unstack()[filter_==1]
indi=size-size.shift(1)
indi.drop_duplicates(keep=False,inplace=True) # 有一些日期全是nan，需要剔除
opnprc=opnprc.loc[indi.index]
clsprc=clsprc.loc[indi.index]
mark_=DataFrame([pd.qcut(indi.iloc[np.maximum(i-15,0):i].mean(),q=percentile,labels=label_) for i in range(1,len(indi))],index=indi.index[1:],columns=indi.columns)
for k,i in enumerate(after_fes_data):
    rtn_tmp=(clsprc.iloc[clsprc.index.get_loc(i[0])+10]-opnprc.iloc[opnprc.index.get_loc(i[0])-5])/opnprc.iloc[opnprc.index.get_loc(i[0])-5]
    for j in range(num_by_):
        rtn_after_fes.iloc[k,j]=rtn_tmp[mark_.iloc[mark_.index.get_loc(i[0])]==j+1].mean()
plt.plot((rtn_after_fes+1).values.cumprod(axis=0))
tmp=rtn_after_fes[10]-rtn_after_fes[1];tmp.mean()/tmp.std()*np.sqrt(len(tmp))

rtn_after_fes.mean(axis=0)






















########### appendix ###################
# 二月份且春节后的第一个交易日
festival_data=[
    (datetime(2005, 2, 16), '2005'),
    (datetime(2006, 2, 6), '2006'),
    (datetime(2007, 2, 26), '2007'),
    (datetime(2008, 2, 13), '2008'),
    (datetime(2009, 2, 2), '2009'),
    (datetime(2010, 2, 22), '2010'),
    (datetime(2011, 2, 9), '2011'),
    (datetime(2012, 2, 2), '2012'),
    (datetime(2013, 2, 18), '2013'),
    (datetime(2014, 2, 7), '2014'),
    (datetime(2015, 2, 25), '2015'),
    (datetime(2016, 2, 15), '2016'),
    (datetime(2017, 2, 3), '2017'),
]
rtn_feb_size['range']=range(len(rtn_feb_size))
rtn_feb_size.index=range(len(rtn_feb_size))

# 春节后首个交易日以及春节前最后一个交易日
after_fes_data=[
    (datetime(2005, 2, 16), '2005'),
    (datetime(2006, 2, 6), '2006'),
    (datetime(2007, 2, 26), '2007'),
    (datetime(2008, 2, 13), '2008'),
    (datetime(2009, 2, 2), '2009'),
    (datetime(2010, 2, 22), '2010'),
    (datetime(2011, 2, 9), '2011'),
    (datetime(2012, 1, 30), '2012'),
    (datetime(2013, 2, 18), '2013'),
    (datetime(2014, 2, 7), '2014'),
    (datetime(2015, 2, 25), '2015'),
    (datetime(2016, 2, 15), '2016'),
    (datetime(2017, 2, 3), '2017'),
]
before_data=[
    (datetime(2005, 2, 3), '2005'),
    (datetime(2006, 1, 26), '2006'),
    (datetime(2007, 2, 16), '2007'),
    (datetime(2008, 2, 4), '2008'),
    (datetime(2009, 1, 23), '2009'),
    (datetime(2010, 2, 12), '2010'),
    (datetime(2011, 1, 31), '2011'),
    (datetime(2012, 1, 20), '2012'),
    (datetime(2013, 2, 8), '2013'),
    (datetime(2014, 1, 30), '2014'),
    (datetime(2015, 2, 17), '2015'),
    (datetime(2016, 2, 3), '2016'),
    (datetime(2017, 1, 26), '2017'),
]

rtn_feb_size.index.get_loc((datetime(2017, 2, 3)))

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
