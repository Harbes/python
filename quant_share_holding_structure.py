import pandas as pd
from pandas import DataFrame
import numpy as np
from pandas.tseries.offsets import YearEnd
from datetime import datetime

#data_path='/Users/harbes/data/xccdata'
data_path='F:/data/xccdata'
after_fes_data = [
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

institution = pd.read_csv(data_path+'/insti_holding.csv').iloc[1:]
institution['Date'] = pd.to_datetime(institution['Date'])
institution = institution.set_index(['Date', 'SCode'])
institution.index.names = ['trddt', 'stkcd']
institution = institution.astype(np.float64)
institution = institution.sum(axis=1)
institution = institution.unstack()
institution[institution > 100] = 100
# institution.head()
# (institution==0).sum(axis=1) # 机构不持股
# (institution==100).sum(axis=1) # 完全由机构持股

government = pd.read_csv(data_path+'/govern_holding.csv').iloc[1:]
government['Year'] = pd.to_datetime(government['Year']) + YearEnd()
government = government.set_index(['Year', 'SCode'])
government = government.unstack()
government.columns = government.columns.get_level_values(1)
government.index.names = ['trddt']
government.columns.names = ['stkcd']
government[government.isnull()] = np.nan
government = government.astype(np.float64)
#government.head()
# (government==0).sum(axis=1) # 政府不持股

# non_individual = government + institution
individual = 100 - institution
# (individual==100).sum(axis=1) # 完全由个人持股


PV_datetime = pd.read_pickle(data_path+'/PV_datetime')[['size_tot', 'adj_close', 'adj_open', 'clsprc']]
filter_ = pd.read_pickle(data_path+'/filter')
size = PV_datetime['size_tot'].unstack()
size.columns = size.columns.str.slice(0, 6)
adj_open = PV_datetime['adj_open'].unstack()
adj_open.columns = adj_open.columns.str.slice(0, 6)
adj_close = PV_datetime['adj_close'].unstack()
adj_close.columns = adj_close.columns.str.slice(0, 6)
clsprc = PV_datetime['clsprc'].unstack()
clsprc.columns = clsprc.columns.str.slice(0, 6)

book = pd.read_pickle(data_path+'/BS')[['ann_dt', 'stkcd', 'tot_assets', 'tot_liab']]  #
book['stkcd'] = book['stkcd'].str.slice(0, 6)
book = book[~np.isnan(book['ann_dt'])]
book.drop_duplicates(book[['ann_dt', 'stkcd']], inplace=True)
book = book.set_index(['ann_dt', 'stkcd']).sort_index()
book['book'] = book['tot_assets'] - book['tot_liab']
book = book['book'].unstack()
book.index = book.index.astype(int).astype(str).to_datetime('%Y%m%d')
book.index.name = 'trddt'
book = book.resample('D').first().ffill()
# book.head()
BM = (book / size / 10000.0).loc[size.index]

institution = institution.resample('D').first().ffill()
institution = institution.loc[size.index & institution.index]
# institution.tail()

individual = individual.resample('D').first().ffill()
individual = individual.loc[size.index & individual.index]
#individual.tail()

# size-insti sort (evenly sort)
num_by_ = 5;
label_1 = [i + 1 for i in range(num_by_)];
percentile_1 = np.linspace(0, 1, num_by_ + 1)
label_2 = [1, 2, 3];
percentile_2 = (0.0, 0.3, 0.7, 1.0)
#label_2 = [i + 1 for i in range(num_by_)];percentile_2 = np.linspace(0, 1, num_by_ + 1)
# rtn=((adj_close-adj_open)/adj_open)[filter_[adj_close.columns]==1]
indicator1 = size.shift(1).loc['2006':'2016'][size.columns & institution.columns];  # 可以改为其他指标，例如BM
indicator2 = institution.loc['2006':'2016'][size.columns & institution.columns]
mark_1 = DataFrame([pd.qcut(indicator1.loc[i], q=percentile_1, labels=label_1) for i in
                    indicator1.index])  # ,index=indicator1.index,columns=indicator1.columns)
mark_2 = DataFrame(np.nan, index=mark_1.index, columns=mark_1.columns)
# mark_2[institution==0]=1
for l_ in label_1:
    tmp = DataFrame(
        [pd.qcut(indicator2.loc[i][mark_1.loc[i] == l_], q=percentile_2, labels=label_2) for i in indicator2.index],
        index=indicator2.index)
    mark_2 = mark_2.combine_first(tmp)


# 计算不同组合institution holding、size的平均值
holding = institution[size.columns & institution.columns]
size = size[size.columns & institution.columns]
average_holding = DataFrame(np.zeros((len(label_1) * len(label_2), 11)),
                            index=pd.MultiIndex.from_product([label_1, label_2]),
                            columns=np.array(after_fes_data)[1:-1, 0])
average_size = DataFrame(np.zeros((len(label_1) * len(label_2), 11)),
                         index=pd.MultiIndex.from_product([label_1, label_2]),
                         columns=np.array(after_fes_data)[1:-1, 0])
# select_stock = {}
for i in label_1:
    for j in label_2:
        for d in np.array(after_fes_data)[1:-1, 0]:
            average_holding.loc[(i, j), d] = holding.loc[d][
                (mark_1.loc[d] == i) & (mark_2.loc[d] == j)].mean()
            average_size.loc[(i, j), d] = size.loc[d][(mark_1.loc[d] == i) & (mark_2.loc[d] == j)].mean()
            #select_stock[((i, j), d)] = np.sort(np.random.choice(size.loc[d][(mark_1.loc[d] == i) & (mark_2.loc[d])].index, 25))

average_holding.mean(axis=1)
average_size.mean(axis=1)

# size-government sort (evenly sort)
government = government.resample('D').first().ffill()
government = government.loc[size.index & individual.index]
num_by_ = 5;
label_1 = [i + 1 for i in range(num_by_)];
percentile_1 = np.linspace(0, 1, num_by_ + 1)
label_2 = (1, 2, 3);
percentile_2 = (0.0, 0.1, 1.0)
indicator1 = size.shift(1).loc['2006':'2016'][size.columns & institution.columns];  # 可以改为其他指标，例如BM
indicator2 = government.loc['2006':'2016'][size.columns & institution.columns]
mark_1 = DataFrame([pd.qcut(indicator1.loc[i], q=percentile_1, labels=label_1) for i in
                    indicator1.index])  # ,index=indicator1.index,columns=indicator1.columns)
mark_2 = DataFrame(np.nan, index=mark_1.index, columns=mark_1.columns)
mark_2[government == 0] = 1
for l_ in label_1:
    tmp = DataFrame(
        [pd.qcut(indicator2.loc[i][mark_2.loc[i] != 1], q=percentile_2, labels=label_2[1:]) for i in indicator2.index],
        index=indicator2.index)
    mark_2 = mark_2.combine_first(tmp)

# 计算不同组合在春节前后的收益
adj_open = adj_open[size.columns & institution.columns]
adj_close = adj_close[size.columns & institution.columns]
average_return_before = DataFrame(np.zeros((len(label_1) * len(label_2), 11)),
                                  index=pd.MultiIndex.from_product([label_1, label_2]),
                                  columns=np.array(after_fes_data)[1:-1, 0])
average_return_after = DataFrame(np.zeros((len(label_1) * len(label_2), 11)),
                                 index=pd.MultiIndex.from_product([label_1, label_2]),
                                 columns=np.array(after_fes_data)[1:-1, 0])
T0 = -15  # 春节前6天股票也在涨，而且，似乎市值较大的股票涨的更多;大约在春节前15-前7天(前21-7天虽然数字更大，但是不够显著)，股票市场在跌，但是相对来说小市值股票表现更好一些
T1 = -7
t1 = 9  # past return
t0 = -6
for i in label_1:
    for j in label_2:
        for d in np.array(after_fes_data)[1:-1, 0]:
            average_return_after.loc[(i, j), d] = \
                ((adj_close.iloc[adj_close.index.get_loc(d) + t1] - adj_open.iloc[adj_open.index.get_loc(d) + t0])
                 / adj_open.iloc[adj_open.index.get_loc(d) + t0])[(mark_1.iloc[mark_1.index.get_loc(d) + t0] == i) & (
                mark_2.iloc[mark_2.index.get_loc(d) + t0] == j)].mean()
            average_return_before.loc[(i, j), d] = \
                ((adj_close.iloc[adj_close.index.get_loc(d) + T1] - adj_open.iloc[adj_open.index.get_loc(d) + T0])
                 / adj_open.iloc[adj_open.index.get_loc(d) + T0])[(mark_1.iloc[mark_1.index.get_loc(d) + T0] == i) & (
                    mark_2.iloc[mark_2.index.get_loc(d) + T0] == j)].mean()
average_return_after.mean(axis=1)
average_return_before.mean(axis=1)
n = 4  # 只有在大市值的股票中，不同institution_holding的股票才表现不同，似乎说明实际上机构和散户都在买，
tmp = average_return_after.loc[n, 1] - average_return_after.loc[n, 3];
tmp.mean() / tmp.std() * np.sqrt(len(tmp))
ttmp = average_return_before.loc[n, 1] - average_return_before.loc[n, 5];
ttmp.mean() / ttmp.std() * np.sqrt(len(ttmp))

n = 1
tmp = average_return_after.loc[1, n] - average_return_after.loc[5, n];
tmp.mean() / tmp.std() * np.sqrt(len(tmp))
tmp = average_return_before.loc[1, n] - average_return_before.loc[5, n];
tmp.mean() / tmp.std() * np.sqrt(len(tmp))

# past-return(-15~-7)与size(或者其他指标的double-sort)
adj_open = adj_open[size.columns & institution.columns]
adj_close = adj_close[size.columns & institution.columns]
num_by_ = 5;
label_1 = [i + 1 for i in range(num_by_)];
percentile_1 = np.linspace(0, 1, num_by_ + 1)
label_2 = [i + 1 for i in range(num_by_)];
percentile_2 = np.linspace(0, 1, num_by_ + 1)
t1 = 9  # past return
t0 = 0
indicator1 = (((adj_close - adj_open.shift(8)) / adj_open.shift(8)).shift(7)).loc['2006':'2016']
indicator2 = size.shift(1).loc['2006':'2016'][size.columns & institution.columns];
mark_1 = DataFrame([pd.qcut(indicator1.loc[i], q=percentile_1, labels=label_1) for i in
                    indicator1.index])  # ,index=indicator1.index,columns=indicator1.columns)
mark_2 = DataFrame(np.nan, index=mark_1.index, columns=mark_1.columns)
for l_ in label_1:
    tmp = DataFrame(
        [pd.qcut(indicator2.loc[i][mark_1.loc[i] == l_], q=percentile_2, labels=label_2) for i in indicator2.index],
        index=indicator2.index)
    mark_2 = mark_2.combine_first(tmp)

average_return = DataFrame(np.zeros((len(label_1) * len(label_2), 11)),
                           index=pd.MultiIndex.from_product([label_1, label_2]),
                           columns=np.array(after_fes_data)[1:-1, 0])
for i in label_1:
    for j in label_2:
        for d in np.array(after_fes_data)[1:-1, 0]:
            average_return.loc[(i, j), d] = \
                ((adj_close.iloc[adj_close.index.get_loc(d) + t1] - adj_open.iloc[adj_open.index.get_loc(d) + t0])
                 / adj_open.iloc[adj_open.index.get_loc(d) + t0])[(mark_1.iloc[mark_1.index.get_loc(d) + t0] == i) & (
                mark_2.iloc[mark_2.index.get_loc(d) + t0] == j)].mean()
average_return.mean(axis=1)
# (average_return.mean(axis=1)-average_return['2015-02-25'])/average_return['2015-02-25']

n = 1;
tmp = average_return.loc[n, 1] - average_return.loc[n, 5];
tmp.mean() / tmp.std() * np.sqrt(len(tmp))




# 观察不同 size-holding 股票的四种类型占比
## 选择股票
num_by_ = 5;
num_selected_per = 30
label_1 = [i + 1 for i in range(num_by_)];
percentile_1 = np.linspace(0, 1, num_by_ + 1)
label_2 = [1, 2, 3];
percentile_2 = (0.0, 0.3, 0.7, 1.0)
holding = institution[size.columns & institution.columns]
size = size[size.columns & institution.columns]
indicator1 = size.shift(1).loc['2006':'2016'][size.columns & institution.columns];  # 可以改为其他指标，例如BM
indicator2 = holding.loc['2006':'2016']
mark_1 = DataFrame([pd.qcut(indicator1.loc[i], q=percentile_1, labels=label_1) for i in
                    indicator1.index])  # ,index=indicator1.index,columns=indicator1.columns)
mark_2 = DataFrame(np.nan, index=mark_1.index, columns=mark_1.columns)
# mark_2[institution==0]=1
for l_ in label_1:
    tmp = DataFrame(
        [pd.qcut(indicator2.loc[i][mark_1.loc[i] == l_], q=percentile_2, labels=label_2) for i in indicator2.index],
        index=indicator2.index)
    mark_2 = mark_2.combine_first(tmp)
select_stock = {}
for i in label_1:
    for j in label_2:
        for d in np.array(after_fes_data)[-3:-2, 0]:
            select_stock[((i, j), d)] = tuple(
                np.sort(np.random.choice(size.loc[d][(mark_1.loc[d] == i) & (mark_2.loc[d])].index, num_selected_per)))
stocks = np.array([select_stock[i] for i in select_stock.keys()])
stocks = stocks.reshape((stocks.size,))
#DataFrame(select_stock).to_pickle(data_path+'/selected_stocks_by_size_institution')

# 从日度数据中提取交易量数据
import h5py
import os
import time

#stocks=pd.read_pickle(data_path+'/selected_stocks_by_size_institution')
t0 = time.time()
rootdir = data_path+'/bid_ask'
li_ = [i for i in os.listdir(rootdir) if not i.endswith('_') and not i.endswith('.h5')][5:]  # 列出文件夹下所有的目录与文件
trade_type = DataFrame(np.nan, index=pd.MultiIndex.from_product([li_, ['indi', 'M', 'insti'], ['buy', 'sell']]),
                       columns=stocks)
institution_standard = 1e4
individual_standard = 1.6e3
for d in li_:
    filename = rootdir + '/' + d
    f = h5py.File(filename, 'r')
    for stk in stocks:
        try:
            data = DataFrame(
                [list(f['stk'][stk]['volume']), list(f['stk'][stk]['trend']), list(f['stk'][stk]['lastPrc'])],
                index=['v', 't', 'p']).T  #
            data['v'] = data['v'] - data['v'].shift(1)
            trade_type.loc[(d, 'indi', 'buy'), stk] = data['v'][(data['v'] < individual_standard / data['p']) & (
            data['t'] > 0)].iloc[3:-3].sum()
            trade_type.loc[(d, 'indi', 'sell'), stk] = data['v'][(data['v'] < individual_standard / data['p']) & (
            data['t'] < 0)].iloc[3:-3].sum()
            trade_type.loc[(d, 'M', 'buy'), stk] = data['v'][
                                                  (data['v'] >= individual_standard / data['p']) & (
                                                  data['v'] < institution_standard / data['p']) & (data['t'] > 0)].iloc[
                                                   3:-3].sum()
            trade_type.loc[(d, 'M', 'sell'), stk] = data['v'][
                                                   (data['v'] >= individual_standard / data['p']) & (
                                                   data['v'] < institution_standard / data['p']) & (
                                                   data['t'] < 0)].iloc[
                                                    3:-3].sum()
            trade_type.loc[(d, 'insti', 'buy'), stk] = data['v'][
                                                      (data['v'] >= institution_standard / data['p']) & (
                                                      data['t'] > 0)].iloc[
                                                       3:-3].sum()
            trade_type.loc[(d, 'insti', 'sell'), stk] = data['v'][
                                                       (data['v'] >= institution_standard / data['p']) & (
                                                       data['t'] < 0)].iloc[
                                                        3:-3].sum()
        except KeyError:
            pass
        else:
            pass
    #f.close()
print(time.time() - t0)
trade_type.to_pickle(data_path+'/trade_type_16_100w')
#trade_type.head(20).iloc[:,:5]

#trade_type=pd.read_pickle(data_path+'/trade_type_16_100w')
tmp=trade_type.copy()
tmp.index.names=['trddt','trader','action']
#tmp=tmp.groupby(by='trddt').apply(lambda x:x/x.sum());tmp
#tmp.index.names=['stkcd']
select_stock=DataFrame(select_stock)
select_stock.columns.get_level_values(0)#iloc[:,0].name[0]
Group_By=DataFrame(np.nan,index=tmp.index,columns=select_stock.columns.get_level_values(0)) #
for c in Group_By.columns:
    Group_By[c]=tmp[select_stock[(c,'2015-02-25')]].mean(axis=1)

G1=Group_By.sort_index().loc(axis=0)[:'20150217',:,:];G1
G2=Group_By.sort_index().loc(axis=0)['20150225':,:,:];G2
G1.groupby(by=['trader','action']).mean().loc(axis=0)[:,'buy']
G2.groupby(by=['trader','action']).mean().loc(axis=0)[:,'buy']
G2.groupby(by=['trader','action']).mean().loc(axis=0)[:,'buy']-G1.groupby(by=['trader','action']).mean().loc(axis=0)[:,'buy']
(G2.groupby(by=['trader','action']).mean().loc(axis=0)[:,'buy']-G1.groupby(by=['trader','action']).mean().loc(axis=0)[:,'buy'])/G1.groupby(by=['trader','action']).mean().loc(axis=0)[:,'buy']

G1.groupby(by=['trader','action']).mean().loc(axis=0)[:,'buy'].reset_index(level=1,drop=True)-G1.groupby(by=['trader','action']).mean().loc(axis=0)[:,'sell'].reset_index(level=1,drop=True)
G1.groupby(by=['trader']).head()#apply(lambda x:x['buy']-x['sell'])#.loc(axis=0)[:,'sell']







### 统计交易数据，以区分散户、中户和机构

t0 = time.time()
rootdir = data_path+'/bid_ask'
amount=[]
for d in ['20150121']:
    filename = rootdir + '/' + d
    f = h5py.File(filename, 'r')
    for stk in select_stock.iloc[[17, 16, 10, 15]].values.reshape((60,)):
        try:
            data = DataFrame(
                [list(f['stk'][stk]['volume']),list(f['stk'][stk]['lastPrc'])],index=['v', 'p']).T  #
            data['v'] = data['v'] - data['v'].shift(1)
            amount.extend(list((data['v']*data['p']).values[3:-3]))
        except KeyError:
            pass
        else:
            pass
    #f.close()
print(time.time() - t0) # 21s左右
from scipy.stats import cumfreq
import matplotlib.pyplot as plt
a=np.array(amount);
pd.Series(np.log(a[a>0])+np.log(100)).hist(cumulative=True,bins=200,normed=1)
plt.plot(b)





import h5py
import os
import time

#stocks=pd.read_pickle(data_path+'/selected_stocks_by_size_institution')
t0 = time.time()
rootdir = data_path+'/bid_ask'
li_ = [i for i in os.listdir(rootdir) if not i.endswith('_') and not i.endswith('.h5')][5:]  # 列出文件夹下所有的目录与文件
trade_type = DataFrame(np.nan, index=pd.MultiIndex.from_product([li_, ['indi', 'M', 'insti'], ['buy', 'sell']]),
                       columns=stocks)
institution_standard = 1e4
individual_standard = 1.6e3
for d in li_:
    filename = rootdir + '/' + d
    f = h5py.File(filename, 'r')
    for stk in stocks:
        try:
            data1 = pd.Series(f['stk'][stk]['volume']);
            data2 = pd.Series(f['stk'][stk]['trend']);
            data3 = pd.Series(f['stk'][stk]['lastPrc'])
            data1 = (data1- data1.shift(1))[3:-3]
            trade_type.loc[(d, 'indi', 'buy'), stk] =data1[(data1 < individual_standard / data3[3:-3]) & (data2[3:-3]> 0)].sum()
            trade_type.loc[(d, 'indi', 'sell'), stk] = data1[(data1 < individual_standard / data3[3:-3]) & (data2[3:-3]< 0)].sum()
            trade_type.loc[(d, 'M', 'buy'), stk] = data1[(data1 >= individual_standard / data3[3:-3]) & (data1< institution_standard / data3[3:-3]) & (data2[3:-3]> 0)].sum()
            trade_type.loc[(d, 'M', 'sell'), stk] = data1[(data1>= individual_standard / data3[3:-3]) & (data1 < institution_standard / data3[3:-3]) & (data2[3:-3] < 0)].sum()
            trade_type.loc[(d, 'insti', 'buy'), stk] = data1[(data1 >= institution_standard / data3[3:-3]) & (data2[3:-3] > 0)].sum()
            trade_type.loc[(d, 'insti', 'sell'), stk] = data1[(data1 >= institution_standard / data3[3:-3]) & (data2[3:-3] < 0)].sum()
        except KeyError:
            pass
        else:
            pass
    #f.close()
print(time.time() - t0)  # 从1.5h到105s

trade_type0=pd.read_pickle(data_path+'/trade_type_16_100w')
(trade_type.sort_index()==trade_type0).sum().sum()
trade_type=trade_type.sort_index()



t0=time.time()
data1=pd.Series(f['stk'][stk]['volume']);data2=pd.Series(f['stk'][stk]['trend']);data3=pd.Series(f['stk'][stk]['lastPrc'])
data1 = (data1- data1.shift(1))[3:-3];data1
data1[(data1 < individual_standard / data3[3:-3]) & (data2[3:-3]> 0)].sum()
data1[(data1 < individual_standard / data3[3:-3]) & (data2[3:-3]< 0)].sum()
data1[(data1 >= individual_standard / data3[3:-3]) & (data1< institution_standard / data3[3:-3]) & (data2[3:-3]> 0)].sum()
data1[(data1>= individual_standard / data3[3:-3]) & (data1 < institution_standard / data3[3:-3]) & (data2[3:-3] < 0)].sum()
data1[(data1 >= institution_standard / data3[3:-3]) & (data2[3:-3] > 0)].sum()
data1[(data1 >= institution_standard / data3[3:-3]) & (data2[3:-3] < 0)].sum()
print(time.time()-t0)


