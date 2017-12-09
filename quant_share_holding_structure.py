import pandas as pd
from pandas import DataFrame
import numpy as np
from pandas.tseries.offsets import YearEnd
from datetime import datetime

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

institution = pd.read_csv('/Users/harbes/data/xccdata/insti_holding.csv').iloc[1:]
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

government = pd.read_csv('/Users/harbes/data/xccdata/govern_holding.csv').iloc[1:]
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

PV_datetime = pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')[['size_tot', 'adj_close', 'adj_open', 'clsprc']]
filter_ = pd.read_pickle('/Users/harbes/data/xccdata/filter')
size = PV_datetime['size_tot'].unstack()
size.columns = size.columns.str.slice(0, 6)
adj_open = PV_datetime['adj_open'].unstack()
adj_open.columns = adj_open.columns.str.slice(0, 6)
adj_close = PV_datetime['adj_close'].unstack()
adj_close.columns = adj_close.columns.str.slice(0, 6)
clsprc = PV_datetime['clsprc'].unstack()
clsprc.columns = clsprc.columns.str.slice(0, 6)

book = pd.read_pickle('/Users/harbes/data/xccdata/BS')[['ann_dt', 'stkcd', 'tot_assets', 'tot_liab']]  #
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
