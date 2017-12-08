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
institution.head()
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
government.head()
# (government==100).sum(axis=1) # 完全由政府持股

non_individual = government + institution
individual = 100 - non_individual
# (individual==100).sum(axis=1) # 完全由个人持股

size = pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')['size_tot'].unstack()
size.columns = size.columns.str.slice(0, 6)

institution = institution.resample('D').ffill()
institution = institution.loc[size.index & institution.index]
institution.tail()

individual = individual.resample('D').ffill()
individual = individual.loc[size.index & individual.index]
individual.tail()

# size-insti sort (evenly sort)
num_by_ = 5;
label_size = [i + 1 for i in range(num_by_)];
percentile_size = np.linspace(0, 1, num_by_ + 1)
label_2 = [1, 2, 3];
percentile_2 = (0.0, 0.3, 0.7, 1.0)
indicator1 = size.shift(1).loc['2006':'2016'][size.columns & institution.columns]
indicator2 = institution.loc['2006':'2016'][size.columns & institution.columns]
mark_1 = DataFrame([pd.qcut(indicator1.loc[i], q=percentile_size, labels=label_size) for i in
                    indicator1.index])  # ,index=indicator1.index,columns=indicator1.columns)
mark_2 = DataFrame(np.nan, index=mark_1.index, columns=mark_1.columns)
# mark_2[institution==0]=1
for l_ in label_size:
    tmp = DataFrame(
        [pd.qcut(indicator2.loc[i][mark_1.loc[i] == l_], q=percentile_2, labels=label_2) for i in indicator2.index],
        index=indicator2.index)
    mark_2 = mark_2.combine_first(tmp)

institution = institution[size.columns & institution.columns]
size = size[size.columns & institution.columns]
average_institution_holding = DataFrame(np.zeros((10, 11)), index=pd.MultiIndex.from_product([label_size, (1, 3)]),
                                        columns=np.array(after_fes_data)[1:-1, 0])
average_size = DataFrame(np.zeros((10, 11)), index=pd.MultiIndex.from_product([label_size, (1, 3)]),
                         columns=np.array(after_fes_data)[1:-1, 0])
select_stock = {}
for i in label_size:
    for j in (1, 3):
        for d in np.array(after_fes_data)[1:-1, 0]:
            average_institution_holding.loc[(i, j), d] = institution.loc[d][
                (mark_1.loc[d] == i) & (mark_2.loc[d] == j)].mean()
            average_size.loc[(i, j), d] = size.loc[d][(mark_1.loc[d] == i) & (mark_2.loc[d] == j)].mean()
            select_stock[((i, j), d)] = np.sort(
                np.random.choice(size.loc[d][(mark_1.loc[d] == i) & (mark_2.loc[d])].index, 25))

average_institution_holding.mean(axis=1)
average_size.mean(axis=1)

np.random.choice(size.loc[d][(mark_1.loc[d] == i) & (mark_2.loc[d])].index,25)
DataFrame(select_stock)
