import pandas as pd
import numpy as np
from pandas.tseries.offsets import YearEnd

insti_holding = pd.read_csv('/Users/harbes/data/xccdata/insti_holding.csv').iloc[1:]
insti_holding['Date'] = pd.to_datetime(insti_holding['Date'])
insti_holding = insti_holding.set_index(['Date', 'SCode'])
insti_holding.index.names = ['trddt', 'stkcd']
insti_holding = insti_holding.sum(axis=1)
# insti_holding=insti_holding.rename('shares')
insti_holding = insti_holding.unstack()
insti_holding.head()

govern_holding = pd.read_csv('/Users/harbes/data/xccdata/govern_holding.csv').iloc[1:]
govern_holding['Year'] = pd.to_datetime(govern_holding['Year']) + YearEnd()
govern_holding = govern_holding.set_index(['Year', 'SCode'])
govern_holding = govern_holding.unstack()
govern_holding.columns = govern_holding.columns.get_level_values(1)
govern_holding.index.names = ['trddt']
govern_holding.columns.names = ['stkcd']
govern_holding.head()
govern_holding[govern_holding.isnull()] = np.nan
govern_holding = govern_holding.astype(np.float64)

non_individual = govern_holding + insti_holding
individual = 100 - non_individual
