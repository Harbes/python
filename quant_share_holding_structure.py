import pandas as pd
from pandas.tseries.offsets import YearEnd

insti_holding = pd.read_csv('/Users/harbes/data/xccdata/insti_holding.csv').iloc[1:].set_index(['Date', 'SCode'])
insti_holding.index.names = ['trddt', 'stkcd']
insti_holding.head()

govern_holding = pd.read_csv('/Users/harbes/data/xccdata/govern_holding.csv').iloc[1:]
govern_holding['Year'] = govern_holding['Year'].to_datetime() + YearEnd()
govern_holding.set_index(['Year', ''])
govern_holding.head()

govern_holding.index.get_level_values(0) = govern_holding.index.get_level_values(0).to_datetime() + YearEnd()
