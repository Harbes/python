import pandas as pd
from pandas import DataFrame
import numpy as np
from pandas.tseries.offsets import YearEnd
from datetime import datetime
import h5py
import os
import time
data_path = '/Users/harbes/data/xccdata'
# data_path='F:/data/xccdata'
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
    for stk in f['stk']:
        try:
            data1 = pd.Series(f['stk'][stk]['volume']);
            data2 = pd.Series(f['stk'][stk]['trend']);
            data3 = pd.Series(f['stk'][stk]['lastPrc'])
            data1 = data1.diff(1)[3:-3]
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
