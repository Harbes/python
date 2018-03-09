# used to accomplish the IVOL in 《empirical asset pricing》
# 1. calculate the market beta
# 2. calculate the returns of two mimicking portfolios--- size & BM

import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt

data_path='/Users/harbes/data/NewData/' #'E:/data/NewData/'
pv=pd.read_pickle(data_path+'PV_datetime')[['adj_close','adj_open']]
close_price=pv['adj_close'].unstack()
#close_price.index=pd.to_datetime(close_price.index.astype(int).astype(str),format='%Y%m%d')
open_price=pv['adj_open'].unstack()
#open_price.index=pd.to_datetime(open_price.index.astype(int).astype(str),format='%Y%m%d')
index_ret=pd.read_pickle(data_path+'index_ret').set_index(['index_code','trddt'])['pctchange'].loc['000016.SH']
index_ret.index=pd.to_datetime(index_ret.index.astype(int).astype(str),format='%Y%m%d')
def cal_beta():
    pass

rtn=(close_price-open_price)/open_price*100
rtn['index']=index_ret['pctchange']
rtn['month']=(rtn.index.year-rtn.index[0].year)*12+rtn.index.month
index_ret=pd.DataFrame({'pctchange':index_ret,'month':(index_ret.index.year-index_ret.index[0].year)*12+index_ret.index.month})
beta=pd.DataFrame(index=pd.date_range('20050101','20180301',freq='M'),columns=rtn.columns)

from time import time
t0=time()
for i in range(3,159):
    beta.iloc[i-1]=rtn.iloc[:,:-1][(i-2<=rtn['month']) & (rtn['month']<=i)].cov(min_periods=20).iloc[-1,:-1]/\
                 np.var(index_ret['pctchange'][(index_ret['month']>=i-2)&(index_ret['month']<=i)])
print(time()-t0)
beta.to_pickle(data_path+'beta_daily_3M')



