import pandas as pd
import numpy as np

data_path='E:/data/NewData/'
pv=pd.read_pickle(data_path+'PV_indexed')
index_ret=pd.read_pickle(data_path+'index_ret').set_index(['index_code','trddt'])['pctchange'].loc['000016.SH']

