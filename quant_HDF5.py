import numpy as np
import pandas as pd
data=pd.read_hdf('F:/data/xccdata/20150105')
# 报错：ValueError: No dataset in HDF5 file.
data=pd.read_hdf('F:/data/xccdata/20150105',key='stk',columns='000001')
# 报错：TypeError: cannot create a storer if the object is not existing nor a value are passed


store=pd.HDFStore('F:/data/xccdata/20150105')
store.keys()
# 结果是：[]
store.root
# / (RootGroup) ''
#  children := ['stk' (Group)]
store.root.stk
store.root['stk']
# 报错:TypeError: 'RootGroup' object is not subscriptable
store.close()


import h5py
filename = 'F:/data/xccdata/20150105'
f = h5py.File(filename, 'r')
a_group_key = list(f.keys())[0];a_group_key
list(f[a_group_key]['002666'])
data = pd.DataFrame(f[a_group_key]['002666'][['volume',
 'askVol_4',
 'askPrc_2',
 'bidPrc_1',
 'trend',
 'numTrades',
 'askPrc_4',
 'bidVol_2',
 'bidVol_3',
 'turnover',
 'bidVol_1',
 'bidVol_5',
 'localtime',
 'bidVol_4',
 'askPrc_1',
 'askPrc_3',
 'askVol_2',
 'askVol_3',
 'askVol_5',
 'askPrc_5',
 'bidPrc_4',
 'bidPrc_5',
 'bidPrc_3',
 'codeint',
 'servertime',
 'bidPrc_2',
 'askVol_1',
 'lastPrc']
])#['bidPrc_1'])#
data
data=pd.DataFrame(f[a_group_key])
# 报错：ValueError: DataFrame constructor not properly called!
f.close()
