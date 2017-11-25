import numpy as np
import pandas as pd
data=pd.read_hdf('F:/data/xccdata/20150105')
# 报错：ValueError: No dataset in HDF5 file.
data=pd.read_hdf('F:/data/xccdata/20150105',key='stk',columns='000001')
# 报错：TypeError: cannot create a storer if the object is not existing nor a value are passed


store=pd.HDFStore('/Users/harbes/data/xccdata/bid_ask/20150225')
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
filename = '/Users/harbes/data/xccdata/bid_ask/20150225'
f = h5py.File(filename, 'r')
#f.name
#list(f.keys())
a_group_key = list(f.keys())[0] ;a_group_key
stkcd=list(f[a_group_key])
stkcd[:10]
data = pd.DataFrame(list(f['stk']['002666']['servertime']),columns=['servertime'])#
data
data=pd.DataFrame(f[a_group_key])
# 报错：ValueError: DataFrame constructor not properly called!
f.close()


import numpy as np
import pandas as pd
import os
import h5py
from time import time

rootdir = '/Users/harbes/data/xccdata/bid_ask'
# rootdir = 'F:/data/xccdata/bid_ask'
li_ = os.listdir(rootdir) #列出文件夹下所有的目录与文件
n_stock=4000
n_obs=20
n_indi=30
now0=time()
for i in li_: # li_[1:]: #
    path = os.path.join(rootdir,i)
    #if os.path.isfile(path):
    f = h5py.File(path, 'r')
    os.mkdir(rootdir + '/' + i + '_')
    for stock in list(f['stk'])[0:n_stock]:
        data=pd.DataFrame((list(f['stk'][stock][indi][-n_obs:]) for indi in list(f['stk'][stock])[:n_indi]),index=list(f['stk'][stock])[:n_indi])
        data.to_pickle(rootdir + '/' + i +'_'+ '/' +stock)
    f.close()
print(time()-now0)


indi_=['volume',
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




























