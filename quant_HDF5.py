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
#filename = '/Users/harbes/data/xccdata/bid_ask/20150225.h5'
filename = 'F:/data/xccdata/bid_ask/20150225.h5'

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
import time
import warnings
warnings.filterwarnings("ignore")

#rootdir = '/Users/harbes/data/xccdata/bid_ask'
 rootdir = 'F:/data/xccdata/bid_ask'
li_ = [i for i in os.listdir(rootdir) if not i.endswith('_') and not i.endswith('.h5')] #列出文件夹下所有的目录与文件
n_stock = 40
n_obs=20
n_indi=30
now0=time()
for i in li_[-1:]:  # li_[1:]: #
    path = rootdir+'/'+i
    #if os.path.isfile(path):
    f = h5py.File(path, 'r')
    os.mkdir(rootdir + '/' + i + '_')
    for stock in list(f['stk'])[0:n_stock]:
        data=pd.DataFrame((list(f['stk'][stock][indi][-n_obs:]) for indi in list(f['stk'][stock])[:n_indi]),index=list(f['stk'][stock])[:n_indi])
        data.to_pickle(rootdir + '/' + i +'_'+ '/' +stock)
        #data.to_hdf(rootdir + '/' + i +'.h5',stock,format='table',append=True) # 速度慢，体积大
    f.close()
print(time()-now0)
data = pd.read_pickle('/Users/harbes/data/xccdata/bid_ask/' + li_[-1] + '_/000017')
data.loc[['bidPrc_1', 'bidVol_1', 'bidPrc_2', 'bidVol_2', 'askPrc_1', 'askVol_1', 'askPrc_2', 'askVol_2']]

indi_s = [['bidPrc_1', 'askPrc_1'], ['bidPrc_2', 'askPrc_2'], ['bidPrc_3', 'askPrc_3'], ['bidPrc_4', 'askPrc_4']]
# 运行之前需要删除已存在的
now0 = time()
for i in li_[-1:]:  # li_[1:]: #
    path = rootdir + '/' + i
    f = h5py.File(path, 'r')
    os.mkdir(rootdir + '/' + i + '_')
    data = 2.0 * pd.DataFrame([[((np.array(f['stk'][stock][indi[1]]) - np.array(f['stk'][stock][indi[0]])) / (
    np.array(f['stk'][stock][indi[1]]) + np.array(f['stk'][stock][indi[0]]))).mean()
                                for indi in indi_s] for stock in list(f['stk'])[0:n_stock]],
                              index=list(f['stk'])[0:n_stock], columns=[str(i) for i in range(1, len(indi_s) + 1)])
    data = data[data > 0]
    data.to_pickle(rootdir + '/' + i + '_' + '/bid_ask_mean')
    f.close()
print(time() - now0)

price = pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')['clsprc'].unstack()
group = pd.cut(price.loc['2015-2-17'], bins=[0, 7.5, 20, 30, 40, 1000], labels=range(1, 6))

size = pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')['size_tot'].unstack()
group = pd.qcut(size.loc['2015-2-17'], q=np.linspace(0, 1, 6), labels=range(1, 6))
data.loc[group[group == 1].index.str.slice(0, 6)].median()

price = pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')['clsprc'].unstack()
price.columns = price.columns.str.slice(0, 6)
rootdir = '/Users/harbes/data/xccdata/bid_ask'
# rootdir = 'F:/data/xccdata/bid_ask'
li_ = [i for i in os.listdir(rootdir) if not i.endswith('_') and not i.endswith('.h5')]  # 列出文件夹下所有的目录与文件
n_stock = 4000
n_obs = 20
n_indi = 30
indi_s = ['bidPrc_1', 'askPrc_1']
now0 = time.time()
bid_ask_bias = [[], []]
for i in li_[1:]:  # Mac要额外注意
    path = rootdir + '/' + i
    # if os.path.isfile(path):
    f = h5py.File(path, 'r')
    # os.mkdir(rootdir + '/' + i + '_')
    data = pd.DataFrame(
        [[np.array(f['stk'][stock]['bidPrc_1'][-20:]).mean(), np.array(f['stk'][stock]['askPrc_1'][-20:]).mean()] for
         stock in list(f['stk'])[0:n_stock]],
        index=list(f['stk'])[0:n_stock], columns=indi_s)
    data = data[data['bidPrc_1'] < data['askPrc_1']]
    # data.to_pickle(rootdir + '/' + i +'_'+ '/' +stock)
    # f.close()
    tmp = price.loc[pd.to_datetime(i, format='%Y%m%d')]
    tmp = tmp[tmp < 7.5]
    bid_ask_bias[0].append((tmp[data.index & tmp.index] <= data['bidPrc_1'][data.index & tmp.index]).sum() / len(
        tmp[data.index & tmp.index]))
    bid_ask_bias[1].append((tmp[data.index & tmp.index] >= data['askPrc_1'][data.index & tmp.index]).sum() / len(
        tmp[data.index & tmp.index]))

print(time.time() - now0)

bid_ask_bias = pd.DataFrame(bid_ask_bias, index=['bid', 'ask'], columns=li_[1:]).T
(bid_ask_bias['ask'] - bid_ask_bias['bid']).plot()
bid_ask_bias.plot()
bid_ask_bias

close_ = pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')['adj_close'].unstack()
open_ = pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')['adj_open'].unstack()
rtn = (close_ - open_) / open_ * 100.0
rtn.columns = rtn.columns.str.slice(0, 6)
rtn[data.index & tmp.index].loc['2015-2']

import time

time.localtime(np.array(f['stk'][list(f['stk'])[0]]['localtime'])[0] / 100000)





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

import numpy as np
import pandas as pd
import os
import h5py
import time
import warnings

warnings.filterwarnings("ignore")
import tushare as ts

price = pd.read_pickle('/Users/harbes/data/xccdata/PV_datetime')['clsprc'].unstack()
price.name='close'
price.columns=price.columns.str.slice(0,6)
price.index=price.index.astype(str)

inti_time=time.time()
trade_price=pd.DataFrame(np.nan,index=price.index,columns=price.columns)
trade_type=pd.DataFrame(np.nan,index=price.index,columns=price.columns)
for d in price.index[:10]:
    for stk in price.columns:
        trade_price.loc[d, stk], trade_type.loc[d, stk] = ts.tick(stk, date=d, conn=ts.get_apis()).iloc[0][
            ['price', 'type']]
#trade_price.to_pickle('F:/data/xccdata/trade_price')
#trade_type.to_pickle('F:/data/xccdata/trade_type')
delta_time=time.time()-inti_time

# trade_price.loc[price.index[0],'600848'],trade_type.loc[price.index[0],'600848']=ts.tick('600848', date=price.index[0],conn=ts.get_apis()).loc[0][['price','type']]


ts.get_tick_data(stk, date=d)

tmp = ts.tick('600848', date=price.index[0], conn=ts.get_apis());
tmp
ts.tick()

ts.xpi()

ts.get_markets()













ts.get_tick_data('600848',date='2014-01-09')