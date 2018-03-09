# 对xcc提供的数据进行预处理
import pandas as pd

data_path='/Users/harbes/data/NewData/'# 'E:/data/NewData/'
def pre_PV(save_data=False):
    dat=pd.read_pickle(data_path+'PV')
    dat= pd.read_pickle(data_path + 'PV')
    dat.drop_duplicates(inplace=True)
    con1=dat.duplicated().sum()
    if con1>0:
        raise ValueError('fail to drop duplicates')
    if save_data:
        dat['trddt']=pd.to_datetime(dat['trddt'].astype(int).astype(str),format='%Y%m%d')
        dat = dat.set_index(['trddt', 'stkcd']).sort_index()
        dat.to_pickle(data_path+'PV_datetime')
    return dat.head()

pre_PV(save_data=True)
