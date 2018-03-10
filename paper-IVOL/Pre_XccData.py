# 对xcc提供的数据进行预处理
import pandas as pd

data_path='E:/data/NewData/' #'/Users/harbes/data/NewData/'#
def pre_PV(save_data=False):
    dat= pd.read_pickle(data_path + 'PV')
    dat['trddt'] = pd.to_datetime(dat['trddt'].astype(int).astype(str), format='%Y%m%d')
    dat = dat.set_index(['trddt', 'stkcd']).sort_index()
    dat.drop(dat.index[dat.index.duplicated(keep='last')],inplace=True)
    #con1=dat.duplicated().sum()
    #if con1>0:
    #    raise ValueError('fail to drop duplicates')
    if save_data:
        dat.to_pickle(data_path+'PV_datetime')
    return dat.head()
pre_PV(save_data=True)
def drop_duplicates(data_name,save_data=False):
    dat=pd.read_pickle(data_path+data_name)
    con1=dat.duplicated().sum()
    if con1>0:
        dat.drop_duplicates(inplace=True)
        print('There are duplicates in '+data_name+', and we have droped them, enjoy ^_^ ')
        if save_data:
            dat.to_pickle(data_path+data_name+'drop')
        return dat
    else:
        print('There are no duplicates in '+data_name)


data=drop_duplicates('BS')
