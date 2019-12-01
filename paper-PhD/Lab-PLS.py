import pandas as pd
import numpy as np
from pandas.tseries.offsets import DateOffset
DPath='/Users/harbes/PycharmProjects/data/CNRDS/' #'E:/data/CNRDS/'
indicators_all=pd.read_pickle(DPath+'indicators_all').replace([np.inf, -np.inf], np.nan)

PV=pd.read_pickle(DPath+'PV');PV.tail(10)
tmp=PV[['Scode','Trddt','Adclsprc']].iloc[1:].set_index(['Trddt','Scode'])
tmp=tmp.astype(float)
adjprc=tmp.loc[~tmp.index.duplicated(keep='last')].unstack() # 剔除duplicates
adjprc.index=pd.to_datetime(adjprc.index,format='%Y-%m-%d')# 调整index格式
adjprc.columns=adjprc.columns.get_level_values(1)# 调整columns格式
GroupBy1=lambda x:x.year*100.0+x.month
#p0=adjprc.groupby(by=GroupBy).first();p0
ret=adjprc.groupby(by=GroupBy1).last().pct_change();ret
ret.index=pd.to_datetime(ret.index.astype(int).astype(str),format='%Y%m')
# 标准化indicators
indicators_all.index.names=['Trddt','Scode']
indi_standardized=indicators_all.groupby(by=['Trddt']).apply(lambda x:(x-x.mean())/x.std())

# step 1: X=a+b*Z
ret1=ret.stack().loc[indi_standardized.index]
ret1.index.names=['Trddt','Scode']
i_set=set(indi_standardized.index.get_level_values(0)) # 日期
indi_set=indi_standardized.columns
lambda_a=pd.DataFrame(np.nan,index=i_set,columns=indi_set).sort_index()
for i in lambda_a.index:
    lambda_a.loc[i]=indi_standardized.loc[i].corrwith(ret1.loc[i])/ret1.loc[i].std()
lambda_a=lambda_a.astype(float)
#ii=lambda_a.index[6]
#indi_standardized.loc[ii,'SC'] # 估计的lambda_a存在missing value
#indicators_all.loc[ii,'SC']


# step 2: X_{i,t-1(t)}=a+b*lambda_a
B0=pd.Series(np.nan,index=indi_standardized.index).unstack().sort_index();
B1=pd.Series(np.nan,index=indi_standardized.index).unstack().sort_index();
for i in B0.index[:-1]:
    B0.loc[i]=indi_standardized.loc[i].corrwith(lambda_a.loc[i],axis=1)/lambda_a.loc[i].std()
    B1.loc[i+DateOffset(months=1)]=indi_standardized.loc[i+DateOffset(months=1)].corrwith(lambda_a.loc[i],axis=1)/lambda_a.loc[i].std()

# step 3: R_{i,t}=a+b*B
ret1=ret1.unstack()
Gamma=pd.Series(np.nan,index=ret1.index)
for i in ret1.index:
    Gamma.loc[i]=ret1.loc[i].cov(B0.loc[i])/B0.loc[i].var()

AFER=B1.shift(-1)#.mul(Gamma.shift(1),axis=0) # 已经是t+1的indicator了，所以不需要调整位置；单因子设定下，第三步似乎是不需要的

# portfolio analysis
group_num=10
p=np.linspace(0,1,group_num+1)
label_AFER=[i+1 for i in range(group_num)] # 1表示流动性较好的组合，另一个极端组合是流动性较差的组合
mark_AFER=pd.DataFrame([pd.qcut(AFER.loc[i],q=p,labels=label_AFER) for i in AFER.index[1:-1]],index=AFER.index[1:-1]).astype(float)
port_ret=pd.DataFrame([[ret1.loc[i][mark_AFER.loc[i]==j].mean() for j in label_AFER] for i in mark_AFER.index],index=mark_AFER.index,columns=label_AFER)
port_ret.shape
hedge_ret=port_ret.iloc[:,-1]-port_ret.iloc[:,0]

def NWest_mean(d,L=None):
    '''
    df不要以Nan开头，会引起误差;
    :param df:
    :param L:
    :return:
    '''
    df=d.copy()
    df-=df.mean()
    T=len(df)
    if L is None:
        L=int(T**0.25)
    w=1.0-np.arange(1,L+1)/(L+1.0)
    return np.sqrt(2.0*pd.DataFrame((df*df.shift(i+1)*w[i]).sum() for i in range(L)).sum()/T+df.var())/np.sqrt(T)
t_port=port_ret.mean()/NWest_mean(port_ret);t_port
t_hedge=hedge_ret.mean()/NWest_mean(hedge_ret);t_hedge
port_ret.mean()
hedge_ret.mean()

### averaging lambda_a
lambda_aX=lambda_a.rolling(1,min_periods=1).mean();#lambda_aX
# step 2: X_{i,t-1(t)}=a+b*lambda_a
#B0=pd.Series(np.nan,index=indi_standardized.index).unstack().sort_index();
B1=pd.Series(np.nan,index=indi_standardized.index).unstack().sort_index();
for i in B0.index[:-1]:
    #B0.loc[i]=indi_standardized.loc[i].corrwith(lambda_aX.loc[i],axis=1)/lambda_aX.loc[i].std()
    B1.loc[i+DateOffset(months=1)]=indi_standardized.loc[i+DateOffset(months=1)].corrwith(lambda_aX.loc[i],axis=1)/lambda_aX.loc[i].std()

AFER=B1#.mul(Gamma.shift(1),axis=0) # 已经是t+1的indicator了，所以不需要调整位置；单因子设定下，第三步似乎是不需要的

# portfolio analysis
#group_num=10
#p=np.linspace(0,1,group_num+1)
#label_AFER=[i+1 for i in range(group_num)] # 1表示流动性较好的组合，另一个极端组合是流动性较差的组合
mark_AFER=pd.DataFrame([pd.qcut(AFER.loc[i],q=p,labels=label_AFER) for i in AFER.index[1:-1]],index=AFER.index[1:-1]).astype(float)
port_ret=pd.DataFrame([[ret1.loc[i][mark_AFER.loc[i]==j].mean() for j in label_AFER] for i in mark_AFER.index],index=mark_AFER.index,columns=label_AFER)
hedge_ret=port_ret.iloc[:,-1]-port_ret.iloc[:,0]
port_ret.mean()*100
hedge_ret.mean()*100
t_port=port_ret.mean()/NWest_mean(port_ret);t_port
t_hedge=hedge_ret.mean()/NWest_mean(hedge_ret);t_hedge