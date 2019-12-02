import numpy as np
import pandas as pd
import sys
import warnings
warnings.filterwarnings("ignore")
#qcut_options = {'q': np.arange(0.0, 1.001, 0.2)};qcut_options['labels'] = list(range(1, len(qcut_options['q'])))

def GetDataPath():
    sys_platform=sys.platform
    if sys_platform =='win32':
        return 'E:/data/CNRDS/' #
    elif sys_platform=='mac':
        return '/Users/harbes/PycharmProjects/data/CNRDS/' #
    elif sys_platform=='linux':
        return '/home/harbes/data/NewData/'
    else:
        raise ValueError('These is no such systerm in your work-station')
        #return None

def Returns_UnivariatePortfolio(chara,ret,qcut_options,weights=None):
    '''
    Given the aligned pairs (chara,ret), generate the characteristic-based portfolio returns
    :param chara: characteristic ; DataFrame[T*N]
    :param ret: stock returns; DataFrame[T*N] , with date aligned with 'chara'
    :param qcut_options: options setting for pd.qcut; 'q'、'labels'
    :param weights: weights for portfolio construction,if None, then equal weights.
    :return: returns for group portfolios
    '''
    if chara.shape != ret.reindex(chara).shape:
        raise ValueError('chara和ret数据没有对齐！')
    #ValidData = (~chara.isnull()) & (~ret.isnull())
    if weights is not None:
        if chara.shape != weights.shape:
            raise ValueError('weights与其他数据没有对齐！')
        #else:
        #    ValidData= ValidData & (~weights.isnull())
        #    weights=weights[ValidData]
    #chara = chara[ValidData]
    #ret = ret[ValidData]
    mark_group=pd.DataFrame(np.nan,index=chara.index,columns=chara.columns)
    for t in chara.index:
        mark_group.loc[t]=pd.qcut(chara.loc[t],q=qcut_options['q'],labels=qcut_options['labels'])
    ret_portfolio=pd.DataFrame(np.nan,index=chara.index,columns=qcut_options['labels'])
    if weights is None:
        for c in qcut_options['labels']:
            ret_portfolio.loc[:,c]=ret[mark_group==c].mean(axis=1)
    else:
        agg_weights=pd.DataFrame(np.nan,index=chara.index,columns=qcut_options['labels'])
        ret_weights = ret * weights
        weights_ = (~ret.isnull()) * weights
        for c in qcut_options['labels']:
            ret_portfolio.loc[:,c]=ret_weights[mark_group==c].sum(axis=1)
            agg_weights.loc[:,c]=weights_[mark_group==c].sum(axis=1)
        ret_portfolio=ret_portfolio/agg_weights
    return ret_portfolio.astype(float)

def Returns_BivariatePortfolio(chara,chara_con,ret,qcut_options,weights=None,independent=False):
    '''

    :param chara:
    :param chara_con: conditional characteristic ; DataFrame[T*N]
    :param ret:
    :param qcut_options: 'q','q_con','labels','labels_con'
    :param weights:
    :return:
    '''

    if (chara.shape != ret.reindex(chara).shape) or (chara_con.shape != ret.reindex(chara).shape):
        raise ValueError('chara、chara_con和ret数据没有对齐！')
    if (weights is not None) & (chara.shape != weights.shape):
        raise ValueError('weights与其他数据没有对齐！')
    mark_group=pd.DataFrame(np.nan,index=chara.index,columns=chara.columns)
    mark_group_con=pd.DataFrame(np.nan,index=chara_con.index,columns=chara_con.columns)
    if independent:
        for t in chara.index:
            mark_group.loc[t]=pd.qcut(chara.loc[t],q=qcut_options['q'],labels=qcut_options['labels'])
            mark_group_con.loc[t]=pd.qcut(chara_con.loc[t],q=qcut_options['q_con'],labels=qcut_options['labels_con'])
    else:
        for t in chara_con.index:
            mark_group_con.loc[t] = pd.qcut(chara_con.loc[t], q=qcut_options['q_con'],
                                            labels=qcut_options['labels_con'])
            for c in qcut_options['labels_con']:
                mark_group.loc[t]=mark_group.loc[t].combine_first(pd.qcut(chara.loc[t,mark_group_con.loc[t]==c],
                                                                                  q=qcut_options['q'],
                                                                                  labels=qcut_options['labels']))
    ret_portfolio = pd.DataFrame(np.nan, index=pd.MultiIndex.from_product(chara.index,qcut_options['labels_con']),
                                 columns=qcut_options['labels'])
    if weights is None:
        for i in qcut_options['labels_con']:
            for c in qcut_options['labels']:
                ret_portfolio.loc[(slice(None),i),c]=ret[(mark_group_con==i) & (mark_group==c)].mean(axis=1)
    else:
        agg_weights = pd.DataFrame(np.nan, index=pd.MultiIndex.from_product(chara.index, qcut_options['labels_con']),
                                     columns=qcut_options['labels'])
        ret_weights=ret*weights
        weights_=(~ret.isnull()) * weights
        for i in qcut_options['labels_con']:
            for c in qcut_options['labels']:
                ret_portfolio.loc[(slice(None),i),c]=ret_weights[(mark_group_con==i) & (mark_group==c)].sum(axis=1)
                agg_weights.loc[(slice(None),i),c]=weights_[(mark_group_con==i) & (mark_group==c)].sum(axis=1)
        ret_portfolio=ret_portfolio/agg_weights
    return ret_portfolio.astype(float)

def describe(df,stats=['skew','kurt']):
    d=df.describe(percentiles=[0.05,0.25,0.5,0.75,0.95])
    return d.append(df.reindex(d.columns, axis=1).agg(stats))
def NWest_mean(dataframe,L=None):
    '''
    df不要以Nan开头，会引起误差;
    求时序均值的标准差估计
    :param df:
    :param L:
    :return:
    '''
    df=dataframe-dataframe.mean()
    T=len(df)
    if L is None:
        L=int(T**0.25)
    w=1.0-np.arange(1,L+1)/(L+1.0)
    return np.sqrt(2.0*pd.DataFrame((df*df.shift(i+1)*w[i]).sum() for i in range(L)).sum()/T+df.var())/np.sqrt(T)
def NWest(e,X,L=None):
    T = len(e)
    if L is None:
        L = int(T ** 0.25) # or : L = 0.75*T**(1.0/3.0)-1
    w = 1.0 - np.arange(1, L + 1) / (L+1.0)
    X.insert(0,'c',np.ones(T))
    S=0.0
    for l in range(1,L+1):
        for i in range(l,T):
            S+=w[l-1]*e[i]*e[i-l]*(X.iloc[i][:,None]*X.iloc[i-l].values+X.iloc[i-l][:,None]*X.iloc[i].values)
    for i in range(T):
        S+=e[i]*e[i]*X.iloc[i][:,None]*X.iloc[i].values
    XX_1=np.linalg.pinv(X.T@X.values)
    X.drop('c', axis=1, inplace=True)
    return np.sqrt((XX_1@S@XX_1)[0,0])

BM.head()
