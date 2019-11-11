import pandas as pd
import numpy as np
import scipy as sp
from scipy.sparse.linalg import svds
from pandas.tseries.offsets import DateOffset
DPath='E:/data/CNRDS/'
indicators_all=pd.read_pickle(DPath+'indicators_all')
PV=pd.read_pickle(DPath+'PV');PV.tail(10)
tmp=PV[['Scode','Trddt','Adclsprc']].iloc[1:].set_index(['Trddt','Scode'])
tmp=tmp.astype(float)
adjprc=tmp.loc[~tmp.index.duplicated(keep='last')].unstack() # 提出duplicates
adjprc.index=pd.to_datetime(adjprc.index,format='%Y-%m-%d')# 调整index格式
adjprc.columns=adjprc.columns.get_level_values(1)# 调整columns格式
GroupBy1=lambda x:x.year*100.0+x.month
#p0=adjprc.groupby(by=GroupBy).first();p0
ret=adjprc.groupby(by=GroupBy1).last().pct_change();
ret.index=pd.to_datetime(ret.index.astype(int).astype(str),format='%Y%m')

# 标准化indicators
indicators_all.index.names=['Trddt','Scode']
indi_standardized=indicators_all.astype(float).groupby(by=['Trddt']).apply(lambda x:(x-x.mean())/x.std())
indi_standardized['ret']=ret.stack().drop_duplicates(keep='last')
X_tmp=indi_standardized.iloc[:,:-1].mul(indi_standardized['ret'],axis=0) # 个股characteristics与return的乘积
Nts=(X_tmp.sum(axis=1)!=0).groupby('Trddt').sum()
## TODO
X=X_tmp.astype(float).groupby('Trddt').sum().div(Nts,axis=0) ## 部分数值为0，原因在于原始predictor可能都是nan【为什么predictor会出现大面积且连续的NaN？】;
# 删与不删对结果影响几何？稳健吗？;上述X的列由X_tmp的67变成了53？？？---> 加上astype(float)
indi_standardized=indi_standardized.loc[X_tmp.sum(axis=1)!=0].fillna(0.0)
W=pd.DataFrame(index=pd.MultiIndex.from_product([Nts.index,indi_standardized.columns[:-1]], names=['Trddt','predictors']),columns=indi_standardized.columns[:-1])
for i in Nts.index:
    W.loc[i]=(indi_standardized.loc[i,:'SI'].T@indi_standardized.loc[i,:'SI']/Nts.loc[i]).values


def num_IPCA_estimate_ALS(Gamma_Old,W,X,Nts,PSF=None):
    '''
    % usage:
    % -[ Gamma_New , F_New ] = num_IPCA_estimate_ALS( Gamma_Old , W , X , Nts )
    % -[ Gamma_New , F_New ] = num_IPCA_estimate_ALS( Gamma_Old , W , X , Nts , PSF )
    % inputs
    % - Gamma_Old : LxK matrix of previous iteration's GammaBeta estimate
    % - W : LxLxT array of Z(:,:,t)'*Z(:,:,t)/Nts(t)
    % - X : LxT array of Z(:,:,t)'*Y(:,t)/Nts(t)
    % - Nts : 1xT vector of cross-sectional size [typically found as sum(LOC)]
    % - (optional) PSF : Kadd x T Pre-Specified Factors
    %
    % outputs
    % - Gamma_New : LxK matrix of this iteration's GammaBeta estimate
    % - F_New : KxT matrix of this iteration's Factor estimate
    %
    % Imposes identification assumption on Gamma_New and F_New: Gamma_New is orthonormal
    % matrix and F_New has positive mean (taken across the T columns)
    %
    % When nargin>4, estimates a LxKtilde Gamma_New and a Ktilde x T F_New, where
    % Ktilde=K+Kadd
    % The first K columns of Gamma_New is GammaBeta.
    % The last Kadd columns of Gamma_New is GammaDelta.
    % We continue to only return the KxT factor estimate in F_New
    %
    % When nargin<=4, Ktilde=K and GammaDelta is not estimated
    '''
    T=len(Nts)
    if PSF is not None:
        _,K_add=np.shape(PSF) # 若PSF的dim=1，则此命令出错
        L,K_tilde=np.shape(Gamma_Old)
        K=K_tilde-K_add
    else:
        L, K_tilde = np.shape(Gamma_Old)
        K=K_tilde
    F_New=pd.DataFrame(np.nan,index=Nts.index)
    if PSF is not None:
        for t in Nts.index:
            F_New.loc[t] = np.linalg.pinv(Gamma_Old[:, :K].T @ W.loc[t].values @ Gamma_Old[:, :K]) @ Gamma_Old[:, :K].T @ \
                          (X.loc[t] - W.loc[t].values @ Gamma_Old[:, K:] * PSF.loc[t])
    else:
        for t in Nts.index:
            F_New.loc[t] = np.linalg.pinv(Gamma_Old.T @ W.loc[t].values @ Gamma_Old) @ Gamma_Old.T @ X.loc[t]

    Numerator=np.zeros(L*K_tilde)
    Denominator=np.zeros(L*K_tilde)
    if PSF is not None:
        for i in Nts.index:
            ff=np.vstack((F_New.loc[t].values,PSF.loc[t].values))
            Numerator+=np.kron(X.loc[t],ff)*Nst.loc[t]
            Denominator+=np.kron(W.loc[t].values,ff@ff.T)*Nst.loc[t]



    return None


#Initial guess
K=2
Gamma_Old,s,v     = sp.sparse.linalg.svds(X.T,K);
Factor_Old          = s*v.T
Gamma_Old.shape
Factor_Old.shape