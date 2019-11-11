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
W=W.astype(float)

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
        K_add=np.size(PSF)/len(PSF)
        L,K_tilde=len(Gamma_Old),np.size(Gamma_Old)/len(Gamma_Old)
        K=K_tilde-K_add # K=0表明除了预设的factors之外，并不需要从工具变量中估计其他的factors
    else:
        L, K_tilde = len(Gamma_Old),int(np.size(Gamma_Old)/len(Gamma_Old))
        K=K_tilde # 从工具变量中估计factors的个数
    F_New=pd.DataFrame(np.empty((T,K)),index=Nts.index,columns=list(str(i) for i in range(1,K+1)))
    if K>0:
        if PSF is not None:
            for t in Nts.index:
                F_New.loc[t] = np.linalg.pinv(Gamma_Old[:, :K].T @ W.loc[t].values @ Gamma_Old[:, :K]) @ Gamma_Old[:,
                                                                                                         :K].T @ \
                               (X.loc[t] - W.loc[t].values @ Gamma_Old[:, K:] * PSF.loc[t])
        else:
            for t in Nts.index:
                F_New.loc[t] = np.linalg.pinv(Gamma_Old.T @ W.loc[t].values @ Gamma_Old) @ Gamma_Old.T @ X.loc[t]

    Numerator=np.zeros(L*K_tilde)
    Denominator=np.zeros((L*K_tilde,L*K_tilde))
    if PSF is not None:
        if K>0:
            for i in Nts.index:
                ff = np.vstack((F_New.loc[t].values, PSF.loc[t].values))
                Numerator += np.kron(X.loc[t].values, ff) * Nst.loc[t]
                Denominator += np.kron(W.loc[t].values, ff @ ff.T) * Nst.loc[t]
        else:
            for i in Nts.index:
                ff = PSF.loc[t].values
                Numerator += np.kron(X.loc[t].values, ff) * Nts.loc[t]
                Denominator += np.kron(W.loc[t].values, ff.reshape(K_add, 1) @ ff.reshape(1, K_add)) * Nts.loc[t]
    else:
        for i in Nts.index:
            ff=F_New.loc[t].values
            Numerator+=np.kron(X.loc[t].values,ff)*Nts.loc[t]
            Denominator+=np.kron(W.loc[t].values,ff.reshape(K_tilde,1)@ff.reshape(1,K_tilde))*Nts.loc[t]
    Gamma_New=(np.linalg.pinv(Denominator)@Numerator).reshape(K_tilde,L).T
    # GammaBeta orthonormal and F_New Orthogonal
    if K>0:
        R1 = sp.linalg.cholesky(Gamma_New.T @ Gamma_New, lower=False)  # np.linalg.cholesky给出的是下三角
        R2, _, _ = np.linalg.svd(R1 @ F_New.values.T @ F_New.values @ R1.T)
        Gamma_New[:, :K] = Gamma_New[:, :K] @ np.linalg.pinv(R1) @ R2
        F_New = (np.linalg.pinv(R2) @ R1 @ F_New.T).T
    # Sign convention on GammaBeta and F_New
    if K>0:
        sg=np.sign(F_New.mean(axis=0))
        sg[sg==0]=1.0
        Gamma_New[:,:K]=Gamma_New[:,:K]*sg.values
        F_New*=sg.values
    return Gamma_New,F_New


#Initial guess
K=2
Gamma_Old,s,v     = sp.sparse.linalg.svds(X.T,K);
Factor_Old          = s*v.T
# Numerical choices
MaxIterations       = 100;
Tolerance           = 1e-6;

# Algorithm
for i in range(MaxIterations):
    # Run regressions using Old estimates
    Gamma_New,Factor_New = num_IPCA_estimate_ALS(Gamma_Old,W,X,Nts)
    # Calculate change in Old and New estimates
    tol  = np.maximum(np.abs(Gamma_New-Gamma_Old).max(),np.abs(Factor_New-Factor_Old).max().max()) # other convergence norms could be used
    # Replace Old estimates for the next iteration
    if tol<=Tolerance:
        break
    Factor_Old  = Factor_New
    Gamma_Old   = Gamma_New

Gamma   = Gamma_New;
Factor  = Factor_New;