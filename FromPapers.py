from my_modules.QuantTools.FactorBase2 import *
import numpy as np
import pandas as pd
from time import time
#from importlib import reload
#import matplotlib.pyplot as plt

def SizeAndValueInChina_LiuStambaughYuan(clip=True):
    '''
    《Size and Value in China》-Jianan Liu, Robert F. Stambaugh, and Yu Yuan
    excule the smallest 30% of firms
    :return:
    '''
    freq='M'
    var_list = ['beta','size','BM','mom','rev','max_ret', 'turnover','skew', 'ivol']
    var_dict=GetVarsFromList(var_list,freq)
    var_dict['ret']=ClipQuantile(var_dict['size'],percentile=[0.0,0.3,1.0],label=[-1.0,1.0])
    res1_EW=SinglePortAnalysis(['ivol'],var_dict=var_dict,freq=freq)
    res1_VW = SinglePortAnalysis(['ivol'], var_dict=var_dict, freq=freq,value_weighted=True)
    res2_EW=DoublePortAnalysis(['turnover'],'ivol',var_dict=var_dict)
    res2_VW = DoublePortAnalysis(['turnover'], 'ivol', var_dict=var_dict,value_weighted=True)
    FM_res=Fama_MacBeth(var_list,freq=freq,var_dict=var_dict)
    return res1_EW,res1_VW,res2_EW,res2_VW,FM_res,var_dict

if __name__ is '__main__':
    #t0=time()
    #res1_EW, res1_VW, res2_EW, res2_VW, FM_res,var_dict=SizeAndValueInChina_LiuStambaughYuan()
    #print(time()-t0) # 450s


    t0 = time()
    freq = 'M'
    var_list = ['beta','size','BM','mom','rev','illiq', 'turnover','skew', 'coskew','iskew','ivol','max_ret','e_ivol']#['beta','size','BM','ivol','turnover','illiq','rev','skew']#
    var_dict = GetVarsFromList(var_list, freq)
    var_d = winsorize_vars(var_dict,(0.01,0.01),var_list=var_list)
    var_d['shell']=var_d['ivol']*pd.DataFrame([pd.qcut(var_d['ivol'].loc[i],q=[0.0,0.65,1.0],labels=[0.0,1.0]) for i in var_d['ivol'].index ])
    var_l = ['turnover','ivol','beta','size','BM','rev','illiq','iskew','shell']#
    # 需要调整cal_SMB_HML的默认参数
    index_ret, SMB, HML = Get_Index_SMB_HML(freq, var_dict)
    #res1_EW = SinglePortAnalysis(var_list, var_dict=var_d, index_ret=index_ret, SMB=SMB, HML=HML)
    #res1_VW = SinglePortAnalysis(var_list, var_dict=var_d, index_ret=index_ret, SMB=SMB, HML=HML,value_weighted=True)
    #mark_ = ClipQuantile(var_d['size'], (0.0, 0.3, 1.0), (-1.0, 1.0), output_indi=False)
    #var_d['ret']=var_d['ret'][mark_]
    res2_EW = DoublePortAnalysis(var_list, 'ivol', var_dict=var_d, index_ret=index_ret, SMB=SMB, HML=HML)
    res2_VW = DoublePortAnalysis(var_list, 'ivol', var_dict=var_d, index_ret=index_ret, SMB=SMB, HML=HML,value_weighted=True)
    FM_res=Fama_MacBeth(var_l,freq=freq,var_dict=var_d);FM_res
    print(time() - t0)
    #var_dict['e_ivol']=pd.read_pickle(GetDataPath()+'Expected_iVol_using_OptimalEgarch').shift(1).loc['2007-09':]
    t_beta=np.empty(12)
    index_ret_d=cal_index_ret(freq='D')
    ret_d=cal_ret(freq='D')
    ret=cal_ret(freq='M')
    for i in range(1,13):
        beta=cal_beta(i,freq='M',index_ret_d=index_ret_d,ret_d=ret_d)
        tmp=cal_mimick_port1(beta,freq='M',ret=ret)
        tmp[len(tmp.columns)+1]=tmp[len(tmp.columns)]-tmp[1.0]
        t_beta[i-1]=(tmp.mean()/NWest_mean(tmp)).iloc[-1]

    res2_EW.to_csv(GetDataPath()+'res2_EW.csv')
    res2_VW.to_csv(GetDataPath() + 'res2_VW.csv')