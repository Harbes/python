from modules.quant.FactorBase import *
import numpy as np
import pandas as pd
from time import time
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
    var_list = ['beta','size','BM','mom','rev','illiq','max_ret', 'turnover','skew', 'coskew','vol','ivol']
    var_dict = GetVarsFromList(var_list, freq)
    var_d = winsorize_vars(var_dict,(0.01,0.01),var_list=var_list)
    index_ret, SMB, HML = Get_Index_SMB_HML(freq, var_dict)

    res1_EW = SinglePortAnalysis(var_list, var_dict=var_d, index_ret=index_ret, SMB=SMB, HML=HML)
    res1_VW = SinglePortAnalysis(var_list, var_dict=var_d, index_ret=index_ret, SMB=SMB, HML=HML,value_weighted=True)
    mark_ = ClipQuantile(var_d['size'], (0.0, 0.3, 1.0), (-1.0, 1.0), output_indi=False)
    var_d['ret']=var_d['ret'][mark_]
    res2_EW = DoublePortAnalysis(var_list, 'ivol', var_dict=var_d, index_ret=index_ret, SMB=SMB, HML=HML)
    res2_VW = DoublePortAnalysis(var_list, 'ivol', var_dict=var_d, index_ret=index_ret, SMB=SMB, HML=HML,value_weighted=True)
    var_l=['beta','size','BM','rev','illiq', 'turnover','ivol','skew']
    FM_res=Fama_MacBeth(var_l,freq=freq,var_dict=var_d);FM_res
    print(time() - t0)