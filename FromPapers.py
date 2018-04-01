from modules.quant.FactorBase import *
import numpy as np
import pandas as pd
from time import time
#import matplotlib.pyplot as plt

def SizeAndValueInChina_LiuStambaughYuan():
    '''
    《Size and Value in China》-Jianan Liu, Robert F. Stambaugh, and Yu Yuan
    excule the smallest 30% of firms
    :return:
    '''
    freq='M'
    var_list = ['beta','size','BM','mom','rev','max_ret', 'turnover','skew', 'ivol']
    var_dict=GetVarsFromList(var_list,freq)
    filter_=FilterSample(var_dict['size'],percentile=[0.0,0.3,1.0],label=[0.0,1.0])
    var_dict['ret']=var_dict['ret'][filter_]
    res1_EW=SinglePortAnalysis(['ivol'],var_dict=var_dict,freq=freq)
    res1_VW = SinglePortAnalysis(['ivol'], var_dict=var_dict, freq=freq,value_weighted=True)
    res2_EW=DoublePortAnalysis(['turnover'],'ivol',var_dict=var_dict)
    res2_VW = DoublePortAnalysis(['turnover'], 'ivol', var_dict=var_dict,value_weighted=True)
    FM_res=Fama_MacBeth(var_list,freq=freq,var_dict=var_dict)
    return res1_EW,res1_VW,res2_EW,res2_VW,FM_res,var_dict


if __name__ is '__main__':
    t0=time()
    res1_EW, res1_VW, res2_EW, res2_VW, FM_res,var_dict=SizeAndValueInChina_LiuStambaughYuan()
    print(time()-t0) # 450s