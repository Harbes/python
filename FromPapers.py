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
    res=DoublePortAnalysis(['turnover'],'ivol',var_dict=var_dict)
    FM_res=Fama_MacBeth(var_list,freq=freq,var_dict=var_dict)
    return res,FM_res



if __name__ is '__main__':
    t0=time()
    res_2,FM_res=SizeAndValueInChina_LiuStambaughYuan()
    print(time()-t0)