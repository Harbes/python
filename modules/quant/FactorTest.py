from modules.quant.FactorBase import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


if __name__ == '__main__':
    #freq = 'M'
    #ret=cal_ret(freq=freq)
    #size=cal_size(freq=freq)
    #port_size=cal_mimick_port1(size,freq=freq,ret=ret)
    #cum_ret=CumRet(port_size)
    #cum_ret.plot()
    #plt.show()
    #np.percentile()

    #from time import time
    #t0 = time()
    #var_list = ['ivol', 'beta', 'size', 'BM', 'mom', 'rev', '']
    #a, b = Fama_MacBeth(var_list)
    #print(time() - t0)

    #mom = cal_mom(12, freq='M')
    #tmp = cal_mimick_port1(mom, freq='M', ret=ret.loc[mom.index])
    #tmp[len(tmp.columns) + 1.0] = tmp.iloc[:, -1] - tmp.iloc[:, 0]
    #tmp.mean() / NWest_mean(tmp)
    #cal_FF_alpha(tmp)


    # portfolio analysis

    from time import time

    t0 = time()
    # 当包含'coskew'时,ValueError: cannot reindex from a duplicate axis
    # 当包含'iskew'时,ValueError: Shape of passed values is (3500, 68564), indices imply (3500, 6286)
    var_list =['beta', 'size', 'mom', 'rev', 'illiq', 'turnover', 'max_ret', 'skew', 'vol','ivol']  #['iskew']#
    var_dict=GetVarsFromList(var_list,freq='M')
    print(time() - t0)
    results1_EW = SinglePortAnalysis(var_list,var_dict=var_dict)
    print(time() - t0)
    results1_VW = SinglePortAnalysis(var_list,var_dict=var_dict,value_weighted=True)
    print(time() - t0)
    results2_EW = DoublePortAnalysis(var_list,'ivol',var_dict=var_dict)
    print(time() - t0)
    results2_VW = DoublePortAnalysis(var_list,'ivol',var_dict=var_dict,value_weighted=True)
    print(time() - t0)