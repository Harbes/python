from modules.quant.QuantBase import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def cal_variables(var_list,start='20050101',end='20180228'):
    #TODO
    pass
def cal_weekday_ret():
    # TODO
    pass

if __name__ == '__main__':
    freq = 'M'
    ret=cal_ret(freq=freq)
    size=cal_size(freq=freq)
    port_size=cal_mimick_port1(size,freq=freq,ret=ret)
    cum_ret=CumRet(port_size)
    cum_ret.plot()
    plt.show()