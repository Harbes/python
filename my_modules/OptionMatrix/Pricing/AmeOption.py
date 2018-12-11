#import numpy as np
from scipy.stats import norm
from math import log,exp,sqrt

def cnd(x,method='default'): # 正态分布累计分布函数
    if method=='default':
        return norm.cdf(x)