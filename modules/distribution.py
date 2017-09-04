
import numpy as np
import scipy.stats as scs
import matplotlib.pyplot as plt
import math
class distribution:
    def __init__(self,arr):
        if type(arr) is not np.ndarray:
            self.arr=np.array(arr)
        else:
            self.arr=arr
    def decile(self,q=[10,20,30,40,50,60,70,80,90]):
        return np.percentile(self.arr,q)

    def normal_test(self):
        '''
        test for normality distribution of given data set.
        :param arr: ndarray,Series,...
            object to generate statistics on
        :return:
        '''
        print('Null hypothesis: it is normally distributed')
        print('======')
        print('     Skew of data set %14.3f' % scs.skew(self.arr))
        print('     Skew test p-value %13.3f' % scs.skewtest(self.arr)[1])
        print('======')
        print('     Kurt of data set %14.3f' % scs.kurtosis(self.arr))
        print('     Kurt test p-value %13.3f' % scs.kurtosistest(self.arr)[1])
        print('======')
        print('     Norm test p-value %13.3f' % scs.normaltest(self.arr)[1])
    def hist(self,bins=40,normed=1,disp_norm=0):
        plt.figure()
        plt.hist(self.arr, bins=bins, normed=normed)
        if disp_norm:
            lb=np.min(self.arr)
            ub=np.max(self.arr)
            plt.plot(np.arange(lb,ub+0.001,0.01),scs.norm(np.mean(self.arr),np.std(self.arr)).pdf(np.arange(lb,ub+0.001,0.01)))
        plt.show()
def NormalCdf(d):
    a1=0.319381530
    a2=-0.356563782
    a3=1.781477937
    a4=-1.821255978
    a5=1.330274429
    gamma=0.2316419
    normalprime=math.exp(-d*d/2.0)/math.sqrt(2*math.pi)
    if d>0:
        k1 = 1.0 / (1.0 + gamma * d)
        return 1-normalprime*(a1*k1+a2*k1**2+a3*k1**3+a4*k1**4+a5*k1**5)
    else:
        k1 = 1.0 / (1.0 - gamma * d)
        return normalprime * (a1 * k1 + a2 * k1 ** 2 + a3 * k1 ** 3 + a4 * k1 ** 4 + a5 * k1 ** 5)