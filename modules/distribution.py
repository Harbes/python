
import numpy as np
import scipy.stats as scs
import matplotlib.pyplot as plt
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
    def hist(self,bins=40,normed=1):
        fig,ax=plt.subplots()
        ax.hist(self.arr,bins=bins,normed=normed)
        plt.show()