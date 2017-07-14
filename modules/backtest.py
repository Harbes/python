#from my_modules import backtest
#cumreturns,SharpeRatio,drawdown=backtest.backtest(data_rtn_group)





import numpy as np
import pandas as pd
def backtest(returns,rtype='simple',benchmark=None,freq='d'):
    '''
    
    :param returns: array
        the stock returns
    :param rtype: string
        "simple" or "compound"
    :return: 
        the cumulative returns
        the corresponding maximum drawback 
        the maximum drawback duration
    '''
#TODO 最好限定输入的 returns 为 Series 或者 DataFrame ，不然整个 backtest 的格式调整会非常麻烦
    if np.size(np.shape(returns))==1:
        if type(returns) != pd.core.frame.Series:
            returns=pd.Series(returns)
        drawdown_duration = pd.Series(np.zeros_like(returns),index=returns.index)
    else:
        if type(returns) !=pd.core.frame.DataFrame:
            returns=pd.DataFrame(returns)
        drawdown_duration = pd.DataFrame(np.zeros_like(returns),index=returns.index,columns=returns.columns)
    if benchmark==None:
        if freq=='d':
            SharpeRatio=np.sqrt(252)*np.mean(returns,axis=0)/np.std(returns,axis=0)
        if freq=='m':
            SharpeRatio = np.sqrt(12) * np.mean(returns, axis=0) / np.std(returns, axis=0)
    elif benchmark != None:
        if len(benchmark)!=len(returns):
            raise ValueError('the benchmark does not match the asset returns')
        else:
            if type(benchmark)!= np.core.frame.Series
                benchmark=pd.Series(benchmark,index=returns.index)
            if freq == 'd':
                SharpeRatio = np.sqrt(252) * np.mean(returns-benchmark, axis=0) / np.std(returns-benchmark, axis=0)
            if freq == 'm':
                SharpeRatio = np.sqrt(12) * np.mean(returns-benchmark, axis=0) / np.std(returns-benchmark, axis=0)


    if rtype == 'simple':
        cumreturns=(returns+1).cumprod(axis=0)-1
    high_watermark=cumreturns.cummax()
    drawdown=1-(1+cumreturns)/(1+high_watermark)
#
#    for j in range(np.shape(drawdown)[1]):
#        for t in range(np.shape(drawdown)[0]-1):
#            if drawdown.ix[t+1,j]==0:
#                drawdown_duration.ix[t+1,j]=0
#            else:
#                drawdown_duration.ix[t+1,j]=drawdown_duration.ix[t,j]+1

    drawdown=drawdown.max(axis=0)*100
#    drawdown_duration=drawdown_duration.max(axis=0)
    return cumreturns,SharpeRatio,drawdown