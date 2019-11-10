# from quantpy import Portfolio as Portf
import pandas as pd
import numpy as np
from pylab import legend, xlabel, ylabel, \
    cov, sqrt, mean, std, plot, show, figure
from numpy import array, zeros, matrix, ones, linspace, hstack
from pandas import Series, DataFrame
from numpy.linalg import inv,pinv
import matplotlib.pyplot as plt


class Portfolio():
    def __init__(self, sample_num=20, select_type='s', start='20070104', end='20151231', bench='000016.SH', \
                 data_source='/Users/harbes/data/xccdata'):
        """
        Parameter
        -------
        sample_num:int
            the number of stocks selected
        select_type:string
            'a':asending; or 'h':head
            'd':descending; or 't':tail
            others:stochastic
        start,end:datetime,string or None
        bench:string
            the benchmark index,like '000016.SH'
        data_source；string
            local data dir

        """

        # Create distionary to hold assets.
        self.asset = {}
        # transform 'start' and 'end' to datetime
        if type(start) is str and type(end) is str:
            start = pd.to_datetime(start)
            end = pd.to_datetime(end)

        # Get Benchmark asset.
        temp = pd.read_pickle(data_source + '/index_ret_datetime') \
            [['index_code', 'trddt', 'clsprc']].set_index(['trddt', 'index_code'])
        self.benchmark = temp.xs(bench, level='index_code').ix[start:end]  # loc[(start:end,bench),slice(None)]
        self.benchmark.columns = ['adj_close']
        self.benchmark['return'] = self.benchmark['adj_close'].pct_change()

        date_len = len(self.benchmark.loc[start:end])

        # Retrieve assets from local data source (~/data/xccdata)
        data = pd.read_pickle(data_source + '/stock_prices_datetime') \
            [['code', 'date', 'adj_close', 'whether_trade']].set_index(['date', 'code'])
        adj_close = data['adj_close'].unstack()
        whether_trade = data['whether_trade'].unstack()

        # select the stock according to sample_num and select_type
        temp = [symbol for symbol in adj_close.columns \
                if whether_trade.loc[start:end, symbol].sum() / date_len > 0.95]
        if sample_num > len(temp):
            sample_num = len(temp)
        if select_type.startswith('a') or select_type.startswith('h'):
            self.stock_pools = temp[:sample_num]
        elif select_type.startswith('d') or select_type.startswith('t'):
            self.stock_pools = temp[-sample_num:]
        else:
            self.stock_pools = np.random.choice(temp, sample_num, replace=False)

        self.asset['adj_close'] = adj_close.loc[start:end, self.stock_pools]


        # Get returns, beta, alpha, and sharp ratio.
        ## get returns
        self.asset['return'] = self.asset['adj_close'].pct_change()
        ## get beta
        self.asset['beta'] = {symbol: self.asset['return'][symbol].cov(self.benchmark['return']) \
                                      / self.benchmark['return'].var() for symbol in self.asset['return'].columns}

        ## get alpha
        temp = np.empty_like(self.asset['return'])
        for i, symbol in enumerate(self.asset['return'].columns):
            temp[:, i] = self.asset['return'][symbol] - self.asset['beta'][symbol] * self.benchmark['return']
        self.asset['alpha'] = DataFrame(temp, index=self.asset['return'].index, columns=self.asset['return'].columns)

        ## get sharpe ratio
        self.asset['sharpe'] = {
        symbol: self.asset['return'][symbol].mean() / self.asset['return'][symbol].std() * sqrt(252) \
        for symbol in self.asset['return'].columns}



    def nplot(self, symbol, color='b', nval=0):
        if symbol == 'bench':
            tmp = self.benchmark['adj_close']
        else:
            tmp = self.asset['adj_close'][symbol]
        tmp = tmp / tmp[nval]
        tmp.plot(color=color, label=symbol)
        legend(loc='best', shadow=True, fancybox=True)

    def betas(self):
        betas = []
        for symbol in self.asset['adj_close'].columns:
            betas.append(self.asset['beta'][symbol])
        return Series(betas, index=self.asset['adj_close'].columns)

    def returns(self):
        returns = []
        for symbol in self.asset['adj_close'].columns:
            returns.append(self.asset['return'][symbol].dropna())
        return Series(returns, index=self.asset['adj_close'].columns)

    def cov(self):
        tmp = self.returns()
        tmpl = []
        for symbol in self.asset['adj_close'].columns:
            #tmpl.append(tmp[symbol])
            tmpl.append(tmp[symbol].values)
        return DataFrame(
            cov(array(tmpl)), index=self.asset['adj_close'].columns, columns=self.asset['adj_close'].columns)


    def get_w(self, kind='sharpe'):
        V = self.cov()
        iV = matrix(inv(V))

        if kind == 'characteristic':
            e = matrix(ones(len(self.asset['adj_close'].columns))).T
        elif kind == 'sharpe':
            suml = []
            for symbol in self.asset['adj_close'].columns:
                suml.append(self.returns()[symbol].sum())
            e = matrix(suml).T
        else:
            print('\n  *Error: There is no weighting for kind ' + kind)
            return

        num = iV * e
        denom = e.T * iV * e
        w = array(num / denom).flatten()
        return Series(w, index=self.asset['adj_close'].columns)

    def efficient_frontier_w(self, fp):
        wc = self.get_w(kind='characteristic')
        wq = self.get_w(kind='sharpe')

        fc = self.ret_for_w(wc).sum()
        fq = self.ret_for_w(wq).sum()

        denom = fq - fc
        w = (fq - fp) * wc + (fp - fc) * wq
        return Series(w / denom, index=self.asset['adj_close'].columns)

    def efficient_frontier(self, xi=0.01, xf=4, npts=100, scale=10):
        frontier = linspace(xi, xf, npts)

        i = 0
        rets = zeros(len(frontier))
        risk = zeros(len(frontier))
        sharpe = zeros(len(frontier))
        for f in frontier:
            w = self.efficient_frontier_w(f)
            tmp = self.ret_for_w(w)
            rets[i] = tmp.sum() * scale
            sharpe[i] = mean(tmp) / std(tmp) * sqrt(len(tmp))
            risk[i] = rets[i] / sharpe[i]
            i += 1
        return Series(rets, index=risk), sharpe.max()

#TODO 感觉画图太慢了，可能是前面的有效前沿计算效率较低，需要改进
    def efficient_frontier_plot(self, xi=0.01, xf=4, npts=100, scale=0.1,
                                col1='b', col2='r', newfig=1, plabel=''):
        eff, m = self.efficient_frontier()

        if newfig == 1:
            figure()

        plot(array(eff.index), array(eff), col1, linewidth=2,
             label="Efficient Frontier " + plabel)
        tmp = zeros(1)
        plot(hstack((tmp, array(eff.index))),
             hstack((0, m * array(eff.index))),
             col2, linewidth=2, label="Max Sharpe Ratio: %6.2g" % m)
        legend(loc='best', shadow=True, fancybox=True)
        xlabel('Risk %', fontsize=16)
        ylabel('Return %', fontsize=16)
        show()

    def min_var_w_ret(self, ret):
        V = self.cov()
        suml = []
        for symbol in self.asset['adj_close'].columns:
            suml.append(self.returns()[symbol].sum())
        e = matrix(suml).T
        iV = matrix(inv(V))
        num = iV * e * ret
        denom = e.T * iV * e
        return Series(array(num / denom).flatten(), index=self.asset['adj_close'].columns)

    def ret_for_w(self, w):
        tmp = self.returns()
        i = 0
        tmpl = []
        for symbol in tmp.keys():
            tmpl.append(tmp[symbol] * w[i])
            i += 1
        return Series(tmpl, index=tmp.keys()).sum()
class MeanVarAnalysis:
    def __init__(self,historical_data,weights=None,is_pandas=True):
        if weights is not None and historical_data.shape[1] != len(weights):
            raise ValueError('权重数目与股票数目不符')
        self.init_data=historical_data
        self.w=weights
        self.pandas=is_pandas
    def Mean(self):
        if self.pandas:
            return np.mean(self.init_data,axis=0).values
        else:
            return np.mean(self.init_data,axis=0)
    def Cov(self):
        return np.cov(self.init_data.T)
    # TODO 暂时添加，便于debug
    def preFrontierPortfolios(self):
        one = np.ones(self.init_data.shape[1])
        V_inv = pinv(self.Cov())
        e = self.Mean()
        A = one @ V_inv @ e
        B = e @ V_inv @ e
        C = one @ V_inv @ one
        D = B * C - A * A
        g = (B * one - A * e) @ V_inv / D
        h = (C * e - A * one) @ V_inv / D
        return g,h
    def FrontierPortfolios(self,Ep,onlyEfficient=False):
        '''

        :param Ep: Number,iterable
        :return: Array
            the p_weights of frontier portfolio
        '''
        one = np.ones(self.init_data.shape[1])
        V_inv = pinv(self.Cov())
        e = self.Mean()
        A = one @ V_inv @ e
        B = e @ V_inv @ e
        C = one @ V_inv @ one
        D = B * C - A * A
        g = (B * one - A * e) @ V_inv / D
        h = (C * e - A * one) @ V_inv / D
        if np.size(Ep) > 1:
            g = np.reshape(g, (len(g), 1))
            h = np.reshape(h, (len(h), 1))
            if onlyEfficient:
                Ep=Ep[Ep>A/C]
            Ep = np.reshape(Ep, (1, len(Ep)))
            return g + h @ Ep
        else:
            return g + h * Ep
    def MiniVarPortfolio(self):
        one = np.ones(self.init_data.shape[1])
        V_inv = np.linalg.pinv(self.Cov())
        C = one @ V_inv @ one
        return one@V_inv/C
    def PlotFrontier(self,Ep_range=None,sigma2=False,onlyEfficient=False):
        if Ep_range is None:
            Ep_range=np.arange(0.1,30.1,0.1)/100
        port_data=self.init_data@self.FrontierPortfolios(Ep_range,onlyEfficient=onlyEfficient)
        m=np.mean(port_data,axis=0)
        if sigma2:
            sigma = np.var(port_data, axis=0)
        else:
            sigma=np.std(port_data,axis=0)
        plt.plot(sigma,m)
        plt.show()
        #return m,sigma







