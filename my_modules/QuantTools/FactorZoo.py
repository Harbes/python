# todo 几点建议：.loc换成.reindex;将SustainableGrowth函数作为模板推广至其他函数；财报数据出现加减时，引入fillna
import pandas as pd
import numpy as np
from pandas.tseries.offsets import Day,DateOffset,MonthEnd#,YearEnd,Week

def AssetsToMarket(tot_assets,market_cap,date_list,annually=True,pub_date=None,most_recent=False):
    ## todo 需要仔细检查
    '''
    Assets-to-market, which is total assets for the fiscal year divided by fiscal-year-end market capitalization.
    for example:
        tot_assets=BS['AT'].unstack()
market_cap=PV['Dmktcap'].unstack()
date_list=pd.date_range('2015-01-01','2017-12-31', freq='Q');
tmp1=AssetsToMarket(tot_assets,market_cap,date_list,annually=True,pub_date=pub_date,most_recent=True)
    :param tot_assets:
    :param market_cap:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp=tot_assets.resample('Q').fillna(method=None)
    if annually:
        tmp=tmp[tmp.index.month==12]
    if pub_date is None:
        if most_recent:
            tmp.index = tmp.index + MonthEnd(6) + Day()
            return (tmp.resample('D').ffill()/market_cap.resample('D').ffill().shift(1)).reindex(date_list)
        else:
            tmp=tmp/market_cap.resample('D').ffill().reindex(tmp.index)
            tmp.index = tmp.index + MonthEnd(6) + Day()
            return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        if most_recent:
            return (tmp2.resample('D').ffill().ffill(limit=365).shift(1) / market_cap.resample('D').ffill().shift(1)).reindex(
                date_list)
        else:
            tmp2= tmp2 / market_cap.resample('D').ffill().reindex(tmp2.index)
            return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def BookToMarket(book_value,market_cap,date_list,annually=True,pub_date=None):
    '''
    Book-to-market equity, which is the book value of equity for fiscal year divided by fiscal-year-end market capitalization
    For example:
        BS=pd.read_pickle(DPath+'BS')
        PV=pd.read_pickle(DPath+'PV')
        book_value=BS['EQU'].unstack()
market_cap=PV['Dmktcap'].unstack()
tmp=BookToMarket(book_value,market_cap,date_list,annually=True,pub_date=None)
    :param book_value:
    :param market_cap:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp = book_value.resample('Q').fillna(method=None)
    if annually:
        tmp = tmp[tmp.index.month == 12]
    if pub_date is None:
        tmp = tmp / market_cap.resample('D').ffill().reindex(tmp.index)
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        tmp2 = tmp2 / market_cap.resample('D').ffill().reindex(tmp2.index)
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def OperatingCashFlowToPrice(oper_cash,market_cap,date_list,pub_date=None):
    '''
    Cash flow-to-price, which is operating cash flows divided by fiscal-year-end market capitalization.
    For example:
        PV=pd.read_pickle(DPath+'PV')
        CF=pd.read_pickle(DPath+'CF')
        oper_cash=CF['NCPOA'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=oper_cash.loc['2004':].index[1:20]
        tmp=OperatingCashFlowToPrice(oper_cash,market_cap,date_list,pub_date=None)
    :param oper_cash:
    :param market_cap:
    :param date_list:
    :param pub_date:
    :return:
    '''
    tmp = oper_cash[oper_cash.index.month == 12]
    if pub_date is None:
        tmp = tmp / market_cap.resample('D').ffill().reindex(tmp.index)
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        tmp2 = tmp2 / market_cap.resample('D').ffill().reindex(tmp2.index)
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def DebtToEquity(tot_lia,market_cap,date_list,annually=True,pub_date=None):
    '''
    Debt-to-equity ratio, which is total liabilities divided by fiscal-year-end market capitalization.
    For example:
        BS=pd.read_pickle(DPath+'BS')
        PV=pd.read_pickle(DPath+'PV')
        tot_lia=BS['LB'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=total_lia.loc['2004':].index[1:20]
        tmp=DebtToEquity(tot_lia,market_cap,date_list,annually=True,pub_date=None)
    :param tot_lia:
    :param market_cap:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp = tot_lia.resample('Q').fillna(method=None)
    if annually:
        tmp = tmp[tmp.index.month == 12]
    if pub_date is None:
        tmp = tmp / market_cap.resample('D').ffill().reindex(tmp.index)
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        tmp2 = tmp2 / market_cap.resample('D').ffill().reindex(tmp2.index)
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def LongTermDebtToMarketEquity(long_debt,market_cap,date_list,annually=True,pub_date=None):
    '''
    Long term debt-to-market equity, which is long term liabilities divided by fiscal-year-end market capitalization.
    For example:
        long_debt=BS['NCL'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=long_debt.loc['2004':].index[1:20]
        tmp=LongTermDebtToMarketEquity(long_debt,market_cap,date_list,annually=True,pub_date=None)
    :param long_debt:
    :param market_cap:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp = long_debt.resample('Q').fillna(method=None)
    if annually:
        tmp = tmp[tmp.index.month == 12]
    if pub_date is None:
        tmp = tmp / market_cap.resample('D').ffill().reindex(tmp.index)
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        tmp2 = tmp2 / market_cap.resample('D').ffill().reindex(tmp2.index)
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def DividendToPrice(dividend,market_cap,date_list,pub_date=None):
    '''
    Dividend-to-price ratio, which is annual total dividends payouts divided by fiscal-year-end market capitalization.
    For example:
        dividend=BS['DVP'].unstack()
market_cap=PV['Dmktcap'].unstack()
tmp=DividendToPrice(dividend,market_cap,date_list,pub_date=None)
    :param dividend:
    :param market_cap:
    :param date_list:
    :return:
    '''
    tmp = dividend[dividend.index.month == 12]*1e6
    if pub_date is None:
        tmp = tmp / market_cap.resample('D').ffill().reindex(tmp.index)
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        tmp2 = tmp2 / market_cap.resample('D').ffill().reindex(tmp2.index)
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def EarningsToPrice(annual_inc_bef_extra,market_cap,date_list,pub_date=None):
    '''
    Earnings-to-price, which is annual income before extraordinary items divided by fiscal-year-end market capitalization.
    For example:
        annual_inc_bef_extra=IS['IBT'].unstack()-IS['TAX'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=annual_inc_bef_extra.loc['2004':].index[1:20]
        tmp=EarningsToPrice(annual_inc_bef_extra,market_cap,date_list)
    :param annual_income_before:
    :param market_cap:
    :param date_list:
    :return:
    '''
    tmp = annual_inc_bef_extra[annual_inc_bef_extra.index.month == 12]
    if pub_date is None:
        tmp = tmp / market_cap.resample('D').ffill().reindex(tmp.index)
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        tmp2 = tmp2 / market_cap.resample('D').ffill().reindex(tmp2.index)
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def LiabilityGrowth(tot_lia,date_list,annually=True,pub_date=None):
    '''
    Liability growth, which is the annual change in total liabilities divided by 1-year-lagged total liabilities.
    For example:
        tot_lia=BS['LB'].unstack()
        tmp2=LiabilityGrowth(total_lia,total_lia.index[15:],annually=False)
    :param total_lia:
    :param date_list:
    :return:
    '''
    if annually:
        tmp = tot_lia[tot_lia.index.month == 12].pct_change(fill_method=None)
    else:
        tmp=tot_lia.resample('Q').fillna(method=None).pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def AdjOperatingCashFlowToPrice(oper_cash,market_cap,date_list,pub_date=None):
    ## todo 与前述的OperatingCashFlowToPrice重复
    '''
    Operating cash flow-to-price, which is operating cash flow(measured as earnings adjusted for depreciation and
        working capital accruals) divided by fiscal-year-end market capitalization.
    For example:
        oper_cash=CF['NCPOA'].unstack()
tmp=OperatingCashFLowToPrice(oper_cash,market_cap,date_list,pub_date=pub_date)
    :param oper_cash:
    :param market_cap:
    :param date_list:
    :param pub_date:
    :return:
    '''
    tmp = oper_cash[oper_cash.index.month == 12]
    if pub_date is None:
        tmp = tmp / market_cap.resample('D').ffill().reindex(tmp.index)
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        tmp2 = tmp2 / market_cap.resample('D').ffill().reindex(tmp2.index)
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def PayoutYield(annual_inc_bef_extra,book_value,market_cap,date_list,pub_date=None):
    '''
    Payout yield, which is annual income before extraordinary items minus the change of book equity divided by
        fiscal year end market capitalization.
    For example:
        annual_inc_bef_extra=IS['IBT'].unstack()-IS['TAX'].unstack()
        book_value=BS['EQU'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=book_value.index[15:]
        tmp=PayoutYield(annual_inc_bef_extra,book_value,market_cap,date_list)
    :param annual_inc_minus_chang_book:
    :param market_cap:
    :param date_list:
    :return:
    '''
    tmp = annual_inc_bef_extra[annual_inc_bef_extra.index.month==12]-book_value[book_value.index.month==12].diff()
    if pub_date is None:
        tmp = tmp / market_cap.resample('D').ffill().reindex(tmp.index)
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        tmp2 = tmp2 / market_cap.resample('D').ffill().reindex(tmp2.index)
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def Reversal_60_13(adj_prc,date_list,fast=True):
    ## todo 本质上与momentum相同
    '''
    Reversal, which is cumulative returns from months t-60 to t-13.
    For example:
        adj_price=PV['Adclsprc'].unstack()
        tmp=Reversal_60_13(adj_price,date_list)
    :param adj_price:
    :param date_list:
    :return:
    '''

    p0 = adj_prc.resample('D').ffill().ffill(limit=7)
    p1 = adj_prc.resample('D').ffill().ffill(limit=7)
    p0.index = p0.index + Day(365*5+1)
    p1.index = p1.index + Day(365)
    return (p1.reindex(date_list) - p0.reindex(date_list)) / p0.reindex(date_list)
def SustainableGrowth(book_value,date_list,annually=True,pub_date=None):
    '''
    Sustainable growth, which is annual growth in book value of equity
    For example:
        book_value=BS['EQU'].unstack()
        pub_date=pd.read_pickle(DPath+'PubDate')['ActlDt']
        tmp=SustainableGrowth(book_value,date_list,annually=False,pub_date=pub_date)
    :param book_value:
    :param date_list:
    :return:
    '''
    if annually:
        tmp=book_value.loc[book_value.index.month==12].pct_change(fill_method=None)
    else:
        tmp=book_value.resample('Q').fillna(method=None).pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['growth']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['growth'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def SalesGrowthMinusInventoryGrowth(sales,inventory,date_list,pub_date=None):
    '''
    Sales growth minus inventory growth, which is annual growth in sales minus annual growth in inventory.
    For example:
        sales=IS['REV'].unstack()
        inventory=BS['INVY'].unstack()
        tmp=SalesGrowthMinusInventoryGrowth(sales,inventory,date_list,pub_date=pub_date)
    :param sales:
    :param inventory:
    :param date_list:
    :return:
    '''
    tmp = sales[sales.index.month==12].pct_change(fill_method=None)-inventory[inventory.index.month==12].pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['growth'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['growth'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def SalesToPrice(sales,market_cap,date_list,pub_date=None):
    '''
    Sales-to-price, which is the annual operating revenue divided by fiscal-year-end market capitalization.
    For example:
        sales=IS['REV'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        tmp=SalesToPrice(sales,market_cap,date_list,pub_date=None)
    :param sales:
    :param market_cap:
    :param date_list:
    :return:
    '''
    tmp = sales[sales.index.month==12]
    if pub_date is None:
        tmp2=tmp / market_cap.resample('D').ffill().reindex(tmp.index)
        tmp2.index = tmp2.index + MonthEnd(6) + Day()
        return tmp2.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['sales'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['sales'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)/market_cap.resample('D').ffill().shift(1).reindex(date_list)
def TaxGrowth(tax,date_list,pub_date=None):
    '''
    Tax growth, which is annual change in taxes payable divided by 1-year-lagged taxes payable.
    For example:
        tax=IS['TAX'].unstack()
        tmp1=TaxGrowth(tax,date_list,pub_date=pub_date)
    :param taxes_payable:
    :param date_list:
    :return:
    '''
    tmp = tax[tax.index.month == 12].pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['growth'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['growth'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def Accruals(accruals,oper_cash,tot_assets,date_list,pub_date=None):
    ## todo accruals的计算与EarningsToPrice、PayoutYield似乎一样？是不是错了？
    '''
    Accruals, which is annual income before extraordinary items minus operating cash flows divided by average total assets.
    For example:
        accruals=IS['IBT'].unstack()-IS['TAX'].unstack()
        oper_cash=CF['NCPOA'].unstack()
        tot_assets=BS['AT'].unstack()
        tmp1=Accruals(accruals,oper_cash,tot_assets,date_list,pub_date=pub_date)
    :param annual_inc_minus_oper_cashflow:
    :param total_assets:
    :param date_list:
    :return:
    '''
    tmp = (accruals[accruals.index.month == 12]-oper_cash[oper_cash.index.month==12])/tot_assets[tot_assets.index.month==12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['accruals'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['accruals'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def PercentAccruals(tot_profit,oper_cash,net_profit,date_list,pub_date=None):
    '''
    Percent accruals, which is total profit minus operating cash flow divided by net profit.
    For example:
        tot_profit=IS['IBT'].unstack()
        oper_cash=CF['NCPOA'].unstack()
        net_profit=IS['NI'].unstack()
        tmp1=PercentAccruals(tot_profit,oper_cash,net_profit,date_list,pub_date=pub_date)
    :param total_profit_minus_oper_cashflow:
    :param net_profit:
    :param date_list:
    :return:
    '''
    tmp = (tot_profit[tot_profit.index.month == 12] - oper_cash[oper_cash.index.month == 12]) / net_profit[
        net_profit.index.month == 12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['per'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['per'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def CapitalExpenditureGrowth(PPE,depre,date_list,pub_date=None):
    '''
    Capital expenditure growth, which is the annual change in capital expenditure divided by 1-
        year-lagged capital expenditure.
    capital expenditure=delta(PPE)+depre
    For example:
        PPE=BS['PPE'].unstack()
        depre=CF['DFAOGA'].unstack()
        tmp1=CapitalExpenditureGrowth(PPE,depre,date_list,pub_date=pub_date)
    :param capital_expenditure:
    :param date_list:
    :return:
    '''
    tmp=(PPE[PPE.index.month==12].diff()+depre[depre.index.month==12]).pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['growth'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['growth'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def ChangeInShareholdersEquity(book_value,tot_assets,date_list,annually=True,pub_date=None):
    '''
    Change in shareholders’ equity, which is the annual change in book equity divided by 1-year-lagged total assets.
    For example:
        book_value=BS['EQU'].unstack()
        tot_assets=BS['AT'].unstack()
        tmp2=ChangeInShareholdersEquity(book_value,tot_assets,date_list,annually=False,pub_date=pub_date)
    :param book_value:
    :param total_asset:
    :param date_list:
    :return:
    '''
    if annually:
        tmp=book_value.loc[book_value.index.month==12].diff()/tot_assets[tot_assets.index.month==12].shift(1)
    else:
        tmp=book_value.diff()/tot_assets.shift(1)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(date_list)
def ChangesInPPEandInventoryToAssets(PPE,inventory,tot_assets,date_list,annually=True,pub_date=None):
    ## todo 对于没有inventory(例如银行)的数据，舍弃(nan)还是取值为0？？？
    '''
    Changes in PPE and inventory-to-assets, which is the annual change in gross property, plant, and equipment
        plus the annual change in inventory scaled by 1-year-lagged total assets.
    For example:
        PPE=BS['PPE'].unstack()
        inventory=BS['INVY'].unstack()
        tot_assets=BS['AT'].unstack()
        tmp2=ChangesInPPEandInventoryToAssets(PPE,inventory.fillna(0),tot_assets,date_list,annually=False,pub_date=pub_date)
    :param PPE:
    :param inventory:
    :param tot_assets:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    if annually:
        tmp=(PPE+inventory).loc[PPE.index.month==12].diff()/tot_assets[tot_assets.index.month==12].shift(1)
    else:
        tmp=(PPE+inventory).diff()/tot_assets.shift(1)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list,columns=tmp.columns)
def InvestmentToAssets(tot_assets,date_list,annually=True,pub_date=None):
    '''
    Investment-to-assets, which is the annual change in total assets divided by 1-year-lagged total assets.
    For example:
        tot_assets=BS['AT'].unstack()
        tmp2=InvestmentToAssets(tot_assets,date_list,annually=False,pub_date=pub_date)
    :param tot_assets:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    if annually:
        tmp=tot_assets[tot_assets.index.month==12].pct_change(fill_method=None)
    else:
        tmp=tot_assets.pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list,columns=tmp.columns)
def InventoryChange(inventory,tot_assets,date_list,annually=True,pub_date=None):
    '''
    Inventory change, which is the annual change in inventory scaled by two-year average of total assets.
    For example:
        inventory=BS['INVY'].unstack()
        tot_assets=BS['AT'].unstack()
        tmp2=InventoryChange(inventory,tot_assets,date_list,annually=False,pub_date=pub_date)
    :param inventory:
    :param tot_assets:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    if annually:
        tmp=2.0*inventory[inventory.index.month==12].diff()/\
            (tot_assets[tot_assets.index.month==12]+tot_assets[tot_assets.index.month==12].shift(1))
    else:
        tmp=2.0*inventory.diff()/(tot_assets+tot_assets.shift(1))
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list,columns=tmp.columns)
def InventoryGrowth(inventory,date_list,annually=True,pub_date=None):
    '''
    Inventory growth, which is the annual change in inventory divided by 1-year-lagged inventory.
    For example:
        inventory=BS['INVY'].unstack()
        tmp2=InventoryGrowth(inventory,date_list,annually=False,pub_date=pub_date)
    :param inventory:
    :param date_List:
    :param annually:
    :param pub_date:
    :return:
    '''
    if annually:
        tmp=inventory[inventory.index.month==12].pct_change(fill_method=None)
    else:
        tmp=inventory.pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list,columns=tmp.columns)
def NetOperatingAssets(fin_assets,fin_lia,book_value,tot_assets,date_list,annually=True,pub_date=None):
    '''
    Net operating assets, which is operating assets minus operating liabilities scaled by total assets
    净经营资产=经营性资产-经营性负债=净金融负债+所有者权益
    For example:
        fin_assets=(BS['HLTFA'].fillna(0)+BS['BBSFA'].fillna(0)+BS['AFSFA'].fillna(0)).unstack()
        fin_lia=(BS['STBL'].fillna(0)+BS['HLTFL'].fillna(0)+BS['INTP'].fillna(0)+BS['DVP'].fillna(0)+BS['LTBL'].fillna(0)).unstack()
        book_value=BS['EQU'].unstack()
        tot_assets=BS['AT'].unstack()
        tmp2= NetOperatingAssets(fin_assets,fin_lia,book_value,tot_assets,date_list,annually=False,pub_date=pub_date)
    :param oper_assets:
    :param oper_lia:
    :param tot_assets:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp = (fin_lia - fin_assets + book_value) / tot_assets
    if annually:
        tmp=tmp[tmp.index.month==12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list,columns=tmp.columns)
def AssetTurnover(fin_assets,fin_lia,book_value,sales,date_list,pub_date=None):
    '''
    Asset turnover, which is sales divided by 1-year-lagged net operating assets
    For example:
        fin_assets=(BS['HLTFA'].fillna(0)+BS['BBSFA'].fillna(0)+BS['AFSFA'].fillna(0)).unstack()
        fin_lia=(BS['STBL'].fillna(0)+BS['HLTFL'].fillna(0)+BS['INTP'].fillna(0)+BS['DVP'].fillna(0)+BS['LTBL'].fillna(0)).unstack()
        book_value=BS['EQU'].unstack()
        sales=IS['REV'].unstack()
        tmp1=AssetTurnover(fin_assets,fin_lia,book_value,sales,date_list,pub_date=pub_date)
    :param fin_assets:
    :param fin_lia:
    :param book_value:
    :param sales:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp =sales/(fin_lia - fin_assets + book_value)
    tmp=tmp[tmp.index.month==12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['change'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['change'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def CashFlowOverAssets(oper_cash,tot_assets,date_list,pub_date=None):
    '''
    Cash flow over assets, which is cash flow from operation scaled by total assets.
    For example:
        oper_cash=CF['NCPOA'].unstack()
        tot_assets=BS['AT']
        tmp2=CashFlowOverAssets(oper_cash,tot_assets,date_list,annually=False,pub_date=pub_date)
    :param oper_cash:
    :param tot_assets:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp = oper_cash[oper_cash.index.month==12]/tot_assets[tot_assets.index.month==12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def CashProductivity(trad_share,long_debt,tot_assets,cash_eq,date_list,annually=True,pub_date=None):
    ## todo 不同频率数据加减时，使用pub_date是否恰当？？？
    '''
    Cash productivity, which is market value of tradable shares plus long-term liabilities minus total
        assets scaled by cash and cash equivalents.
    For example:
        long_debt=BS['NCL'].unstack()
        tot_assets=BS['AT'].unstack()
        cash_eq=CF['CEEPBL'].unstack()
        trad_share=PV['Domktcap'].unstack()
        tmp1=CashProductivity(trad_share,long_debt,tot_assets,cash_eq,date_list,annually=True,pub_date=pub_date)
    :param trad_share:
    :param long_debt:
    :param tot_assets:
    :param cash_eq:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp = (trad_share.resample('D').ffill().ffill(limit=365).loc[long_debt.index]+long_debt-tot_assets)/cash_eq
    if annually:
        tmp = tmp[tmp.index.month == 12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def CashToAssets(cash_eq,tot_assets,date_list,annually=True,pub_date=None):
    '''
    Cash-to-assets, which is cash and cash equivalents divided by the two-year average of total assets.
    For example:
        tot_assets=BS['AT'].unstack()
        cash_eq=CF['CEEPBL'].unstack()
        tmp2=CashToAssets(cash_eq,tot_assets,date_list,annually=False,pub_date=pub_date)
    :param cash_eq:
    :param tot_assets:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    if annually:
        tmp=2.0*cash_eq[cash_eq.index.month==12]/(tot_assets[tot_assets.index.month==12]+tot_assets[tot_assets.index.month==12].shift(1))
    else:
        tmp=2.0*cash_eq/(tot_assets+tot_assets.shift(1))
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def CapitalTurnover(sales,tot_assets,date_list,pub_date=None):
    '''
    Capital turnover, which is sales divided by 1-year-lagged total assets.
    For example:
        sales=IS['REV'].unstack()
        tot_assets=BS['AT'].unstack()
        tmp1=CapitalTurnover(sales,tot_assets,date_list,pub_date=pub_date)
    :param sales:
    :param tot_assets:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp=sales[sales.index.month==12]/tot_assets[tot_assets.index.month==12].shift(1)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def ReturnOnCapital(inc_bef_tax,fin_exp,net_working,net_fixed,date_list,pub_date=None):
    ## todo inc_bef_tax 与 前述的tot_profit 是否相等？？？
    '''
    Return on Capital = EBIT/(Net Working Capital + Net Fixed Assets)
    For example:
        inc_bef_tax=IS['IBT'].unstack()
        fin_exp=IS['FINEXP'].unstack()
        net_working=(BS['CA']-BS['CL']).unstack()
        net_fixed=(BS['PPE']-CF['DFAOGA']).unstack()
        tmp1=ReturnOnCapital(inc_bef_tax,fin_exp,net_working,net_fixed,pub_date=pub_date)
    :param inc_bef_tax:
    :param fin_exp:
    :param net_working:
    :param net_fixed:
    :param pub_date:
    :return:
    '''
    tmp = (inc_bef_tax+fin_exp)/(net_fixed+net_working)
    tmp=tmp[tmp.index.month==12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().shift(1).reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def EarningsYield(inc_bef_tax,fin_exp,market_cap,tot_lia,cash_eq,date_list):
    '''
    Earnings yield, which is earnings before interests and taxes divided by enterprise value.
    (Enterprise value formula:EV = Market capitalization + Total debt − Cash and cash equivalents)
    For example:
        inc_bef_tax=IS['IBT'].unstack()
        fin_exp=IS['FINEXP'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        tot_lia=BS['LB'].unstack()
        cash_eq=CF['CEEPBL'].unstack()
        tmp=EarningsYield(inc_bef_tax,fin_exp,market_cap,tot_lia,cash_eq,date_list)
    :param inc_bef_tax:
    :param fin_exp:
    :param market_cap:
    :param tot_lia:
    :param cash_eq:
    :param pub_date:
    :return:
    '''
    tmp = (inc_bef_tax + fin_exp) / (market_cap.resample('D').ffill().ffill(limit=365).loc[tot_lia.index]+tot_lia-cash_eq)
    tmp = tmp[tmp.index.month == 12]
    tmp.index = tmp.index + MonthEnd(6)
    return tmp.resample('D').ffill().shift(1).reindex(date_list)
def GrossMargins(oper_rev,oper_exp,date_list,pub_date=None):
    '''
    Gross margins, which is operating revenue minus operating expenses divided by 1-year-lagged operating revenue.
    For example:
        oper_rev=IS['OPREV'].unstack()
        oper_exp=IS['COGS'].fillna(0).unstack()
        tmp1=GrossMargins(oper_rev,oper_exp,date_list,pub_date=pub_date)
    :param oper_rev:
    :param oper_exp:
    :param date_list:
    :param pub_date:
    :return:
    '''
    tmp = (oper_rev[oper_rev.index.month==12]-oper_exp[oper_exp.index.month==12])/oper_rev[oper_rev.index.month==12].shift(1)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def GrossProfitability(oper_rev,oper_exp,tot_assets,date_list,pub_date=None):
    '''
    Gross profitability ratio, which is the quarterly operating revenue minus quarterly operating expenses
        divided by the average of current quarterly total assets and 1-quarter-lagged total assets.
    For example:
        oper_rev=IS['OPREV'].unstack()
        oper_exp=IS['COGS'].fillna(0).unstack()
        tot_assets=BS['AT'].unstack();tot_assets[tot_assets<=0.0]=np.nan
        tmp1=GrossProfitability(oper_rev,oper_exp,tot_assets,date_list,pub_date=None)
    :param oper_rev_q:
    :param oper_exp_q:
    :param tot_assets_q:
    :param date_list:
    :param pub_date:
    :return:
    '''
    oper_rev_q=oper_rev.resample('Q').fillna(method=None).diff()
    #oper_rev_q=oper_rev_q-oper_rev_q.shift(1)
    oper_rev_q[oper_rev_q.index.month==3]=oper_rev[oper_rev.index.month==3]
    oper_exp_q = oper_exp.resample('Q').fillna(method=None).diff()
    #oper_exp_q=oper_exp-oper_exp.shift(1)
    oper_exp_q[oper_exp_q.index.month == 3] = oper_exp[oper_exp.index.month == 3]
    tmp = 2.0*(oper_rev_q-oper_exp_q)/(tot_assets+tot_assets.shift(1))
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def NetPayoutOverProfit(net_inc,book_value,tot_profit,date_list,pub_date=None):
    '''
    Net payout over profits, which is the sum of total net payout (net income minus changes in book equity)
        divided by total profits.
    For example:
        net_inc=IS['NI'].unstack()
        book_value=BS['EQU'].unstack()
        tot_profit=IS['IBT'].unstack()
        date_list=book_value.loc['2005':].index
        tmp1=NetPayoutOverProfit(net_inc,book_value,tot_profit,date_list,pub_date=pub_date)
    :param net_inc:
    :param book_value:
    :param tot_profit:
    :param date_list:
    :param pub_date:
    :return:
    '''
    tmp = (net_inc[net_inc.index.month==12]-book_value[book_value.index.month==12].diff())/\
          tot_profit[tot_profit.index.month==12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def ReturnOnOperatingAsset(oper_inc,depre,fin_assets,fin_lia,date_list,pub_date=None):
    '''
    Return on net operating assets, which is operating income after depreciation divided by 1-year lagged net operating assets.
    For example:
        oper_inc=IS['OI'].unstack()
depre=CF['DFAOGA'].unstack()
fin_assets = (BS['HLTFA'].fillna(0) + BS['BBSFA'].fillna(0) + BS['AFSFA'].fillna(0)).unstack()
fin_lia = (BS['STBL'].fillna(0) + BS['HLTFL'].fillna(0) + BS['INTP'].fillna(0) + BS['DVP'].fillna(0) + BS['LTBL'].fillna(0)).unstack()
tmp1=ReturnOnOperatingAsset(oper_inc,depre,fin_assets,fin_lia,date_list,pub_date=pub_date)
    :param oper_inc:
    :param depre:
    :param oper_asset:
    :param oper_lia:
    :param date_list:
    :param pub_date:
    :return:
    '''
    net_oper_asset=fin_assets-fin_lia
    tmp=(oper_inc[oper_inc.index.month==12]-depre[depre.index.month==12])/net_oper_asset[net_oper_asset.index.month==12].shift(1)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def ReturnOnAssets(oper_inc,tot_assets,date_list,pub_date=None):
    '''
    Return on assets, which is quarterly total operating profit divided by the average of current quarterly total assets
        and 1-quarter-lagged total assets.
    For example:
        oper_inc=IS['OI'].unstack()
tot_assets=BS['AT'].unstack()
tmp1=ReturnOnAssets(oper_inc,tot_assets,date_list,pub_date=pub_date)
    :param oper_inc:
    :param tot_assets:
    :param date_list:
    :param pub_date:
    :return:
    '''
    oper_inc_q=oper_inc.resample('Q').fillna(method=None).diff()
    #oper_inc_q=oper_inc_q-oper_inc_q.shift(1)
    oper_inc_q[oper_inc_q.index.month==3]=oper_inc[oper_inc.index.month==3]
    tmp=2.0*oper_inc_q/(tot_assets+tot_assets.shift(1))
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def ReturnOnEquity(net_inc,book_value,date_list,pub_date=None):
    '''
    Return on equity, which is quarterly net income divided by the average of current quarterly total shareholders’ equity and 1-quarter-lagged shareholders’ equity.
    For example:
        net_inc=IS['NI'].unstack()
book_value=BS['EQU'].unsatck()
tmp1=ReturnOnEquity(net_inc,book_value,date_list,pub_date=pub_date)
    :param net_inc:
    :param book_value:
    :param date_list:
    :param pub_date:
    :return:
    '''
    net_inc_q=net_inc.resample('Q').fillna(method=None).diff()
    #net_inc_q=net_inc_q-net_inc_q.shift(1)
    net_inc_q[net_inc_q.index.month==3]=net_inc[net_inc.index.month==3]
    tmp=net_inc_q*2.0/(book_value+book_value.resample('Q').fillna(method=None).shift(1))
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def ReturnOnInvestedCapital(inc_bef_tax,fin_exp,non_oper_inc,market_cap,tot_lia,cash_eq,date_list):
    '''
    Return on invested capital, which is t annual earnings before interest and taxes minus non- operating income divided by non-cash enterprise value.
    For example:
        inc_bef_tax=IS['IBT'].unstack()
fin_exp=IS['FINEXP'].unstack().fillna(0)
non_oper_inc=IS['NOINC'].unstack().fillna(0)
market_cap=PV['Dmktcap'].unstack()
tot_debt=BS['LB'].unstack()
cash_eq=CF['CEEPBL'].unstack().fillna(0)
tmp=ReturnOnInvestedCapital(inc_bef_tax,fin_exp,non_oper_inc,market_cap,tot_debt,cash_eq,date_list)
    :param inc_bef_tax:
    :param fin_exp:
    :param market_cap:
    :param tot_debt:
    :param cash_eq:
    :param date_list:
    :return:
    '''
    tmp = (inc_bef_tax + fin_exp - non_oper_inc) / (market_cap.resample('D').ffill().ffill(limit=365).loc[tot_lia.index] + tot_lia - cash_eq)
    tmp = tmp[tmp.index.month == 12]
    tmp.index = tmp.index + MonthEnd(6) + Day()
    return tmp.resample('D').ffill().reindex(date_list)
def TexableIncomeToBookIncome(inc_bef_tax,net_inc,date_list,pub_date=None):
    '''
    Taxable income-to-book income, which is pretax income divided by net income.
    For example:
        inc_bef_tax=IS['IBT'].unstack()
net_inc=IS['NI'].unstack()
tmp1=TexableIncomeToBookIncome(inc_bef_tax,net_inc,date_list,pub_date=pub_date)
    :param inc_bef_tax:
    :param net_inc:
    :param date_list:
    :param pub_date:
    :return:
    '''
    tmp=inc_bef_tax/net_inc
    tmp=tmp[tmp.index.month==12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def ZScore(net_working,tot_assets,retained_earnings,ebit,market_cap,tot_lia,sales,date_list):
    '''
    Z-score, we follow Dichev (1998) to construct Z-score = 1.2 × (working capital / total assets) +
        1.4 × (retained earnings / total assets) + 3.3 × (EBIT / total assets) +
        0.6 × (market value of equity / book value of total liabilities) + (sales / total assets).
    For example:
        net_working=BS['CA'].unstack().fillna(0)-BS['CL'].unstack().fillna(0)
retained_earnings=BS['SR'].unstack().fillna(0)+BS['UDP'].unstack().fillna(0)
ebit=IS['IBT'].unstack()-IS['FINEXP'].unstack().fillna(0)
market_cap=PV['Dmktcap'].unstack()
tot_lia=BS['LB'].unstack()
sales=IS['REV'].unstack()
tmp=ZScore(net_working,tot_assets,retained_earnings,ebit,market_cap,tot_lia,sales,date_list)
    :param net_working:
    :param tot_assets:
    :param ebit:
    :param market_cap:
    :param tot_lia:
    :param sales:
    :param date_list:
    :return:
    '''
    tmp = 1.2*net_working/tot_assets+1.4*retained_earnings/tot_assets+3.3*ebit/tot_assets+\
          0.6*market_cap.resample('D').ffill().ffill(limit=365).loc[tot_lia.index]/tot_lia+sales/tot_assets
    tmp = tmp[tmp.index.month == 12]
    tmp.index = tmp.index + MonthEnd(6) +Day()
    return tmp.resample('D').ffill().reindex(date_list)
def ChangeIn6MonthMomentum(adj_prc,date_list):
    '''
    Change in 6-month momentum, which is cumulative returns from months t-6 to t-1 minus months t-12 to t-7.
    For example:
        adj_prc=PV['Adclsprc'].unstack()
tmp=ChangeIn6MonthMomentum(adj_prc,date_list)
    :param adj_prc:
    :param date_list:
    :return:
    初始Code【运行较慢】
    ch_mom6=pd.DataFrame(index=date_list,columns=adj_prc.columns)
    for i in date_list:
        ch_mom6.loc[i]=adj_prc.loc[i-DateOffset(months=6):i-Day()].iloc[[0,-1]].pct_change().iloc[1]-\
                       adj_prc.loc[i-DateOffset(months=12):i-DateOffset(months=6)].iloc[[0,-1]].pct_change().iloc[1]
    return ch_mom6
    '''
    p0 = adj_prc.resample('D').ffill().ffill(limit=7)
    p1 = adj_prc.resample('D').ffill().ffill(limit=7)
    p0.index = p1.index + Day(183)#pd.Series(p0.index) + DateOffset(months=6)
    p1.index = p1.index + Day(1)
    tmp1 = (p1.reindex(date_list) - p0.reindex(date_list)) / p0.reindex(date_list)
    p0.index = p0.index + Day(183)
    p1.index = p1.index + Day(183)
    return tmp1 - (p1.reindex(date_list) - p0.reindex(date_list)) / p0.reindex(date_list)
def IndustryMomentum(adj_prc,sector,sec_sign,date_list,period_start=Day(183),period_end=Day()):
    ## 将个股return替换为行业return
    '''
    sector=pd.read_pickle(DPath+'sector_datetime')
adj_prc=pd.read_pickle(DPath+'PVd')['Adclsprc'].unstack()
sec_sign=np.arange(29.0)
    :param adj_prc:
    :param sector:
    :param sec_sign:
    :param date_list:
    :param weight:
    :return:
    '''
    #sector.loc[date_list[-1]]=np.nan
    #sector_tmp=sector.resample('D').ffill().ffill()
    sector_tmp = sector.iloc[-1]
    #mom = pd.DataFrame(index=date_list, columns=adj_prc.columns&sector.columns)
    #for i in date_list:
    #    mom.loc[i] = adj_prc.loc[i - period_start:i - period_end].iloc[[0, -1]].pct_change().iloc[1]
    p0 = adj_prc.resample('D').ffill().ffill(limit=7)
    p1 = adj_prc.resample('D').ffill().ffill(limit=7)
    p0.index = p0.index + period_start
    p1.index = p1.index + period_end
    mom=((p1.reindex(date_list) - p0.reindex(date_list)) / p0.reindex(date_list))[adj_prc.columns&sector.columns]
    for i in mom.index:
        for j in sec_sign:
            #mom.loc[i][sector_tmp.loc[i] == j]=mom.loc[i][sector_tmp.loc[i]==j].mean()
            mom.loc[i][sector_tmp == j] = mom.loc[i][sector_tmp == j].mean()
    return mom
def Momentum(adj_prc,date_list,period_start=Day(30),period_end=Day()):
    '''
    1-month momentum, which is one-month cumulative returns.(short-term reversal)
    For example:
        adj_prc = PV['Adclsprc'].unstack()
tmp=Momentum1Month(adj_prc,date_list)
    :param adj_prc:
    :param date_list:
    :return:
    '''
    p0 = adj_prc.resample('D').ffill().ffill(limit=7)
    p1 = adj_prc.resample('D').ffill().ffill(limit=7)
    p0.index = p0.index + period_start
    p1.index = p1.index + period_end
    return (p1.reindex(date_list) - p0.reindex(date_list)) / p0.reindex(date_list)
    #mom1=pd.DataFrame(index=date_list, columns=adj_prc.columns)
    #for i in date_list:
    #    mom1.loc[i] = adj_prc.loc[i - period_start:i - period_end].iloc[[0, -1]].pct_change().iloc[-1]
    #return mom1
def VolumeMomentum(adj_prc,volume,date_list,period_start=DateOffset(months=1),period_end=Day()):
    ## todo two-way sorted portfolio
    '''
    文献：over the next 12 months, price momentum is more pronounced among high volume stocks
    Volume Momentum, which is buy- and- hold returns from t-6 through t-1. We limit the sample to high trading volume
        stocks, i.e., stocks in the highest quintile of average monthly trading volume measured over the past six months.

    :param adj_prc:
    :param volume:
    :param date_list:
    :return:
    '''
    return np.nan
def VolumeTrend(volume,date_list):
    ## todo 仅使用CNRDS数据，适用性不够
    '''
    volume trend, which is five-year trend in monthly trading volume scaled by average trading volume during the same
        five-year period.
    For example:
        volume=pd.read_pickle(DPath+'PVm')['Mtrdvol'].unstack()
t0=time.time()
tmp=VolumeTrend(volume,date_list)
dt=time.time()-t0
    :param volume:
    :param date_list:
    :return:
    '''
    volume_tmp=volume.reset_index()
    volume_tmp['groupby']=volume.index.year*100.0+volume.index.month
    volume_tmp_groupby=volume_tmp.groupby('groupby').nth(-1)
    volume_tmp=volume_tmp_groupby.set_index('Trdmt')
    volume_trend=pd.DataFrame(index=date_list,columns=volume_tmp.columns)
    for i in date_list:
        tmp=volume_tmp.loc[i-DateOffset(years=3)-Day(10):i-Day(1)].apply(lambda x:x/x.mean()-1.0)# -Day(10)防止最后交易日不在月末，从而造成遗漏
        X=np.arange(1,len(tmp)+1)-(len(tmp)+1.0)*0.5
        volume_trend.loc[i] =X@tmp.values/(X*X).sum()
    return volume_trend
def BetaDimson(ret_d,market_ret,date_list):
    '''
    The Dimson beta, we follow Dimson (1979) to use the lead and the lag of the market return along with the current
        market return to estimate the Dimson beta.
    For example:
        ret_d=pd.read_pickle(DPath+'PVd')['Adret'].unstack()
market_ret=pd.read_pickle(DPath+'market_ret_d_tot')
t0=time.time()
tmp=BetaDimson(ret_d,market_ret,date_list)
time.time()-t0
    :param ret_d:
    :param market_ret:
    :param date_list:
    :return:
    '''
    beta_dimson=pd.DataFrame(index=date_list,columns=ret_d.columns)
    X=pd.concat([market_ret.shift(2),market_ret.shift(1),market_ret,market_ret.shift(-1),market_ret.shift(-2)],axis=1)
    for i in date_list:
        tmpX=X.loc[i-DateOffset(years=1):i-Day()].iloc[:-2].apply(lambda x: x - x.mean())
        tmpY = ret_d.loc[tmpX.index].apply(lambda x: x - x.mean()).fillna(0)
        beta_dimson.loc[i]=(np.linalg.pinv(tmpX.values.T@tmpX.values)@tmpX.values.T@tmpY).sum()
    return beta_dimson[beta_dimson!=0.0]
def BetaDownside(ret_d,market_ret,date_list):
    '''
    Downside beta, we follow Ang, Chen, and Xing (2006) to estimate downside beta as the conditional covariance between
        a stock’s excess return and market excess return, divided by the conditional variance of market excess return,
        which is on condition that market excess return is lower than the average of market excess return.
    For example:
    :param ret_d:
    :param market_ret:
    :param date_list:
    :return:
    '''
    beta_down = pd.DataFrame(index=date_list, columns=ret_d.columns)
    for i in date_list:
        tmpX = market_ret.loc[i - DateOffset(years=1):i - Day()]
        tmpX = tmpX - tmpX.mean()
        tmpX[tmpX>0.0]=0.0 # key
        tmpY = ret_d.loc[tmpX.index].apply(lambda x: x - x.mean()).fillna(0)
        beta_down.loc[i] = tmpX.values.T @ tmpY.values/(tmpX.values.T @ tmpX.values)
    return beta_down[beta_down != 0.0]
def BetaMarket(ret_d,market_ret,date_list):
    '''
    Market beta, which is the estimated market beta from weekly returns and equal weighted market returns for 3 years
        ending month t-1 with at least 52 weeks of returns.###本代码使用的仍然是日交易数据
    :param ret_d:
    :param market_ret:
    :param date_list:
    :return:
    '''
    beta_ = pd.DataFrame(index=date_list, columns=ret_d.columns)
    for i in date_list:
        tmpX = market_ret.loc[i - DateOffset(months=1):i - Day()]
        tmpX = tmpX - tmpX.mean()
        tmpY = ret_d.loc[tmpX.index].apply(lambda x: x - x.mean()).fillna(0)
        beta_.loc[i] = tmpX.values.T @ tmpY.values/(tmpX*tmpX).sum()
    return beta_[beta_ != 0.0]
def BetaSquared(ret_d,market_ret,date_list):
    '''
    Beta squared, which is market beta squared.
    :param ret_d:
    :param market_ret:
    :param date_list:
    :return:
    '''
    tmp=BetaMarket(ret_d,market_ret,date_list)
    return tmp**2.0
def BetaFP(ret_d,market_ret,date_list):
    '''
    Fama and French (1992) beta, we follow Fama and French (1992) to estimate individual stocks’ betas by regressing
        monthly return on the current and recent lag of the market return with a five-year rolling window.
    For example:
        ret_d=pd.read_pickle(DPath+'PVd')['Adret'].unstack()*0.01
        market_ret=pd.read_pickle(DPath+'market_ret_d_tot')*0.01
t0=time.time()
tmp=BetaFP(ret_d*0.01,market_ret*0.01,date_list)
time.time()-t0 # 68.86
    :param ret_d:
    :param market_ret:
    :param date_list:
    :return:
    '''
    #ret_d_ln=np.log(ret_d+1.0)
    ret_3d=np.log(ret_d+1.0)+np.log(ret_d+1.0).shift(-1)+np.log(ret_d+1.0).shift(-2)
    #market_ret_ln=np.log(market_ret+1.0)
    market_ret_3d=np.log(market_ret+1.0)+np.log(market_ret+1.0).shift(-1)+np.log(market_ret+1.0).shift(-2)
    beta_FP=pd.DataFrame(index=date_list,columns=ret_d.columns)
    for i in date_list:
        beta_FP.loc[i]=ret_d.loc[i-DateOffset(months=3):i-Day()].std()/\
                       market_ret.loc[i-DateOffset(months=3):i-Day()].std()*\
                       ret_3d.loc[i-DateOffset(years=1):i-Day()].corrwith(market_ret_3d.loc[i-DateOffset(years=1):i-Day()])
    return beta_FP
def IdiosyncraticVolatility(ret_d,market_ret,date_list):
    '''
    Idiosyncratic return volatility, which is standard deviation of residuals of weekly returns on weekly equal
        weighted market returns for 3 years prior to month end.
    For example:
        ret_d=pd.read_pickle(DPath+'PVd')['Adret'].unstack()
market_ret=ret_d.mean(axis=1) # market_ret=pd.read_pickle('market_ret_d_eq')
t0=time.time()
tmp=IdiosyncraticVolatility(ret_d,market_ret,date_list)
time.time()-t0 # 66.45
    :param ret_d:
    :param market_ret:
    :param date_list:
    :return:
    '''
    ivol = pd.DataFrame(index=date_list, columns=ret_d.columns)
    for i in date_list:
        tmpX = market_ret.loc[i - DateOffset(months=1):i - Day()]
        tmpX=tmpX-tmpX.mean()
        tmpY = ret_d.loc[tmpX.index].apply(lambda x: x - x.mean()).fillna(0)
        tmpBeta = tmpX.values.T @ tmpY.values/(tmpX.values.T @ tmpX.values)
        ivol.loc[i]=tmpY.sub(np.matrix(tmpX).T@np.matrix(tmpBeta)).std()
    return ivol[ivol>0.0]
def Illiquidity(ret_d,tra_amount,date_list):
    '''
    Illiquidity, which is the average of absolute daily return divided by daily RMB trading volume over the past 12
        months ending on June 30 of year t+1.
    For example:
        ret_d=pd.read_pickle(DPath+'PVd')['Adret'].unstack()
tra_amount=pd.read_pickle(DPath+'PVd')['Dtrdamt'].unstack()
tmp=Illiquidity(ret_d,tra_amount,date_list)
    :param ret_d:
    :param amount:
    :param date_list:
    :return:
    '''
    illiq=pd.DataFrame(index=date_list,columns=ret_d.columns)
    illiq_d=ret_d.abs()/tra_amount*1.0e6
    for i in date_list:
        illiq.loc[i]=illiq_d.loc[i - DateOffset(months=1):i - Day()].mean()
    return illiq[illiq>0.0]
def MaxDailyReturn(ret_d,date_list):
    '''
    Maximum daily returns, which is the maximum daily return from returns during calendar month t-1.
    For example:
        ret_d=pd.read_pickle(DPath+'PVd')['Adret'].unstack()
tmp=MaxDailyReturn(ret_d,date_list)
    :param ret_d:
    :param date_list:
    :return:
    '''
    max_ret=pd.DataFrame(index=date_list,columns=ret_d.columns)
    for i in date_list:
        max_ret.loc[i]=ret_d.loc[i - DateOffset(months=1):i - Day()].max()
    return max_ret[max_ret!=0.0]
def Price(cls_prc,date_list):
    '''
    Price, which is the share price at the end of month t-1
    For example:
        cls_prc=pd.read_pickle(DPath+'PVd')['Dclsprc'].unstack()
tmp=Price(cls_prc,date_list)
    :param cls_prc:
    :param date_list:
    :return:
    '''
    return cls_prc.resample('D').ffill().shift(1).reindex(date_list)
##todo def PriceDelay():
def TradingAmount(tra_amount,date_list):
    ## todo 平均值还是last observation
    '''
    RMB trading volume, which is the natural log of RMB trading volume times price per share from month t-2
    :param tra_amount:
    :param date_list:
    :return:
    '''
    tra_amt=pd.DataFrame(index=date_list,columns=tra_amount.columns)
    tmp=np.log(tra_amount) #
    for i in date_list:
        tra_amt.loc[i]=tmp.loc[i - DateOffset(months=2):i - Day()].mean()
    return tra_amt
def Capitalization(market_cap,date_list):
    '''
    Frim size, which is market value of tradable shares at the end of each month.
    For example:
        market_cap=pd.read_pickle(DPath+'PVd')['Dmktcap'].unstack() #'Domktcap'
tmp=Size(market_cap,date_list)
    :param market_cap:
    :param date_list:
    :return:
    '''
    return np.log(market_cap[market_cap>0.0].resample('D').ffill().shift(1).reindex(date_list))
def VolatilityOfTradingAmount(tra_amount,date_list):
    '''
    Volatility of RMB trading volume, which is monthly standard deviation of daily dollar trading volume.
    For example:
        tra_amount=pd.read_pickle(DPath+'PVd')['Dtrdamt'].unstack()
tmp=VolatilityOfTradingAmount(tra_amount,date_list)
    :param tra_amount:
    :param date_list:
    :return:
    '''
    vol_tra_amt=pd.DataFrame(index=date_list,columns=tra_amount.columns)
    tmp=np.log(tra_amount)
    for i in date_list:
        vol_tra_amt.loc[i]=tmp.loc[i - DateOffset(months=1):i - Day()].std()
    return vol_tra_amt[vol_tra_amt>0.0]
def VolatilityOfTurnover(turnover,date_list):
    '''
    Volatility of turnover, which is monthly standard deviation of daily share turnover.
    For example:
        turnover=pd.read_pickle(DPath+'PVd')['Dtnor'].unstack()
tmp=VolatilityOfTurnover(turnover[turnover>0.0],date_list)
    :param turnover:
    :param date_list:
    :return:
    '''
    vol_turnover = pd.DataFrame(index=date_list, columns=turnover.columns)
    #tmp = np.log(turnover)
    for i in date_list:
        vol_turnover.loc[i] = turnover.loc[i - DateOffset(months=1):i - Day()].std()
    return vol_turnover[vol_turnover>0.0]
def VolatilityOfReturns(ret_d,date_list):
    '''
    Return volatility, which is standard deviation of daily returns from month t-1
    For example:
        ret_d=pd.read_pickle(DPath+'PVd')['Adret'].unstack();
tmp=VolatilityOfReturns(ret_d,date_list)
    :param ret_d:
    :param date_list:
    :return:
    '''
    vol_ret = pd.DataFrame(index=date_list, columns=ret_d.columns)
    for i in date_list:
        vol_ret.loc[i] = ret_d.loc[i - DateOffset(months=1):i - Day()].std()
    return vol_ret[vol_ret>0.0]
def Turnover(turnover,date_list):
    '''
    Share turnover, which is average monthly trading volume for most recent 3 months scaled by number of shares
        outstanding in current month.
    For example:
        turnover=pd.read_pickle(DPath+'PVd')['Dtnor'].unstack()
tmp=Turnover(turnover[turnover>0.0],date_list)
    :param turnover:
    :param date_list:
    :return:
    '''
    tnor=pd.DataFrame(index=date_list, columns=turnover.columns)
    for i in date_list:
        tnor.loc[i]=turnover.loc[i - DateOffset(months=2):i - Day()].mean()
    return tnor[tnor>0.0]
def CashflowToDebt(oper_cash,tot_lia,date_list,pub_date=None):
    '''
    The cash flow-to-debt ratio is the ratio of a company’s cash flow from operations to its total debt.
    现金债务总额比《偿债能力《CNRDS
    For example:
        oper_cash=CF['NCPOA'].unstack()
tot_lia=BS['LB'].unstack()
tmp1=CashflowToDebt(oper_cash,tot_lia,date_list,pub_date=pub_date)
    :param oper_cash:
    :param tot_lia:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp=2.0*oper_cash[oper_cash.index.month==12]/(tot_lia[tot_lia.index.month==12]+tot_lia[tot_lia.index.month==12].shift(1))
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def CurrentRatio(curr_assets,curr_lia,date_list,annually=True,pub_date=None):
    '''
    Current ratio, which is current assets divided by current liabilities.
    For example:
        curr_assets=BS['CA'].unstack()
curr_lia=BS['CL'].unstack()
pub_date=pd.read_pickle(DPath+'PubDate')['ActlDt']
tmp1=CurrentRatio(curr_assets,curr_lia,date_list,annually=True,pub_date=pub_date)
    :param curr_assets:
    :param curr_lia:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp=curr_assets/curr_lia
    if annually:
        tmp=tmp[tmp.index.month==12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) +Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['tmp']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['tmp'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list,columns=tmp.columns)
def CurrentRatioGrowth(curr_assets,curr_lia,date_list,annually=True,pub_date=None):
    '''
    Current ratio growth, which is annual growth in current ratio.
    For example:
        curr_assets=BS['CA'].unstack()
curr_lia=BS['CL'].unstack()
pub_date=pd.read_pickle(DPath+'PubDate')['ActlDt']
tmp1=CurrentRatioGrowth(curr_assets,curr_lia,date_list,annually=True,pub_date=pub_date)
    :param curr_assets:
    :param curr_lia:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp=curr_assets/curr_lia
    if annually:
        tmp=tmp[tmp.index.month==12].pct_change(fill_method=None)
    else:
        tmp=tmp.pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['tmp']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['tmp'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list,columns=tmp.columns)
def QuickRatio(curr_assets,curr_lia,inventory,date_list,annually=True,pub_date=None):
    '''
    Quick ratio, which is current assets minus inventory, divided by current liabilities.
    For example:
        inventory = BS['INVY'].unstack().fillna(0)
curr_assets=BS['CA'].unstack()
curr_lia=BS['CL'].unstack()
pub_date=pd.read_pickle(DPath+'PubDate')['ActlDt']
tmp1=QuickRatio(curr_assets,curr_lia,inventory,date_list,annually=True,pub_date=pub_date)
    :param curr_assets:
    :param curr_lia:
    :param inventory:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp = (curr_assets - inventory) / curr_lia
    if annually:
        tmp = tmp[tmp.index.month == 12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def QuickRatioGrowth(curr_assets,curr_lia,inventory,date_list,annually=True,pub_date=None):
    '''
    Quick ratio growth, which is annual growth in quick ratio.
    For example:
        inventory = BS['INVY'].unstack().fillna(0)
curr_assets=BS['CA'].unstack()
curr_lia=BS['CL'].unstack()
pub_date=pd.read_pickle(DPath+'PubDate')['ActlDt']
tmp1=QuickRatioGrowth(curr_assets,curr_lia,inventory,date_list,annually=True,pub_date=pub_date)
    :param curr_assets:
    :param curr_lia:
    :param inventory:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    tmp = (curr_assets - inventory) / curr_lia
    if annually:
        tmp = tmp[tmp.index.month == 12].pct_change(fill_method=None)
    else:
        tmp=tmp.pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
def SalesToCash(sales,cash_eq,date_list,pub_date=None):
    '''
    Sales-to-cash, which is sales divided by cash and cash equivalents.
    For example:
        sales=IS['REV'].unstack()
cash_eq=CF['CEEPBL'].unstack()
tmp1=SalesToCash(sales,cash_eq,date_list,pub_date=None)
    :param sales:
    :param cash_eq:
    :param date_list:
    :param pub_date:
    :return:
    '''
    tmp=sales[sales.index.month==12]/cash_eq[cash_eq.index.month==12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6)+Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2=pd.DataFrame()
        tmp2['tmp']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['tmp'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list,columns=tmp.columns)
def SalesToInventory(sales,inventory,date_list,pub_date=None):
    '''
    Sales-to-inventory, which is sales divided by total inventory.
    For example:
        sales=IS['REV'].unstack()
inventory = BS['INVY'].unstack()
tmp1=SalesToInventory(sales,inventory,date_list,pub_date=pub_date)
    :param sales:
    :param inventory:
    :param date_list:
    :param pub_date:
    :return:
    '''
    tmp = sales[sales.index.month == 12] / inventory[inventory.index.month == 12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6) + Day()
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill(limit=365).shift(1).reindex(index=date_list, columns=tmp.columns)
