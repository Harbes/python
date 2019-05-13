# todo 几点建议：.loc换成.reindex;将SustainableGrowth函数作为模板推广至其他函数；财报数据出现加减时，引入fillna
import sys
import pandas as pd
from pandas import DataFrame,Series,qcut
import numpy as np
from pandas.tseries.offsets import MonthEnd,YearEnd,Week,Day,DateOffset
import statsmodels.api as sm
from copy import deepcopy
from scipy.stats import mstats
from dateutil.parser import parse
import time
#import warnings
#warnings.filterwarnings("ignore")
DPath='/Users/harbes/data/CNRDS/'
BS=pd.read_pickle(DPath+'BS')
#BS=BS.iloc[:,:-1].astype(float)
#BS.iloc[0,0]
#BS.to_pickle(DPath+'BS')
PV=pd.read_pickle(DPath+'PV')
#PV=PV.astype(float)
#PV.iloc[0,0]
#PV.to_pickle(DPath+'PV')
CF=pd.read_pickle(DPath+'CF')
IS=pd.read_pickle(DPath+'IS')
annual_inc_bef_extra=IS['IBT'].unstack()-IS['TAX'].unstack()
market_cap=PV['Dmktcap'].unstack()
date_list=annual_inc_bef_extra.loc['2004':].index[1:20]
tmp=EarningsToPrice(annual_inc_bef_extra,market_cap,date_list)
market_cap.index[market_cap.index.month==12]

tot_assets=BS['AT'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=tot_assets.index[10:]
        tmp21=AssetsToMarket(tot_assets,market_cap,date_list)
def AssetsToMarket(total_assets,market_cap,date_list,most_recent=False,annually=True):
    '''
    Attention:数据时间index影响指标是否是instantly updated
    Assets-to-market, which is total assets for the fiscal year divided by fiscal-year-end market capitalization.
    for example:
        BS=pd.read_pickle(DPath+'BS')
        PV=pd.read_pickle(DPath+'PV')
        tot_assets=BS['AT'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=tot_assets.index[10:]
        tmp3=AssetsToMarket(tot_assets,market_cap,date_list,most_recent=True)
    :param tot_assets:
    :param market_cap:
    :param date_list:
    :return:
    '''
    common_columns=total_assets.columns&market_cap.columns
    if most_recent:
        return total_assets.resample('D').ffill().loc[(date_list, common_columns)].div(
            market_cap.resample('D').ffill().loc[(date_list, common_columns)])
    elif annually:
        #GroupByFunc = lambda x: x.year
        #tmp = market_cap.groupby(GroupByFunc).last()
        #tmp.index = pd.to_datetime(tmp.index, format='%Y') + YearEnd()
        #tmp2 = total_assets.loc[(total_assets.index.month == 12, common_columns)] / tmp[common_columns]
        tmp = total_assets.loc[(total_assets.index.month == 12, common_columns)]
        tmp2 = tmp / market_cap[common_columns].resample('D').ffill().loc[tmp.index]
        tmp2.index = tmp2.index + MonthEnd(6)
        return tmp2.resample('D').ffill().loc[date_list]
    else:
        tmp3= total_assets[common_columns] / market_cap[common_columns].resample('D').ffill().loc[total_assets.index]
        tmp3.index = tmp3.index + MonthEnd(6)
        return tmp3.resample('D').ffill().loc[date_list]
def BookToMarket(book_value,market_cap,date_list,most_recent=False,annually=True):
    '''
    Book-to-market equity, which is the book value of equity for fiscal year divided by fiscal-year-end market capitalization
    For example:
        BS=pd.read_pickle(DPath+'BS')
        PV=pd.read_pickle(DPath+'PV')
        book_value=BS['EQU'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=book_value.loc['2004':].index[1:20]
        tmp=BookToMarket(book_value,market_cap,date_list)

    :param book_value:
    :param market_cap:
    :param date_list:
    :return:
    '''
    common_columns=book_value.columns&market_cap.columns
    if most_recent:
        return book_value.resample('D').ffill().loc[(date_list, common_columns)].div(
            market_cap.resample('D').ffill().loc[(date_list, common_columns)])
    elif annually:
        #GroupByFunc=lambda x:x.year
        #tmp=market_cap.groupby(GroupByFunc).last()
        #tmp.index=pd.to_datetime(tmp.index,format='%Y')+YearEnd()
        #tmp2=book_value.loc[(book_value.index.month==12,common_columns)]/tmp[common_columns]
        tmp = book_value.loc[(book_value.index.month == 12, common_columns)]
        tmp2 = tmp / market_cap[common_columns].resample('D').ffill().loc[tmp.index]
        tmp2.index=tmp2.index+MonthEnd(6)
        return tmp2.resample('D').ffill().loc[date_list]
    else:
        tmp3 = book_value[common_columns] / market_cap[common_columns].resample('D').ffill().loc[book_value.index]
        tmp3.index = tmp3.index + MonthEnd(6)
        return tmp3.resample('D').ffill().loc[date_list]
def OperatingCashFlowToPrice(oper_cash,market_cap,date_list,most_recent=False,annually=True):
    '''
    Cash flow-to-price, which is operating cash flows divided by fiscal-year-end market capitalization.
    For example:
        PV=pd.read_pickle(DPath+'PV')
        CF=pd.read_pickle(DPath+'CF')
        oper_cash=CF['NCPOA'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=oper_cash.loc['2004':].index[1:20]
        tmp=OperatingCashFlowToPrice(oper_cash,market_cap,date_list)

    :param cash_flow:
    :param market_cap:
    :param date_list:
    :return:
    '''

    common_columns=oper_cash.columns&market_cap.columns
    if most_recent:
        return oper_cash.resample('D').ffill().loc[(date_list, common_columns)].div(
            market_cap.resample('D').ffill().loc[(date_list, common_columns)])
    elif annually:
        tmp = oper_cash.loc[(oper_cash.index.month == 12, common_columns)]
        tmp2 = tmp / market_cap[common_columns].resample('D').ffill().loc[tmp.index]
        tmp2.index=tmp2.index+MonthEnd(6)
        return tmp2.resample('D').ffill().loc[date_list]
    else:
        tmp3 = oper_cash[common_columns] / market_cap[common_columns].resample('D').ffill().loc[oper_cash.index]
        tmp3.index = tmp3.index + MonthEnd(6)
        return tmp3.resample('D').ffill().loc[date_list]
def DebtToEquity(total_lia,market_cap,date_list,most_recent=False,annually=True):
    '''
    Debt-to-equity ratio, which is total liabilities divided by fiscal-year-end market capitalization.
    For example:
        BS=pd.read_pickle(DPath+'BS')
        PV=pd.read_pickle(DPath+'PV')
        total_lia=BS['LB'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=total_lia.loc['2004':].index[1:20]
        tmp=DebtToEquity(total_lia,market_cap,date_list)
    :param total_lia:
    :param market_cap:
    :param date_list:
    :return:
    '''
    common_columns = total_lia.columns & market_cap.columns
    if most_recent:
        return total_lia.resample('D').ffill().loc[(date_list, common_columns)].div(
            market_cap.resample('D').ffill().loc[(date_list, common_columns)])
    elif annually:
        tmp = total_lia.loc[(total_lia.index.month == 12, common_columns)]
        tmp2 = tmp / market_cap[common_columns].resample('D').ffill().loc[tmp.index]
        tmp2.index=tmp2.index+MonthEnd(6)
        return tmp2.resample('D').ffill().loc[date_list]
    else:
        tmp3 = total_lia[common_columns] / market_cap[common_columns].resample('D').ffill().loc[total_lia.index]
        tmp3.index = tmp3.index + MonthEnd(6)
        return tmp3.resample('D').ffill().loc[date_list]
def LongTermDebtToMarketEquity(long_debt,market_cap,date_list,most_recent=False,annually=True):
    '''
    Long term debt-to-market equity, which is long term liabilities divided by fiscal-year-end market capitalization.
    For example:
        long_debt=BS['NCL'].unstack()
        market_cap=PV['Dmktcap'].unstack()
        date_list=long_debt.loc['2004':].index[1:20]
        tmp=LongTermDebtToMarketEquity(long_debt,market_cap,date_list)
    :param long_debt:
    :param market_cap:
    :param date_list:
    :return:
    '''
    common_columns = long_debt.columns & market_cap.columns
    if most_recent:
        return long_debt.resample('D').ffill().loc[(date_list, common_columns)].div(
            market_cap.resample('D').ffill().loc[(date_list, common_columns)])
    elif annually:
        tmp = long_debt.loc[(long_debt.index.month == 12, common_columns)]
        tmp2 = tmp / market_cap[common_columns].resample('D').ffill().loc[tmp.index]
        tmp2.index=tmp2.index+MonthEnd(6)
        return tmp2.resample('D').ffill().loc[date_list]
    else:
        tmp3 = long_debt[common_columns] / market_cap[common_columns].resample('D').ffill().loc[long_debt.index]
        tmp3.index = tmp3.index + MonthEnd(6)
        return tmp3.resample('D').ffill().loc[date_list]
def DividendToPrice(annual_div,market_cap,date_list,most_recent=False,annually=True):
    '''
    Dividend-to-price ratio, which is annual total dividends payouts divided by fiscal-year-end market capitalization.

    :param dividend:
    :param market_cap:
    :param date_list:
    :return:
    '''
    common_columns = annual_div.columns & market_cap.columns
    if most_recent:
        return annual_div.resample('D').ffill().loc[(date_list, common_columns)].div(
            market_cap.resample('D').ffill().loc[(date_list, common_columns)])
    elif annually:
        tmp = annual_div.loc[(annual_div.index.month == 12, common_columns)]
        tmp2 = tmp / market_cap[common_columns].resample('D').ffill().loc[tmp.index]
        tmp2.index = tmp2.index + MonthEnd(6)
        return tmp2.resample('D').ffill().loc[date_list]
    else:
        tmp3 = annual_div[common_columns] / market_cap[common_columns].resample('D').ffill().loc[annual_div.index]
        tmp3.index = tmp3.index + MonthEnd(6)
        return tmp3.resample('D').ffill().loc[date_list]
def EarningsToPrice(annual_inc_bef_extra,market_cap,date_list,most_recent=False,annually=True):
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
    common_columns = annual_inc_bef_extra.columns & market_cap.columns
    if most_recent:
        return annual_inc_bef_extra.resample('D').ffill().loc[(date_list, common_columns)].div(
            market_cap.resample('D').ffill().loc[(date_list, common_columns)])
    elif annually:
        tmp = annual_inc_bef_extra.loc[(annual_inc_bef_extra.index.month == 12, common_columns)]
        tmp2 = tmp / market_cap[common_columns].resample('D').ffill().loc[tmp.index]
        tmp2.index = tmp2.index + MonthEnd(6)
        return tmp2.resample('D').ffill().loc[date_list]
    else:
        tmp3 = annual_inc_bef_extra[common_columns] / market_cap[common_columns].resample('D').ffill().loc[annual_inc_bef_extra.index]
        tmp3.index = tmp3.index + MonthEnd(6)
        return tmp3.resample('D').ffill().loc[date_list]
def LiabilityGrowth(total_lia,date_list,most_recent=False,annually=True):
    '''
    Liability growth, which is the annual change in total liabilities divided by 1-year-lagged total liabilities.
    For example:
        total_lia=BS['LB'].unstack()
        tmp2=LiabilityGrowth(total_lia,total_lia.index[15:],annually=False)
    :param total_lia:
    :param date_list:
    :return:
    '''
    if most_recent:
        tmp=total_lia.resample('D').ffill().pct_change(fill_method=None)
        tmp[tmp==0.0]=np.nan
        return tmp.loc[date_list]
    elif annually:
        tmp=total_lia[total_lia.index.month==12].pct_change(fill_method=None)
        tmp.index=tmp.index+MonthEnd(6)
        return tmp.resample('D').ffill().ffill().loc[date_list]
    else:
        tmp = total_lia.pct_change(fill_method=None)
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
def AdjOperatingCashFLowToPrice(adj_operating,market_cap,date_list):
    ##todo adj_opeating???
    '''
    Operating cash flow-to-price, which is operating cash flow(measured as earnings adjusted for depreciation and
        working capital accruals) divided by fiscal-year-end market capitalization.
    :param adj_operating:
    :param market_cap:
    :param date_list:
    :return:
    '''
    common_columns = adj_operating.columns & market_cap.columns
    return adj_operating.resample('D').ffill().loc[(date_list, common_columns)].div(
        market_cap.resample('D').ffill().loc[(date_list, common_columns)])
def PayoutYield(annual_inc_bef_extra,book_value,market_cap,date_list):
    ##todo
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
    common_columns = annual_inc_bef_extra.columns& book_value.columns & market_cap.columns
    annual_inc_minus_change_book=annual_inc_bef_extra[annual_inc_bef_extra.index.month==12]-\
                                 book_value[book_value.index.month==12]+\
                                 book_value[book_value.index.month==12].shift(1)
    tmp=annual_inc_minus_change_book[common_columns]/market_cap[common_columns].resample('D').ffill().ffill().loc[annual_inc_minus_change_book.index]
    tmp.index=tmp.index+MonthEnd(6)
    return tmp.resample('D').ffill().ffill().loc[date_list]
def Reversal_60_13(adj_price,date_list):
    ##TODO
    '''
    Reversal, which is cumulative returns from months t-60 to t-13.
    For example:
        adj_price=PV['Adclsprc'].unstack()
        tmp=Reversal_60_13(adj_price,date_list)
    :param adj_price:
    :param date_list:
    :return:
    '''
    rev=pd.DataFrame(index=date_list,columns=adj_price.columns)
    for i in date_list:
        tmp=adj_price.loc[i-DateOffset(months=60):i-DateOffset(months=12)]
        rev.loc[i]=(tmp.iloc[-1]-tmp.iloc[0])/tmp.iloc[0]
    return rev

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
        tmp=book_value.pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2=pd.DataFrame()
        tmp2['growth']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['growth'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(date_list)
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
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['growth'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['growth'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(date_list)
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
        tmp2=tmp / market_cap.resample('D').ffill().ffill().reindex(tmp.index)
        tmp2.index = tmp2.index + MonthEnd(6)
        return tmp2.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['sales'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['sales'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(date_list)/market_cap.resample('D').ffill().reindex(date_list)
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
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['growth'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['growth'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(date_list)
def Accruals(accruals,oper_cash,tot_assets,date_list,pub_date=None):
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
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['accruals'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['accruals'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(date_list)
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
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['per'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['per'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(date_list)
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
    tmp=(PPE[PPE.index.month==12]-PPE[PPE.index.month==12].shift(1)+depre[depre.index.month==12]).pct_change(fill_method=None)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().reindex(date_list)
    else:
        tmp2 = pd.DataFrame()
        tmp2['growth'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['growth'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(date_list)
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
        tmp=(book_value.loc[book_value.index.month==12]-book_value[book_value.index.month==12].shift(1))/\
            tot_assets[tot_assets.index.month==12].shift(1)
    else:
        tmp=(book_value-book_value.shift(1))/tot_assets.shift(1)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(date_list)
def ChangesInPPEandInventoryToAssets(PPE,inventory,tot_assets,date_list,annually=True,pub_date=None):
    ## todo 对于没有inventory(例如银行)的数据，舍弃(nan)还是取值为0？？？
    '''
    Changes in PPE and inventory-to-assets, which is the annual change in gross property, plant, and equipment
        plus the annual change in inventory scaled by 1-year-lagged total assets.
    For example:
        PPE=BS['PPE'].unstack()
        inventory=BS['INVY'].unstack()
        tot_assets=BS['AT'].unstack()
        tmp2=ChangesInPPEandInventoryToAssets(PPE,inventory,tot_assets,date_list,annually=False,pub_date=pub_date)
    :param PPE:
    :param inventory:
    :param tot_assets:
    :param date_list:
    :param annually:
    :param pub_date:
    :return:
    '''
    if annually:
        tmp=((PPE+inventory).loc[PPE.index.month==12]-(PPE+inventory)[PPE.index.month==12].shift(1))/\
            tot_assets[tot_assets.index.month==12].shift(1)
    else:
        tmp=((PPE+inventory)-(PPE+inventory).shift(1))/tot_assets.shift(1)
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(index=date_list,columns=tmp.columns)
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
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(index=date_list,columns=tmp.columns)

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
        tmp=2.0*(inventory[inventory.index.month==12]-inventory[inventory.index.month==12].shift(1))/\
            (tot_assets[tot_assets.index.month==12]+tot_assets[tot_assets.index.month==12].shift(1))
    else:
        tmp=2.0*(inventory-inventory.shift(1))/(tot_assets+tot_assets.shift(1))
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(index=date_list,columns=tmp.columns)
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
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(index=date_list,columns=tmp.columns)
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
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2=pd.DataFrame()
        tmp2['change']=tmp.stack()
        tmp2['date']=pub_date
        tmp2=tmp2.reset_index()
        tmp2=tmp2.loc[tmp2['date'].notnull()].set_index(['date','Scode'])['change'].sort_index()
        tmp2=tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(index=date_list,columns=tmp.columns)
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
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2 = pd.DataFrame()
        tmp2['change'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['change'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(index=date_list, columns=tmp.columns)
def CashFlowOverAssets(oper_cash,tot_assets,date_list,annually=True,pub_date=None):
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
    tmp = oper_cash/tot_assets
    if annually:
        tmp = tmp[tmp.index.month == 12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(index=date_list, columns=tmp.columns)
def CashProductivity(trad_share,long_debt,tot_assets,cash_eq,date_list,annually=True,pub_date=None):
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
    tmp = (trad_share.resample('D').ffill().ffill().loc[long_debt.index]+long_debt-tot_assets)/cash_eq
    if annually:
        tmp = tmp[tmp.index.month == 12]
    if pub_date is None:
        tmp.index = tmp.index + MonthEnd(6)
        return tmp.resample('D').ffill().loc[date_list]
    else:
        tmp2 = pd.DataFrame()
        tmp2['tmp'] = tmp.stack()
        tmp2['date'] = pub_date
        tmp2 = tmp2.reset_index()
        tmp2 = tmp2.loc[tmp2['date'].notnull()].set_index(['date', 'Scode'])['tmp'].sort_index()
        tmp2 = tmp2.loc[~tmp2.index.duplicated()].unstack()
        return tmp2.resample('D').ffill().ffill().reindex(index=date_list, columns=tmp.columns)
