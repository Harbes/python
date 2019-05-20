from my_modules.QuantTools.FactorZoo import *
import pandas as pd
import numpy as np
from pandas.tseries.offsets import MonthEnd,YearEnd,Week,Day,DateOffset
import time
#import sys
#from pandas import DataFrame,Series,qcut
#import statsmodels.api as sm
#from copy import deepcopy
#from scipy.stats import mstats
#from dateutil.parser import parse
#import warnings
#warnings.filterwarnings("ignore")
DPath='/Users/harbes/data/CNRDS/'
BS=pd.read_pickle(DPath+'BS')
PV=pd.read_pickle(DPath+'PVd')
CF=pd.read_pickle(DPath+'CF')
IS=pd.read_pickle(DPath+'IS')
##
market_cap=PV['Dmktcap'].unstack()
tot_assets=BS['AT'].unstack()
book_value=BS['EQU'].unstack()
oper_cash=CF['NCPOA'].unstack()
tot_lia=BS['LB'].unstack()
long_debt=BS['NCL'].unstack()
dividend=BS['DVP'].unstack()
annual_inc_bef_extra=IS['IBT'].unstack()-IS['TAX'].unstack()
adj_prc=PV['Adclsprc'].unstack()
sales=IS['REV'].unstack()
inventory=BS['INVY'].unstack()
tax=IS['TAX'].unstack()
tot_profit=IS['IBT'].unstack()
net_profit=IS['NI'].unstack()
PPE=BS['PPE'].unstack()
depre=CF['DFAOGA'].unstack()
fin_assets = (BS['HLTFA'].fillna(0) + BS['BBSFA'].fillna(0) + BS['AFSFA'].fillna(0)).unstack()
fin_lia = (BS['STBL'].fillna(0) + BS['HLTFL'].fillna(0) + BS['INTP'].fillna(0) + BS['DVP'].fillna(0) + BS['LTBL'].fillna(0)).unstack()
cash_eq=CF['CEEPBL'].unstack()
trad_share=PV['Domktcap'].unstack()
inc_bef_tax=IS['IBT'].unstack()
fin_exp=IS['FINEXP'].unstack()
net_working=(BS['CA']-BS['CL']).unstack()
net_fixed=BS['PPE'].unstack()-CF['DFAOGA'].unstack().fillna(0)
oper_rev=IS['OPREV'].unstack()
oper_exp=IS['COGS'].fillna(0).unstack()
net_inc=IS['NI'].unstack()
oper_inc=IS['OI'].unstack()
non_oper_inc=IS['NOINC'].unstack()
ebit=IS['IBT'].unstack()-IS['FINEXP'].unstack().fillna(0)
sector=pd.read_pickle(DPath+'sector_datetime')
sec_sign=np.arange(29.0)
volume=pd.read_pickle(DPath+'PVm')['Mtrdvol'].unstack()
ret_d=pd.read_pickle(DPath+'PVd')['Adret'].unstack()
market_ret=pd.read_pickle(DPath+'market_ret_d_tot')
tra_amount=PV['Dtrdamt'].unstack()
cls_prc=PV['Dclsprc'].unstack()
turnover=PV['Dtnor'].unstack()
curr_assets=BS['CA'].unstack()
curr_lia=BS['CL'].unstack()
retained_earnings=BS['SR'].unstack().fillna(0)+BS['UDP'].unstack().fillna(0)
del PV
del BS
del IS
del CF
##
date_list=pd.date_range('2004-12-31','2018-9-1',freq='M')+Day()
#pub_date=pd.read_pickle(DPath+'PubDate')['ActlDt']

##指标计算
t0=time.time()
AM=AssetsToMarket(tot_assets,market_cap,date_list,annually=True,pub_date=None,most_recent=False)
BM=BookToMarket(book_value,market_cap,date_list,annually=True,pub_date=None)
CFP=OperatingCashFlowToPrice(oper_cash,market_cap,date_list,pub_date=None)
DER=DebtToEquity(tot_lia,market_cap,date_list,annually=True,pub_date=None)
DLME=LongTermDebtToMarketEquity(long_debt,market_cap,date_list,annually=True,pub_date=None)
DP=DividendToPrice(dividend,market_cap,date_list,pub_date=None)
EP=EarningsToPrice(annual_inc_bef_extra,market_cap,date_list,pub_date=None)
LG=LiabilityGrowth(tot_lia,date_list,annually=True,pub_date=None)
PY=PayoutYield(annual_inc_bef_extra,book_value,market_cap,date_list,pub_date=None)
Rev1=Reversal_60_13(adj_prc,date_list)
SG=SustainableGrowth(book_value,date_list,annually=True,pub_date=None)
SMI=SalesGrowthMinusInventoryGrowth(sales,inventory,date_list,pub_date=None)
SP=SalesToPrice(sales,market_cap,date_list,pub_date=None)
TG=TaxGrowth(tax,date_list,pub_date=None)
ACC=Accruals(annual_inc_bef_extra,oper_cash,tot_assets,date_list,pub_date=None)
PACC=PercentAccruals(tot_profit,oper_cash,net_profit,date_list,pub_date=None)
CAPXG=CapitalExpenditureGrowth(PPE,depre,date_list,pub_date=None)
dBe=ChangeInShareholdersEquity(book_value,tot_assets,date_list,annually=True,pub_date=None)
dPIA=ChangesInPPEandInventoryToAssets(PPE,inventory.fillna(0),tot_assets,date_list,annually=True,pub_date=None)
IA=InvestmentToAssets(tot_assets,date_list,annually=True,pub_date=None)
IVC=InventoryChange(inventory,tot_assets,date_list,annually=True,pub_date=None)
IVG=InventoryGrowth(inventory,date_list,annually=True,pub_date=None)
NOA=NetOperatingAssets(fin_assets,fin_lia,book_value,tot_assets,date_list,annually=True,pub_date=None)
ATO=AssetTurnover(fin_assets,fin_lia,book_value,sales,date_list,pub_date=None)
CFOA=CashFlowOverAssets(oper_cash,tot_assets,date_list,pub_date=None)
CP=CashProductivity(trad_share,long_debt,tot_assets,cash_eq,date_list,annually=True,pub_date=None)
CTA=CashToAssets(cash_eq,tot_assets,date_list,annually=True,pub_date=None)
CTO=CapitalTurnover(sales,tot_assets,date_list,pub_date=None)
EBIT=ReturnOnCapital(inc_bef_tax,fin_exp,net_working,net_fixed,date_list,pub_date=None)
EY=EarningsYield(inc_bef_tax,fin_exp,market_cap,tot_lia,cash_eq,date_list)
GM=GrossMargins(oper_rev,oper_exp,date_list,pub_date=None)
GP=GrossProfitability(oper_rev,oper_exp,tot_assets,date_list,pub_date=None)
NPOP=NetPayoutOverProfit(net_inc,book_value,tot_profit,date_list,pub_date=None)
RNA=ReturnOnOperatingAsset(oper_inc,depre,fin_assets,fin_lia,date_list,pub_date=None)
ROA=ReturnOnAssets(oper_inc,tot_assets,date_list,pub_date=None)
ROE=ReturnOnEquity(net_inc,book_value,date_list,pub_date=None)
ROIC=ReturnOnInvestedCapital(inc_bef_tax,fin_exp,non_oper_inc.fillna(0),market_cap,tot_lia,cash_eq,date_list)
TBI=TexableIncomeToBookIncome(inc_bef_tax,net_inc,date_list,pub_date=None)
Z_score=ZScore(net_working,tot_assets,retained_earnings,ebit,market_cap,tot_lia,sales,date_list)
CHMOM=ChangeIn6MonthMomentum(adj_prc,date_list)
INDMOM=IndustryMomentum(adj_prc,sector,sec_sign,date_list,period_start=DateOffset(months=3),period_end=Day())
MOM1M=Momentum(adj_prc,date_list,period_start=Day(30),period_end=Day())
MOM6M=Momentum(adj_prc,date_list,period_start=Day(183),period_end=Day(30))
MOM12M=Momentum(adj_prc,date_list,period_start=Day(365),period_end=Day(30))
MOM36M=Momentum(adj_prc,date_list,period_start=Day(365*3),period_end=Day(365))
VOLT=VolumeTrend(volume,date_list)
B_DIM=BetaDimson(ret_d,market_ret,date_list)
B_DN=BetaDownside(ret_d,market_ret,date_list)
BETA=BetaMarket(ret_d,market_ret,date_list)
BETASQ=BetaSquared(ret_d,market_ret,date_list)
B_FP=BetaFP(ret_d*0.01,market_ret*0.01,date_list)
IVOL=IdiosyncraticVolatility(ret_d,market_ret,date_list)
ILLIQ=Illiquidity(ret_d,tra_amount,date_list)
MAXRET=MaxDailyReturn(ret_d,date_list)
PRC=Price(cls_prc,date_list).log()
RVOL=TradingAmount(tra_amount,date_list).log()
SIZE=Capitalization(market_cap,date_list)
STD_RVOL=VolatilityOfTradingAmount(tra_amount,date_list).log() #＃TODO 调整量纲还是用对数函数
STD_TURN=VolatilityOfTurnover(turnover,date_list)
RETVOL=VolatilityOfReturns(ret_d,date_list)
TURN=Turnover(turnover,date_list)
CFD=CashflowToDebt(oper_cash,tot_lia,date_list,pub_date=None)
CR=CurrentRatio(curr_assets,curr_lia,date_list,annually=True,pub_date=None)
CRG=CurrentRatioGrowth(curr_assets,curr_lia,date_list,annually=True,pub_date=None)
QR=QuickRatio(curr_assets,curr_lia,inventory,date_list,annually=True,pub_date=None)
QRG=QuickRatioGrowth(curr_assets,curr_lia,inventory,date_list,annually=True,pub_date=None)
SC=SalesToCash(sales,cash_eq,date_list,pub_date=None)
SI=SalesToInventory(sales,inventory,date_list,pub_date=None)
t1=time.time()-t0
##
signal_list=[AM,BM,CFP,DER,DLME,DP,EP,LG,PY,SG,SMI,SP,TG,ACC,PACC,CAPXG,dBe,dPIA,IA,IVC,IVG,NOA,ATO,CFOA,CP,CTA,CTO,EBIT,
             EY,GM,GP,NPOP,RNA,ROA,ROE,ROIC,TBI,Z_score,
             CHMOM,INDMOM,MOM1M,MOM6M,MOM12M,MOM36M,VOLT,
             B_DIM,B_DN,BETA,BETASQ,B_FP,
             IVOL,ILLIQ,MAXRET,PRC,RVOL,SIZE,STD_RVOL,STD_TURN,RETVOL,TURN,
             CFD,CR,CRG,QR,QRG,SC,SI]

