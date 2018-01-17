import numpy as np
import pandas as pd
from pandas.tseries.offsets import MonthBegin
from datetime import datetime


def import_data():
    global data1, data2, price0, price1, size
    # data_path='F:/data/xccdata'
    data_path = '/Users/harbes/data/xccdata'
    data1 = pd.read_pickle(data_path + '/PV_datetime')[['adj_open', 'adj_close', 'size_tot']]
    data2 = pd.read_pickle(data_path + '/BS')[['fin_year', 'stkcd', 'tot_assets', 'tot_liab']].set_index(
        ['fin_year', 'stkcd']).sort_index()
    # R_free


def process_and_calculate():
    global size, BM
    size = data1['size_tot'].unstack()
    book = (data2['tot_assets'] - data2['tot_liab']).unstack()
    book.index = book.index.astype(str).to_datetime('%Y%m%d')
    book.index.name = 'trddt'
    # data.to_pickle('./five_factor/book')
    book = book.resample('D').first().ffill()
    book.index = book.index + MonthBegin()


if __name__ == '__main__':
    import data

    ()
