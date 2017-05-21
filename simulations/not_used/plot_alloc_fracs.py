#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np

def date_converter(*args): 
    s = str(int(float(args[0]))) + " " + str(int(float(args[1])))
    return dt.datetime.strptime(s, '%Y %j')

def build_df(fname):
    return pd.read_csv(fname, parse_dates=[[0,1]], skiprows=3, index_col=0, 
                       sep=",", keep_date_col=True, 
                       date_parser=date_converter)

def get_alloc_fracs(df):
    total_growth = (df.groupby("YEAR").GL.sum() + \
                    df.groupby("YEAR").GW.sum() + \
                    df.groupby("YEAR").GCR.sum() + \
                    df.groupby("YEAR").GR.sum())
    
    
    
    #af = df["GL"] / df["NPP"]
    #aw = df["GW"] / df["NPP"]
    #ar = df["GR"] / df["NPP"]
    
    af = df.groupby("YEAR").GL.sum() / total_growth
    aw = df.groupby("YEAR").GW.sum() / total_growth
    afr = df.groupby("YEAR").GR.sum() / total_growth
    acr = (df.groupby("YEAR").GCR.sum() ) / total_growth
    
    return af, aw, afr, acr
                
fname = "../outputs/D1GDAYEUCAMBVAR.csv"
df_amb = build_df(fname)
(af_amb, aw_amb, afr_amb, acr_amb) = get_alloc_fracs(df_amb)

fname = "../outputs/D1GDAYEUCELEVAR.csv"
df_ele = build_df(fname)
(af_ele, aw_ele, afr_ele, acr_ele) = get_alloc_fracs(df_ele)

years = np.unique(df_amb.YEAR)

plt.subplot(111)
plt.plot(years, af_amb, "b-", label="Af")
plt.plot(years, aw_amb, "g-", label="Aw")
plt.plot(years, afr_amb, "r-", label="Afr")
plt.plot(years, acr_amb, "y-", label="Acr")
plt.plot(years, af_ele, "b--")
plt.plot(years, aw_ele, "g--")
plt.plot(years, afr_ele, "r--")
plt.plot(years, acr_ele, "y--")
plt.legend(numpoints=1, loc="best")
plt.ylim(0, 1)
plt.show()
                  
