#!/usr/bin/env python

""" run GDAY, plot LAI, GPP, shoot NC """

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import plot_settings as ps
import pandas as pd
import datetime as dt

__author__  = "Martin De Kauwe"
__version__ = "1.0 (16.05.2013)"
__email__   = "mdekauwe@gmail.com"

def date_converter(*args): 
    return dt.datetime.strptime(str(int(float(args[0]))) + " " +\
                                str(int(float(args[1]))), '%Y %j')

def read_data(fname):
    df = pd.read_csv(fname, parse_dates=[[0,1]], index_col=0, sep=",", 
                     keep_date_col=True, date_parser=date_converter, 
                     na_values=["-999.9"], skiprows=2)
    return df
    

def main(site, experiment):
    # run new simulations
    os.system("eucface_simulations.py")

    # load data
    amb = read_data("../outputs/D1GDAY%sAMB%s.csv" % (site, experiment))
    ele = read_data("../outputs/D1GDAY%sELE%s.csv" % (site, experiment))

    
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['legend.fontsize'] = 9
    fig = plt.figure()



    ax = fig.add_subplot(311)
    ax.plot(amb["NCAN"]/amb["CL"], "b-", label="Amb")
    ax.plot(ele["NCAN"]/ele["CL"], "r-", label="Ele")
    ax.legend(numpoints=1, loc="best")
    ax.set_ylabel("Shoot N:C")
    ax.set_ylim(0, 0.05)
    
    ax = fig.add_subplot(312)
    ax.plot(amb["GPP"], "r-", label="Amb")
    ax.plot(ele["GPP"], "g-", label="Ele")
    ax.legend(numpoints=1, loc="best")
    ax.set_ylabel("GPP")
    #ax.set_ylim(0, 16)

    ax = fig.add_subplot(313)
    ax.plot(amb["LAI"], "r-", label="Amb")
    ax.plot(ele["LAI"], "g-", label="Ele")
    ax.legend(numpoints=1, loc="best")
    ax.set_ylabel("LAI")
   

    plt.show()
    
if __name__ == "__main__":
    
    site = "EUC"
    experiment = "VAR"
    main(site, experiment)
