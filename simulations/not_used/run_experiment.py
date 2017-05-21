#!/usr/bin/env python

""" run GDAY, plot LAI, GPP, shoot NC """

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import plot_settings as ps
import subprocess
import pandas as pd
import datetime as dt

__author__  = "Martin De Kauwe"
__version__ = "1.0 (16.05.2013)"
__email__   = "mdekauwe@gmail.com"


def main(SEND_TO_DROPBOX=False):
    os.system("eucface_simulations.py")

    # load data
    amb = read_data("../outputs/D1GDAYEUCAMBVAR.csv")
    ele = read_data("../outputs/D1GDAYEUCELEVAR.csv")

    plt.rcParams['figure.subplot.hspace'] = 0.15
    plt.rcParams['figure.subplot.wspace'] = 0.15
    plt.rcParams['font.size'] = 10
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 10.0
    plt.rcParams['ytick.labelsize'] = 10.0
    plt.rcParams['axes.labelsize'] = 10.0
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.style'] = "normal"
    plt.rcParams['font.serif'] = "Helvetica"
    fig = plt.figure()

    years = np.unique(amb.YEAR)


    ax1 = fig.add_subplot(311)

    ax1.plot(years, amb.groupby("YEAR").GPP.sum(), "b-", label="Amb")
    ax1.plot(years, ele.groupby("YEAR").GPP.sum(), "r-", label="Ele")
    ax1.legend(numpoints=1, loc="best")
    ax1.set_ylabel("GPP")
    plt.setp(ax1.get_xticklabels(), visible=False)
    #ax1.set_ylim(0., 2000)
    ax1.set_xlim(float(years[0]), float(years[-1]))

    ax2 = fig.add_subplot(312)
    ax2.plot(years, amb.groupby("YEAR").LAI.max(), "b-", label="Amb")
    ax2.plot(years, ele.groupby("YEAR").LAI.max(), "r-", label="Ele")
    ax2.set_ylabel("LAI")
    ax2.set_xlabel("Year")
    #ax2.set_ylim(0., 8)
    ax2.set_xlim(int(float(years[0])), int(float(years[-1])))

    shoot_amb = amb.groupby("YEAR").CL.max()
    shoot_ele = ele.groupby("YEAR").CL.max()
    shootn_amb = amb.groupby("YEAR").NCAN.max()
    shootn_ele = ele.groupby("YEAR").NCAN.max()

    gpp_amb = amb.groupby("YEAR").GPP.sum()
    gpp_ele = ele.groupby("YEAR").GPP.sum()

    trans_amb = amb.groupby("YEAR").T.sum()
    trans_ele = ele.groupby("YEAR").T.sum()

    wue_amb = gpp_amb / trans_amb
    wue_ele = gpp_ele / trans_ele
    
    co2_amb = amb.groupby("YEAR").CO2.mean()
    co2_ele = ele.groupby("YEAR").CO2.mean()
    
    gpp_response = ((gpp_ele / gpp_amb) - 1.0) * 100.0
    trans_response = ((trans_ele / trans_amb) - 1.0) * 100.0
    wue_response = ((wue_ele / wue_amb) - 1.0) * 100.0
    co2_response = ((co2_ele / co2_amb) - 1.0) * 100.0

    ax3 = fig.add_subplot(313)
    ax3.plot(years, gpp_response, "r-", label="GPP")
    ax3.plot(years, trans_response, "g-", label="T")
    ax3.plot(years, wue_response, "b-",label="WUE")

    ax3.legend(numpoints=1, loc="best")
    ax3.set_ylabel("CO$_2$ response")
    ax3.set_xlabel("Year")
    ax3.set_ylim(-30., 50)
    ax3.set_xlim(float(years[0]), float(years[-1]))
    plt.show()

    print "GPP: %.2f +/- %.2f" % (gpp_response.mean(), gpp_response.std())
    print "T: %.2f +/- %.2f" % (trans_response.mean(), trans_response.std())
    print "WUE: %.2f +/- %.2f" % (wue_response.mean(), wue_response.std())
    print "CO2: %.2f +/- %.2f" % (co2_response.mean(), co2_response.std())
    print
    print "leaf NC amb: %.5f +/- %.5f; leaf NC ele: %.5f +/- %.5f" % (np.mean(shootn_amb/shoot_amb), np.std(shootn_amb/shoot_amb, ddof=1), np.mean(shootn_ele/shoot_ele), np.std(shootn_ele/shoot_ele, ddof=1))

   
    if SEND_TO_DROPBOX:
        cwd_path = os.getcwd()
        path = "../outputs"
        os.chdir(path)
        subprocess.call("send_files_to_dropbox.BASH")
        print
        print "...files sent to dropbox..."
        os.chdir(cwd_path)
    
    
def date_converter(*args): 
    return dt.datetime.strptime(str(int(float(args[0]))) + " " +\
                                str(int(float(args[1]))), '%Y %j')

def read_data(fname):
    return pd.read_csv(fname, parse_dates=[[0,1]], index_col=0, sep=",", 
                    keep_date_col=True, date_parser=date_converter, 
                    na_values=["-999.9"], skiprows=3)

        
if __name__ == '__main__':
    
    user_arg = sys.argv[1:]
    if len(user_arg) == 1 and user_arg[0] == "send":
        SEND_TO_DROPBOX = True
    else:
        SEND_TO_DROPBOX = False
    
    main(SEND_TO_DROPBOX)


















