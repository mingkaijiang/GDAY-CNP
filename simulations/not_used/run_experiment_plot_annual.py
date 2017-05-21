#!/usr/bin/env python

""" run GDAY, plot LAI, GPP, shoot NC """

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import plotsettings as ps
import pandas as pd
import datetime as dt

__author__  = "Martin De Kauwe"
__version__ = "1.0 (16.05.2013)"
__email__   = "mdekauwe@gmail.com"

def main(site, experiment):

    # run new simulations
    os.system("eucface_simulations.py")

    # load data
    amb = read_data("../outputs/D1GDAY%sAMB%s.csv" % (site, experiment))
    ele = read_data("../outputs/D1GDAY%sELE%s.csv" % (site, experiment))

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

    ax1.plot(years, amb.groupby("YEAR").NPP.sum(), "b-", label="Amb")
    ax1.plot(years, ele.groupby("YEAR").NPP.sum(), "r-", label="Ele")

    #ax1.legend(numpoints=1, loc="best")
    ax1.set_ylabel("NPP")
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylim(0., 1200)
    ax1.set_xlim(2012, 2024)

    ax2 = fig.add_subplot(312)
    ax2.plot(years, amb.groupby("YEAR").LAI.max(), "b^", label="Amb")
    ax2.plot(years, ele.groupby("YEAR").LAI.max(), "r^", label="Ele")

    #ax2.legend(numpoints=1, loc="best")
    ax2.set_ylabel("LAI")
    ax2.set_xlabel("Year")
    ax2.set_ylim(0., 3)
    ax2.set_xlim(2012, 2024)

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

    gpp_response = ((gpp_ele / gpp_amb) - 1.0) * 100.0
    trans_response = ((trans_ele / trans_amb) - 1.0) * 100.0
    wue_response = ((wue_ele / wue_amb) - 1.0) * 100.0


    ax3 = fig.add_subplot(313)
    ax3.plot(years, gpp_response, "r-", label="GPP")
    ax3.plot(years, trans_response, "g-", label="T")
    ax3.plot(years, wue_response, "b-",label="WUE")

    ax3.legend(numpoints=1, loc="best")
    ax3.set_ylabel("CO$_2$ response")
    ax3.set_xlabel("Year")
    ax3.set_ylim(-30., 50)
    ax3.set_xlim(2012, 2024)
    plt.show()

    print("GPP: %.2f +/- %.2f" % (gpp_response.mean(), gpp_response.std()))
    print("T: %.2f +/- %.2f" % (trans_response.mean(), trans_response.std()))
    print("WUE: %.2f +/- %.2f" % (wue_response.mean(), wue_response.std()))
    print("\n")
    print("leaf NC amb: %.2f; leaf NC ele: %.2f" % (np.mean(shootn_amb/shoot_amb), np.mean(shootn_ele/shoot_ele)))

def date_converter(*args):
    return dt.datetime.strptime(str(int(float(args[0]))) + " " +\
                                str(int(float(args[1]))), '%Y %j')

def read_data(fname):
    df = pd.read_csv(fname, parse_dates=[[0,1]], index_col=0, sep=",",
                     keep_date_col=True, date_parser=date_converter,
                     na_values=["-999.9"], skiprows=3)

    #df = pd.read_csv(fname, header=0, sep=",")

    return df




if __name__ == "__main__":

    site = "EUC"
    experiment = "VAR"
    main(site, experiment)
