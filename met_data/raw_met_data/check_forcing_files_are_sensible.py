#!/usr/bin/env python

"""
Make some integrated plots to make sure the data looks OK.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (04.06.2013)"
__email__ = "mdekauwe@gmail.com"

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import sys
import matplotlib.mlab as mlab


def main(fname, yr=1990, end_year=2011):
   
    names = ["year","doy","hour","sw","par","lw","tair","rf","sf","qair",\
             "vpd","rh","wind","psurf","co2","tsoil"]
    var_dict = dict(zip(names, np.arange(16)))
    
    names="year,doy,sw,par,lw,tair,rf,sf,qair,vpd,rh,wind,psurf,co2,tsoil"
    data = np.loadtxt(fname, delimiter=",", skiprows=1)

    # Make daily data (not DAYLIGHT!!)
    K_to_C = 272.15
    PASCALS_TO_KPA = 1E-03
    
    of = open(fname.split(".")[0] + "_daily.csv", "w")    
    print >> of, names
    
    for i in xrange(len(data)):
        if yr >= end_year:
            break
       
        yrs_data = data[np.where(data[:,0]==yr)]
        ndays = int(len(yrs_data) / 48.)
        
        for d in xrange(ndays):
            days_data = yrs_data[np.where(yrs_data[:,1]==d+1)]
            
            print >> of, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f" % (yr, d+1, \
                                                                               days_data[:,var_dict["sw"]].mean(),   \
                                                                               days_data[:,var_dict["par"]].mean(),   \
                                                                               days_data[:,var_dict["lw"]].mean(),   \
                                                                               days_data[:,var_dict["tair"]].mean()-K_to_C,\
                                                                               np.sum(days_data[:,var_dict["rf"]]*1800.),   \
                                                                               days_data[:,var_dict["sf"]].sum(),   \
                                                                               days_data[:,var_dict["qair"]].mean(), \
                                                                               np.mean(days_data[:,var_dict["vpd"]]*PASCALS_TO_KPA), \
                                                                               days_data[:,var_dict["rh"]].mean(), \
                                                                               days_data[:,var_dict["wind"]].mean(), \
                                                                               days_data[:,var_dict["psurf"]].mean(),\
                                                                               days_data[:,var_dict["co2"]].mean(),  \
                                                                               days_data[:,var_dict["tsoil"]].mean()-K_to_C)
        yr += 1
    of.close()
                                                                                    
def make_plot(data, treat, VAR_TYPE, yr=1990, end_year=2011):
    
    names = ["year","doy","sw","par","lw","tair","rf","sf","qair",\
             "vpd","rh","wind","psurf","co2","tsoil"]
    var_dict = dict(zip(names, np.arange(15)))
    
    for key,var_index in var_dict.iteritems():
        
        fig = plt.figure(figsize=(10,10))
        fig.subplots_adjust(hspace=0.9)
        fig.subplots_adjust(wspace=0.9)
    
        for i in np.arange(yr, end_year):
            yrs_data = data[np.where(data[:,0]==i)]
    
            ax = fig.add_subplot(3,5,(i-yr)+1)
            ax.plot(yrs_data[:,1], yrs_data[:,var_index], "r-")
            ax.set_title(str(i))
            ax.set_xticks(np.array([0, 182, 365]))
            xlabels = ax.get_xticklabels()
            # rotate ticks
            for label in xlabels:
                label.set_rotation(30)
        plt.savefig("plots/" + treat + "_" + VAR_TYPE + "_" + key + ".png", dpi=150)

    
def plot_co2_ndep(data):
    
    fig = plt.figure(figsize=(10,10))    
    ax = fig.add_subplot(211)
    ax.plot(data[:,2], "b-", label="Amb")
    ax.plot(data[:,3], "r-", label="Ele")
    ax.set_title("CO$_2$")
    ax.legend(numpoints=1, loc="best")
    
    ax = fig.add_subplot(212)
    ax.plot(data[:,4], "b-", label="Amb")
    ax.set_title("Ndep (kg N ha$^{-1}$ yr$^{-1}$)")
    plt.savefig("plots/daily_CO2_ndep.png")

if __name__ == "__main__":
    
    
    fname = "EucFACE_forcing_2012-2023_AMBVAR.csv"
    make_daily = True
    if make_daily:
        main(fname, yr=2012, end_year=2023)
    
    # Read back in, make some plots
    fname = fname.split(".")[0] + "_daily.csv"
    data = np.loadtxt(fname, delimiter=",", skiprows=1)
    
    make_plot(data, "AMB", "VAR", yr=2012, end_year=2023)
    
    fname = "EucFACE_forcing_2012-2023_AMBAVG.csv"
    make_daily = True
    if make_daily:
        main(fname, yr=2012, end_year=2023)
    
    # Read back in, make some plots
    fname = fname.split(".")[0] + "_daily.csv"
    data = np.loadtxt(fname, delimiter=",", skiprows=1)
    
    make_plot(data, "AMB", "AVG", yr=2012, end_year=2023)
    
    
    fname = "EucFACE_forcing_2012-2023_ELEVAR.csv"
    make_daily = True
    if make_daily:
        main(fname, yr=2012, end_year=2023)
    
    # Read back in, make some plots
    fname = fname.split(".")[0] + "_daily.csv"
    data = np.loadtxt(fname, delimiter=",", skiprows=1)
    
    make_plot(data, "ELE", "VAR", yr=2012, end_year=2023)
    
    fname = "EucFACE_forcing_2012-2023_ELEAVG.csv"
    make_daily = True
    if make_daily:
        main(fname, yr=2012, end_year=2023)
    
    # Read back in, make some plots
    fname = fname.split(".")[0] + "_daily.csv"
    data = np.loadtxt(fname, delimiter=",", skiprows=1)
    
    make_plot(data, "ELE", "AVG", yr=2012, end_year=2023)
    
    
    fname = "EucFACE_forcing_daily_CO2NDEP_1750-2023.dat"
    data = np.loadtxt(fname, delimiter=" ", skiprows=1)
    plot_co2_ndep(data)