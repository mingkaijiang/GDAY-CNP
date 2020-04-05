#!/usr/bin/env python

"""
Create GDAY met forcing files for the EucFACE experiment.
Note we are outputting "daytime" averages not whole day averages...

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (25.01.2011)"
__email__ = "mdekauwe@gmail.com"

import sys
import os
import csv
import math
import numpy as np
from datetime import date
import datetime
import calendar
import csv

import sys
sys.path.append(os.getcwd())
import generate_daily_met_data as mk_daily_forcing

def main(site, daily_data, yr_sequence, ofname=None, out_yrs=None,
         vary_co2=False, co2=None, vary_ndep=False, ndep=None):
    """ Generate an equilibrium input forcing file, something like 4000 yrs
    we are using random repeated years...
    """
    start_sim = int(daily_data[:,0][0])
    end_sim = int(daily_data[:,0][-1])
    year = str(start_sim)
    month = str(01)
    day = str(01)
    date = datetime.datetime.strptime((year + month + day), "%Y%m%d")

    data_yr, data_doy = [], []
    for i in xrange(len(daily_data)):
        data_yr.append(date.year)
        data_doy.append(int(date.strftime('%j')))
        date += datetime.timedelta(days=1)

    ovar_names = ['#year', 'doy', 'sw_rad', 'tair', 'rain', 'tsoil',
                  'tam', 'tpm', 'vpd_am', 'vpd_pm', 'vpd_avg', 'co2',
                  'ndep', 'wind', 'atmos_press', 'par', 'wind_am',
                  'wind_pm', 'sw_rad_am', 'sw_rad_pm']
    ounits = ['#--', '--', 'mj/m2/day', 'c', 'mm', 'c', 'c',
              'c', 'kPa', 'kPa', 'kPa', 'ppm', 't/ha/year', 'm/s',
              'kPa', 'umol/m2/d', 'm/s', 'm/s', 'mj/m2/am', 'mj/m2/pm']


    try:
        ofp = open(ofname, 'wb')
        wr = csv.writer(ofp, delimiter=',', quoting=csv.QUOTE_NONE,
                        escapechar=None, dialect='excel')
        wr.writerow(['# %s daily met forcing' % (site)])
        wr.writerow(['# Data from %s-%s' % (start_sim, end_sim)])
        wr.writerow(['# Created by Martin De Kauwe: %s' % date.today()])
        wr.writerow([var for i, var in enumerate(ounits)])
        wr.writerow([var for i, var in enumerate(ovar_names)])

    except IOError:
        raise IOError('Could not write met file: "%s"' % ofname)

    yrs = np.arange(start_sim, end_sim + 1)

    # make a sequence longer than the number of years we actually want
    # as we may have supplied an odd number, there is a neater way to do this
    # see grassFACE runs, but I didn't use that here!
    num_seq = np.ones(out_yrs * len(yrs))
    num_years = len(yrs) * len(num_seq)

    # shuffle the years so we get "random" weather for our equilibrium
    # simulations
    shuff_years = (yrs * num_seq[:,None]).reshape(num_years)
    np.random.shuffle(shuff_years)

    # now just take the number of years the user requested, we can just
    # take the first lot as these should be random anyway
    shuff_years = shuff_years[:out_yrs]
    year_doy = np.array([ data_yr, data_doy ])

    for yr_index, y in enumerate(shuff_years):

        # get the index of a given shuffled year
        yrs_index = np.asarray(np.where(year_doy[0,:] == y))

        for i in yrs_index[0,:]:
            wr.writerow([yr_sequence[yr_index], \
                         daily_data[i,1], daily_data[i,2], \
                         daily_data[i,3], daily_data[i,4], \
                         daily_data[i,5], daily_data[i,6], \
                         daily_data[i,7], daily_data[i,8], \
                         daily_data[i,9], daily_data[i,10],\
                         co2[yr_index], ndep[yr_index],\
                         daily_data[i,13], daily_data[i,14],\
                         daily_data[i,15], daily_data[i,16],\
                         daily_data[i,17],daily_data[i,18],\
                         daily_data[i,19]])
    ofp.close()

if __name__ == "__main__":

    KG_2_TONNES = 1E-3
    YR_TO_DAY = 365.25
    site = "EUC"
    met_dir = "raw_met_data/model_forcing"

    #=============================================
    # make the experimental run daily driving file -AMB/ELE, var/avg
    #=============================================
    met_fname = os.path.join(met_dir, "EucFACE_forcing_2012-2023_AMBAVG.csv")
    ndep_fname = os.path.join(met_dir, "EucFACE_forcing_daily_CO2NDEP_1750-2023.dat")
    ofname = "%s_met_data_amb_avg_co2.csv" % (site)
    exp_years = np.arange(2012, 2024)
    mk_daily_forcing.main(met_fname, ndep_fname, site, ofname, exp_years)

    met_fname = os.path.join(met_dir, "EucFACE_forcing_2012-2023_AMBVAR.csv")
    ndep_fname = os.path.join(met_dir, "EucFACE_forcing_daily_CO2NDEP_1750-2023.dat")
    ofname = "%s_met_data_amb_var_co2.csv" % (site)
    exp_years = np.arange(2012, 2024)
    mk_daily_forcing.main(met_fname, ndep_fname, site, ofname, exp_years)

    met_fname = os.path.join(met_dir, "EucFACE_forcing_2012-2023_ELEAVG.csv")
    ndep_fname = os.path.join(met_dir, "EucFACE_forcing_daily_CO2NDEP_1750-2023.dat")
    ofname = "%s_met_data_ele_avg_co2.csv" % (site)
    exp_years = np.arange(2012, 2024)
    mk_daily_forcing.main(met_fname, ndep_fname, site, ofname, exp_years)

    met_fname = os.path.join(met_dir, "EucFACE_forcing_2012-2023_ELEVAR.csv")
    ndep_fname = os.path.join(met_dir, "EucFACE_forcing_daily_CO2NDEP_1750-2023.dat")
    ofname = "%s_met_data_ele_var_co2.csv" % (site)
    exp_years = np.arange(2012, 2024)
    mk_daily_forcing.main(met_fname, ndep_fname, site, ofname, exp_years)

    #=============================================
    # make the experimental run daily driving file -AMB
    #=============================================

    # read daily file back in...can just use the ambient one.
    ofname = "%s_met_data_amb_var_co2.csv" % (site)
    daily_data = np.loadtxt(ofname, skiprows=5, delimiter=",")


    #==========================
    # Generate Equilibrium FILE
    #==========================
    ofname = "%s_met_data_equilibrium_50_yrs.csv" % (site)
    # Make a forcing file to cover a run to equilibrium
    num_yrs = 50
    years = np.unique(daily_data[:,0])
    yr_sequence = np.arange(years[0],years[0]+num_yrs)
    co2 = [276.84] * num_yrs
    #2.25 kg N hectare-1 year-1
    ndep = [0.00225 / YR_TO_DAY] * num_yrs #  t/ha/yr --> t/ha/day

    main(site, daily_data, yr_sequence, ofname=ofname, out_yrs=num_yrs,
         vary_co2=False, co2=co2, vary_ndep=False, ndep=ndep)


    #================================================================
    # Make a forcing file from 1750 years with increasing co2/ndep
    #================================================================
    end_yr = int(years[0])-1
    num_yrs = (end_yr-1750)+1 # 1750-year before experiment
    yr_sequence = np.arange(years[0],years[0]+num_yrs)



    # read global co2 for equilibrium values
    fname = "EucFACE_forcing_daily_CO2NDEP_1750-2023.dat"
    fname = os.path.join(met_dir, fname)
    co2_ndep_data = np.loadtxt(fname, delimiter=" ", skiprows=1)
    # Need the additional index for DOY =1, otherwise we want get a different
    # CO2 per year.
    global_co2_data = co2_ndep_data[(co2_ndep_data[:,0]>=1750.0) & \
                                    (co2_ndep_data[:,1]==1) & \
                                    (co2_ndep_data[:,0]<=end_yr)][:,2]

    # read varying ndep data
    # Need the additional index for DOY =1, otherwise we want get a different
    # NDEP per year.
    ndep_data = co2_ndep_data[(co2_ndep_data[:,0]>=1750.0) & \
                              (co2_ndep_data[:,1]==1) & \
                              (co2_ndep_data[:,0]<=end_yr)][:,4]
    ndep_data = ndep_data * KG_2_TONNES / YR_TO_DAY


    ofname = "%s_met_data_industrial_to_present_1750_%d.csv" % (site, end_yr)
    main(site, daily_data, yr_sequence, ofname=ofname, out_yrs=num_yrs,
         vary_co2=True, co2=global_co2_data, vary_ndep=True, ndep=ndep_data)
