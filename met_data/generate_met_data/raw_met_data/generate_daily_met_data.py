#!/usr/bin/env python
"""
The idea is this should be called by generate_FORCING_files.py, unless you just
want forcing for the experimental period.

Output relevant met drivers for GDAY model at the daily timestep using the
30 min raw data
"""

import calendar
import numpy as np
import datetime as dt
from datetime import date
import sys
import os
import solar_elevation as se
import csv

__author__  = "Martin De Kauwe"
__version__ = "1.0 (05.06.2013)"
__email__   = "mdekauwe@gmail.com"


def main(met_fname, ndep_fname, site, ofname, exp_years):

    # McCree, K. J. 1981. Photosynthetically active radiation. pp 41-55.
    # In: O. L. Lange, P. S. Nobel, C. B. Osmond and H. Ziegler (Ed).
    # Encyclopedia of Plant Physiology. New Series Vol 12-A. Springer-Verlag.
    # Berlin.
    # umol/m2/sec = 4.6 umol/J
    # Solar radiation to PAR ~ 0.5 (0.5 MJ PAR / MJ SRAD)
    # SRAD in Mj/m2 d to par in mol/m2 s = 4.6 * 0.5 ~ 2

    UMOL_TO_J = 1.0 / 4.6 # coversion from umol quanta to J, 118.708 / 550nm
    KG_2_TONNES = 1E-3
    YR_TO_DAY = 365.25
    SEC_TO_HFHR = 60.0 * 30.0 # sec/half hr
    PA_TO_KPA = 0.001
    J_TO_MJ = 1.0E-6
    PAR_TO_SW = 2.3 # assuming a conversion factor of 2.3 umol/m^2/s per W/m^2
    K_to_C = 273.15
    PASCALS_TO_KPA = 1E-03

    names = ["year","doy","hour","sw","par","lw","tair","rf","sf","qair",\
             "vpd","rh","wind","psurf","co2","tsoil"]
    var_dict = dict(zip(names, np.arange(16)))

    data = np.loadtxt(met_fname, delimiter=",", skiprows=1)

    # read varying ndep data original units kg N ha-1 yr-1
    ndep_data = np.loadtxt(ndep_fname, delimiter=" ", skiprows=1)
    ndep = ndep_data[ndep_data[:,0]>=2012][:,4] * KG_2_TONNES / YR_TO_DAY

    ovar_names = ['#year', 'doy', 'sw_rad', 'tair', 'rain', 'tsoil',
                  'tam', 'tpm', 'vpd_am', 'vpd_pm', 'vpd_avg', 'co2',
                  'ndep', 'wind', 'atmos_press', 'par', 'wind_am',
                  'wind_pm', 'sw_rad_am', 'sw_rad_pm']
    ounits = ['#--', '--', 'mj/m2/day', 'c', 'mm', 'c', 'c',
              'c', 'kPa', 'kPa', 'kPa', 'ppm', 't/ha/year', 'm/s',
              'kPa', 'umol/m2/d', 'm/s', 'm/s', 'mj/m2/am', 'mj/m2/pm']

    start_sim = 2012
    end_sim = 2023
    try:
        ofp = open(ofname, 'wb')
        wr = csv.writer(ofp, delimiter=',', quoting=csv.QUOTE_NONE,
                        escapechar=None, dialect='excel')
        wr.writerow(['# %s daily met forcing' % (site)])
        wr.writerow(['# Data from %s-%s' % (start_sim, end_sim)])
        wr.writerow(['# Created by Martin De Kauwe: %s' % date.today()])
        #wr.writerow(["#" if i==0 else var for i, var in enumerate(ounits)])
        wr.writerow([var for i, var in enumerate(ounits)])
        wr.writerow([var for i, var in enumerate(ovar_names)])

    except IOError:
        raise IOError('Could not write met file: "%s"' % ofname)

    non_lp_yr_length = 365 * 48
    numvars = len(names)
    numvars2 = 11 # co2
    prjday = 0
    daily_day_count = 0 # used for data read from daily file
    for yr in exp_years:
        yr_data = data[np.where(data[:,var_dict["year"]] == yr)]

        if yr_data.shape[0] == non_lp_yr_length:
            ndays = 365
            yr_data = yr_data.reshape(ndays, 48, numvars)

        else:
            ndays = 366
            yr_data = yr_data.reshape(ndays, 48, numvars)

        for day in xrange(1, ndays+1):

            days_data = yr_data[np.where(yr_data[:,:,var_dict["doy"]] == day)]

            lat = -33.59999847
            lon = 150.75

            sun_up = np.array([se.calc_solar_elev(lat, float(hr), day)\
                               for hr in xrange(24)])
            sun_up_index = np.array([np.where(sun_up > 0.0)])

            morn_index = sun_up_index[np.where(sun_up_index<=12.0)]
            aftern_index = sun_up_index[np.where(sun_up_index>12.0)]

            morning = days_data[np.where(np.logical_and(
                            days_data[:,var_dict["hour"]] >= morn_index.min(),
                            days_data[:,var_dict["hour"]] <= morn_index.max()))]
            afternoon = days_data[np.where(np.logical_and(
                            days_data[:,var_dict["hour"]] >= aftern_index.min(),
                            days_data[:,var_dict["hour"]] <= aftern_index.max()))]


            morning = morning[np.where(morning[:,var_dict["par"]] >= 0.0)]
            afternoon = afternoon[np.where(afternoon[:,var_dict["par"]] >= 0.0)]


            # temp -> degC
            tmean = np.mean(np.hstack((morning[:,var_dict["tair"]]-K_to_C, \
                                       afternoon[:,var_dict["tair"]]-K_to_C)))
            tam = np.mean(morning[:,var_dict["tair"]]-K_to_C)
            tpm = np.mean(afternoon[:,var_dict["tair"]]-K_to_C)


            # vpd -> kPa
            if len(morning) == 0:
                vpd_am = 0.05
            else:
                vpd_am = np.mean(morning[:,var_dict["vpd"]] * PA_TO_KPA)
                if vpd_am < 0.05:
                    vpd_am = 0.05

            if len(afternoon) == 0:
                vpd_pm = 0.05
            else:
                vpd_pm = np.mean(afternoon[:,var_dict["vpd"]] * PA_TO_KPA)
                if vpd_pm < 0.05:
                    vpd_pm = 0.05


            vpd_avg = np.mean(np.hstack((morning[:,var_dict["vpd"]]* PA_TO_KPA, \
                              afternoon[:,var_dict["vpd"]]* PA_TO_KPA)))
            if vpd_avg < 0.05:
                vpd_avg = 0.05

            # SW_down [W/m2]=[J m-2 s-1] -> MJ m-2 d-1
            conv = J_TO_MJ * SEC_TO_HFHR
            swdown_morning = morning[:,var_dict["sw"]] * conv
            swdown_afternoon = afternoon[:,var_dict["sw"]] * conv
            sw_rad = np.sum(np.hstack((swdown_morning, swdown_afternoon)))

            # convert PAR [umol m-2 s-1] -> umol m-2 d-1
            conv = SEC_TO_HFHR
            par_morning = morning[:,var_dict["par"]] * conv
            par_afternoon = afternoon[:,var_dict["par"]] * conv
            par_day = np.sum(np.hstack((par_morning, par_afternoon)))

            tsoil = np.mean(np.hstack((morning[:,var_dict["tsoil"]]-K_to_C, \
                                       afternoon[:,var_dict["tsoil"]]-K_to_C)))

            # rain -> mm
            # coversion 1800 seconds to half hours and summed gives day value.
            #rain = np.sum(np.hstack((morning[:,var_dict["rf"]]*1800., \
            #                           afternoon[:,var_dict["rf"]]*1800.)))
            # Going to use the whole day including the night data
            rain = np.sum(days_data[:,var_dict["rf"]]*1800.)

            # wind speed -> m/s
            wind_sp = np.mean(np.hstack((morning[:,var_dict["wind"]], \
                                       afternoon[:,var_dict["wind"]])))

            # odd occasions whne there is no data, so set a very small wind speed.
            wind_am = np.mean(morning[:,var_dict["wind"]])
            if len(morning[:,var_dict["wind"]]) == 0:
                wind_am = 0.1
            wind_pm = np.mean(afternoon[:,var_dict["wind"]])
            if len(afternoon[:,var_dict["wind"]]) == 0:
                wind_pm = 0.1

            # air pressure -> kPa
            atpress = np.mean(np.hstack((morning[:,var_dict["psurf"]]* PA_TO_KPA, \
                                       afternoon[:,var_dict["psurf"]]* PA_TO_KPA)))

            # co2 -> [ppm]
            # Need to take the whole day mean for CO2 to get the correct
            # forcing data, not the daylight only data.
            co2 = np.mean(days_data[:,var_dict["co2"]])


            wr.writerow([yr, day, sw_rad, tmean, rain, tsoil, tam, tpm, vpd_am, \
                         vpd_pm, vpd_avg, co2, ndep[prjday], \
                         wind_sp, atpress, par_day, wind_am, wind_pm,\
                         np.sum(swdown_morning), np.sum(swdown_afternoon)])

            prjday += 1
        daily_day_count+=1
    ofp.close()




if __name__ == "__main__":

    met_fname = "EucFACE_forcing_2012-2023_AMBVAR.csv"
    ndep_fname = "EucFACE_forcing_daily_CO2NDEP_1750-2023.dat"
    ofname = "euc_run_amb_var_2012-2023.csv"
    exp_years = np.arange(2012, 2024)
    main(met_fname, ndep_fname, ofname, exp_years)

    print ofname
