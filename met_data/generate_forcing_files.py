#!/usr/bin/env python
"""
Create a G'DAY spinup & driving met forcing file for the EUCFACE experiment.
This is essentially to replicate the Medlyn et al. paper but adjusting the
met file script to match the latest incarnation of GDAY

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (22.08.2016)"
__email__ = "mdekauwe@gmail.com"

import sys
import glob
import os
import csv
import math
import numpy as np
from datetime import date
import calendar
import pandas as pd
import datetime as dt
import netCDF4 as nc


class CreateMetData(object):

    def __init__(self, site, soilT=True):

        self.site = site
                       
        self.ovar_names = ['#year', 'doy', 'tair', 'rain', 'tsoil',
                           'tam', 'tpm', 'tmin', 'tmax', 'tday', 'vpd_am',
                           'vpd_pm', 'co2', 'ndep', 'nfix', 'pdep', 'wind', 'pres',
                           'wind_am', 'wind_pm', 'par_am', 'par_pm']
                           
        self.ounits = ['#--', '--', 'degC', 'mm/d', 'degC','degC', 'degC',
                       'degC','degC', 'degC', 'kPa', 'kPa', 'ppm', 't/ha/d',
                       't/ha/d', 't/ha/d', 'm/s', 'kPa', 'm/s', 'm/s', 'mj/m2/d',
                       'mj/m2/d']

        self.soilT = soilT # if True model it using Tair & moving window

        # unit conversions
        self.SEC_TO_HFHR = 60.0 * 30.0 # sec/half hr
        self.PA_TO_KPA = 0.001
        self.J_TO_MJ = 1.0E-6
        self.K_to_C = 273.15
        self.G_M2_to_TONNES_HA = 0.01
        self.MM_S_TO_MM_HR = 3600.0
        self.SW_2_PAR = 2.3
        self.SEC_2_HLFHR = 1800.0
        self.J_TO_UMOL = 4.57
        self.UMOL_TO_J = 1.0 / self.J_TO_UMOL

    def write_spinup_file(self, df, yr_sequence, ofname=None, vary_co2=False,
                          co2=None, vary_ndep=False, ndep=None, nfix=None, 
                          vary_pdep=False, pdep=None):

        start_sim = yr_sequence[0]
        end_sim = yr_sequence[-1]
        year = str(start_sim)
        (ofp, wr) = self.write_hdr(yr_sequence, ofname, spinup=True)
        self.write_data(df, yr_sequence, ofp, wr, ofname, vary_co2, co2,
                        vary_ndep, ndep, nfix, vary_pdep, pdep)

    def write_daily_met_file(self, df, yr_sequence, ofname=None, vary_co2=False,
                             co2=None, vary_ndep=False, ndep=None, nfix=None,
                             vary_pdep=False, pdep=None):

        start_sim = yr_sequence[0]
        end_sim = yr_sequence[-1]
        year = str(start_sim)
        (ofp, wr) = self.write_hdr(yr_sequence, ofname, spinup=False)
        self.write_data(df, yr_sequence, ofp, wr, ofname, vary_co2, co2,
                        vary_ndep, ndep, nfix, vary_pdep, pdep)

    def write_hdr(self, yr_sequence, ofname, spinup=True):
        start_sim = yr_sequence[0]
        end_sim = yr_sequence[-1]
        year = str(start_sim)

        if spinup:
            tag = "spinup"
        else:
            tag = "forcing"

        try:
            ofp = open(ofname, 'w')
            wr = csv.writer(ofp, delimiter=',', quoting=csv.QUOTE_NONE,
                            escapechar=None, dialect='excel')
            wr.writerow(['# %s daily met %s' % (self.site, tag)])
            wr.writerow(['# Data from %s-%s' % (start_sim, end_sim)])
            wr.writerow(['# Created by Martin De Kauwe: %s' % date.today()])
            wr.writerow([var for i, var in enumerate(self.ounits)])
            wr.writerow([var for i, var in enumerate(self.ovar_names)])
        except IOError:
            raise IOError('Could not write met file: "%s"' % ofname)

        return (ofp, wr)

    def write_data(self, df, yr_sequence, ofp, wr, ofname, vary_co2, co2,
                   vary_ndep, ndep_data, nfix, vary_pdep, pdep_data):

        # use running mean - slow!!
        # or use day average
        if self.soilT:
            # There is no Tsoil data, so we are going to use the day average
            # of air temperature and a 7-day running mean to remove some of
            # the Tair=Tsoil
            tsoil = []
            dates = []
            for i, yr in enumerate(yr_sequence):
                days = np.unique(df[df.YEAR == yr].DOY)
                for j, doy in enumerate(days):
                    days_data = df[(df.YEAR == yr) & (df.DOY == doy)]
                    tsoil.append( np.mean(days_data["Tair"]-self.K_to_C) )

            # for the spinup stuff the years not being in order will mess the
            # dates up, so make up some random date series. It doesn't really
            # matter as we aren't using the actual date for anything other
            # than to interface with the pandas lib
            st = dt.datetime.strptime("01/01/80 00:00:00", '%d/%m/%y %H:%M:%S')
            nintervals = len(tsoil)
            dates = pd.date_range(st, periods=nintervals, freq='D')

            D = pd.Series(tsoil, dates)
            window_size = 7
            d_mva = D.rolling(window_size).mean()

            # The first few values will be nans, so we will use the 24-hr tair
            # values as replacements here
            for i in range(window_size-1):
                d_mva[i] = tsoil[i]
            tsoil_data = d_mva.values

        cnt = 0
        for i, yr in enumerate(yr_sequence):
            days = np.unique(df[df.YEAR == yr].DOY)
            for j, doy in enumerate(days):

                days_data = df[(df.YEAR == yr) & (df.DOY == doy)]

                morning = days_data[(days_data.HOUR <= 11.5) &
                                    (days_data.PAR >= 5.0)]
                afternoon = days_data[(days_data.HOUR >= 12.0) &
                                      (days_data.PAR >= 5.0)]
                day_light = days_data[days_data.PAR >= 5.0]

                # PAR data is crap, use the rest of the met data,
                # the PAR will be set below
                if len(morning) == 0 or len(afternoon):
                    morning = days_data[(days_data.HOUR >= 6.0) &
                                        (days_data.HOUR <= 11.5)]
                    afternoon = days_data[(days_data.HOUR >= 12.0) &
                                          (days_data.HOUR <= 19.0)]
                    day_light = days_data[(days_data.HOUR >= 6.0) &
                                          (days_data.HOUR <= 19.0)]


                #if yr == 2015 and doy == 31:
                #    print(np.mean(day_light["Tair"]-self.K_to_C))
                #    print(len(morning) , len(afternoon), len(day_light) )
                #    sys.exit()

                tair = np.mean(day_light["Tair"]-self.K_to_C)

                tam = np.mean(morning["Tair"]-self.K_to_C)
                tpm = np.mean(afternoon["Tair"]-self.K_to_C)
                if self.soilT:
                    tsoil = tsoil_data[cnt]
                else:
                    tsoil = np.mean(day_light["Tsoil"]-self.K_to_C)

                # daytime min/max temp
                tmin = np.min(days_data["Tair"]-self.K_to_C)
                tmax = np.max(days_data["Tair"]-self.K_to_C)
                tday = np.mean(days_data["Tair"]-self.K_to_C)


                if len(morning) == 0:
                    vpd_am = 0.05
                else:
                    vpd_am = np.mean(morning.VPD * self.PA_TO_KPA)
                    if vpd_am < 0.05:
                        vpd_am = 0.05

                if len(afternoon) == 0:
                    vpd_pm = 0.05
                else:
                    vpd_pm = np.mean(afternoon.VPD * self.PA_TO_KPA)
                    if vpd_pm < 0.05:
                        vpd_pm = 0.05

                # convert PAR [umol m-2 s-1] -> mj m-2 30min-1
                conv = self.UMOL_TO_J * self.J_TO_MJ * self.SEC_TO_HFHR
                par_am = np.sum(morning.PAR * conv)
                par_pm = np.sum(afternoon.PAR * conv)

                # rain -> mm
                # coversion 1800 seconds to half hours and summed gives day value.
                # Going to use the whole day including the night data
                rain = max(0.0, np.sum(days_data.Rainf * 1800.))

                # wind speed -> m/s
                wind = np.mean(day_light.Wind)
                if wind <= 0.0:
                    wind = 0.1 # set v.small speed but not zero

                # odd occasions when there is no data, so set a very small
                # wind speed.
                wind_am = np.mean(morning.Wind)
                if len(morning.Wind) == 0:
                    wind_am = 0.1

                wind_pm = np.mean(afternoon.Wind)
                if len(afternoon.Wind) == 0:
                    wind_pm = 0.1

                # air pressure -> kPa
                press = np.mean(day_light.PSurf * self.PA_TO_KPA)


                if vary_co2:
                    co2 = np.mean(days_data.CO2)

                if vary_ndep:
                    ndep = ndep_data[i]
                else:
                    ndep = ndep_data
                    
                if vary_pdep:
                    pdep = pdep_data[i]
                else:
                    pdep = 0.0

                wr.writerow([yr, doy, tair, rain, tsoil, tam, tpm, tmin, tmax,\
                             tday, vpd_am, vpd_pm, co2, ndep, nfix, pdep, wind, \
                             press, wind_am, wind_pm, par_am, par_pm])

                cnt += 1
        ofp.close()

    def get_random_year_sequence(self, start_yr, end_yr, out_yrs,
                                 preserve_leap=False):

        # Set the seed so we can repeat this if required
        np.random.seed(42)
        yrs = np.arange(start_yr, end_yr+1)

        if preserve_leap:

            # preserve leap yrs, so find them first
            leapyrs = np.zeros(0)
            for yr in yrs:
                if calendar.isleap(yr):
                    leapyrs = np.append(leapyrs, yr)


            # However we don't want the leapyrs in the sequence, so exclude them
            yrs = np.array([yrs[i] for i, yr in enumerate(yrs) \
                                    if yr not in leapyrs])

            shuff_years = self.randomise_array(out_yrs, yrs)
            shuff_years_leap = self.randomise_array(out_yrs, leapyrs)

            sequence = []
            i = 0
            for yr in np.arange(start_yr, end_yr+1):

                if i == 0:
                    prev_yr_leap = 1666 # anything not in the sequence
                    prev_yr = 1666 # anything not in the sequence

                if calendar.isleap(yr):
                    out_yr = shuff_years_leap[i]

                    # Make sure we don't get the same year twice
                    while prev_yr_leap == int(out_yr):
                        i += 1
                        out_yr = shuff_years_leap[i]

                    sequence.append(out_yr)
                    prev_yr_leap = shuff_years_leap[i]
                else:
                    out_yr = shuff_years[i]

                    # Make sure we don't get the same year twice
                    while prev_yr == int(out_yr):
                        i += 1
                        out_yr = shuff_years[i]

                    sequence.append(out_yr)
                    prev_yr = shuff_years[i]

                i += 1
        else:
            yrs = np.array([yrs[i] for i, yr in enumerate(yrs)])
            shuff_years = self.randomise_array(out_yrs, yrs)
            yr_list = np.random.uniform(start_yr, end_yr,
                                        out_yrs).astype(np.int)
            sequence = []
            i = 0
            for yr in yr_list:

                if i == 0:
                    prev_yr = 1666 # anything not in the sequence

                out_yr = shuff_years[i]

                # Make sure we don't get the same year twice
                while prev_yr == int(out_yr):
                    i += 1
                    out_yr = shuff_years[i]

                sequence.append(out_yr)
                prev_yr = shuff_years[i]

                i += 1

        return sequence

    def randomise_array(self, out_yrs, yrs):
        # make a sequence longer than the number of years we actually want
        num_seq = np.ones(out_yrs * len(yrs))
        num_years = len(yrs) * len(num_seq)
        shuff_years = (yrs * num_seq[:,None]).reshape(num_years)
        np.random.shuffle(shuff_years)

        return shuff_years

    def qair_to_vpd(self, qair, tair, press):

        # convert back to Pa
        press /= self.PA_TO_KPA

        # saturation vapor pressure
        es = 100.0 * 6.112 * np.exp((17.67 * tair) / (243.5 + tair))

        # vapor pressure
        ea = (qair * press) / (0.622 + (1.0 - 0.622) * qair)

        vpd = (es - ea) * self.PA_TO_KPA

        return vpd


def rename_colnames(df, frame_type="met"):

    if frame_type == "met":
        df = df.rename(columns={'SWdown [W/m^2]': 'SW',
                                'PAR [umol/m^2/s]': 'PAR',
                                'Tair [K]': 'Tair',
                                'Rainf [kg/m^2/s]': 'Rainf',
                                'Qair [kg/kg]': 'Qair',
                                'VPD [Pa]': 'VPD',
                                'Wind [m/s]': 'Wind',
                                'PSurf [Pa]': 'PSurf',
                                'CO2air [ppmv]': 'CO2',
                                'SoilTemp [K]': 'Tsoil'})
    elif frame_type == "ndep":
        df = df.rename(columns={'Year': 'YEAR',
                                'Doy': 'DOY',
                                'ambient CO2 [ppmv]': 'aCO2',
                                'elevated CO2 [ppmv]': 'eCO2',
                                'Ndep [kg N ha-1 yr-1]': 'ndep'})
    return df


if __name__ == "__main__":

    KG_2_TONNES = 1E-3
    G_M2_to_TONNES_HA = 0.01
    MM_2_CM = 0.1
    KG_HA_2_G_M2 = 0.1
    YR_TO_DAY = 365.25

    site = "EUC"
    met_dir = "raw_met_data"
    C = CreateMetData(site, soilT=False)

    exp_years = np.arange(1992, 2011)

    met_fname = "EucFACE_forcing_1992-2011.csv"
    df_met_varA = pd.read_csv(os.path.join(met_dir, met_fname), sep=",")
    df_met_varA = rename_colnames(df_met_varA, frame_type="met")


    # Following Cleveland et al. (1999) we calculate BNF as a function of ET
    # Similarly to Smith et al. (2014; LPJ-GUESS) and Wieder et al (2015; CLM)
    # we choose the conservative N fixation equation (Fig. 1).
    #
    # For ET we are using the sum of canopy evap and transpiration as using
    # total ET leads to very high BNF values in arid regions
    # (see Wieder et al. pg 3)
    #
    et = 300.0  #mm yr-1
    bn1 = 0.102
    bn2 = 0.524
    # bnf (kg N ha-1 yr-1) = (bn1 * (ET * mm_2_cm) + bn2) * kg_ha_2_g_m2
    nfix = (bn1 * (et * MM_2_CM) + bn2) * KG_HA_2_G_M2
    nfix *= G_M2_to_TONNES_HA / YR_TO_DAY
    
    pdep = 0.00000049

    # 2.25 kg N hectare-1 year-1 -> t/ha/day
    ndep = 0.00225 / YR_TO_DAY
    
    start_yr = exp_years[0]
    end_yr = exp_years[-1]
    num_yrs = 50
    yr_sequence = C.get_random_year_sequence(start_yr, end_yr, num_yrs,
                                             preserve_leap=False)
    
    #==========================
    # Generate Equilibrium FILE
    #==========================
    co2 = 400.0

    ofname = "%s_met_data_amb_var_co2.csv" % (site)
    C.write_daily_met_file(df_met_varA, yr_sequence, ofname=ofname, vary_co2=True,
                           co2=co2, vary_ndep=False, ndep=ndep, nfix=nfix, pdep=pdep)

    #==========================
    # Generate Elevated CO2 FILE
    #==========================
    co2 = 550.0

    ofname = "%s_met_data_ele_var_co2.csv" % (site)
    C.write_daily_met_file(df_met_varA, yr_sequence, ofname=ofname, vary_co2=False,
                           co2=co2, vary_ndep=False, ndep=ndep, nfix=nfix, pdep=pdep)
