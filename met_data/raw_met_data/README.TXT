!--------------------------------------------------------------------------------!
! This is the readme file for                                                    !
! meteorolgical forcing files for the EucFACE MIP activity                       !
! created by S. Zaehle (szaehle@bgc-jena.mpg.de) on 04/06/2013                   !
! -- revised 16/07/2013
!--------------------------------------------------------------------------------!

File naming convention:
*.csv -> comma separated ascii files
*_ALMA.nc -> ALMA convention like netcdf files
*_CDF.nc -> CDF compliant netcdf file

The period 1992-2011 is the meteorological data for the spin-ip
The period 2012-2023 is the meteorological data for the actual simulation 
- 4 version exist for the simulation period 
  AMBAVG: annual CO2 increase as in the RCP3PD scenario
          meteorological data from 1998 (corrected for leap-years where necessary)
  AMBVAR: annual CO2 increase as in the RCP3PD scenario
          meteorological data from 1998-2010 
  ELEAVG: annual CO2 increase as in the RCP3PD scenario + 150ppm from Sept 2012 
          meteorological data from 1998 (corrected for leap-years where necessary)
  ELEVAR: annual CO2 increase as in the RCP3PD scenario + 150ppm from Sept 2012
          meteorological data from 1998-2010 

The files *CO2NDEP_1860-2023.dat contain annual/daily CO2 and N deposition 
  values for spin-up and model simulation

The data were generate using
- hourly data from Princton forcing data set, disaggregated in time using
  the weather generator of the CABLE model. Comparison of daily data with local 
  site meteorology suggests an acceptable fit.
- annual CO2 are rom the RCP historic period and future RCP3PD scenario (2005-2023)
- eCO2 treatment follows the experimental protocol of the site (pers.comm.)
- N dep follows the relative increase in atm. CO2 and reaches 3 kg N / ha / yr
  in the year 2000. 

Modifications to original data:
- Air pressure has been set to 1013 hPa
- PAR has been generated from incident downward shortwave radiation, 
  assuming a conversion factor of 2.3 umol/m^2/s per W/m^2
- VPD and RH have been calculated from Qair using the approach of 
  Montheith & Unsworth 2008
- Wind data have been taken from Bernard Pak, derived with Cable from the 
  Sheffield climate data 1992-2008. Values for 2009-2011 correspond to those
  of 1992-1994
- solar noon was corrected to be at 12 hours!
!--------------------------------------------------------------------------------!

