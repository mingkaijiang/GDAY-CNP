##############################################################################
#### To analyze GDAY simulation output
#### Author: Mingkai Jiang
#### Created on: 25-May-2017
#### Modified on: 23-Jun-2017

########################## Main program ######################################
##############################################################################


### Read in files
## DE relationship, N only
de_n_DF1 <- read.csv("outputs/EUC_amb_01_ellsworth_N.csv", skip=1)
de_n_DF2 <- read.csv("outputs/EUC_ele_02_ellsworth_N.csv", skip=1)
de_n_DF3 <- read.csv("outputs/EUC_ele_03_ellsworth_N.csv", skip=1)

## DE relationship, N P 
de_np_DF1 <- read.csv("outputs/EUC_amb_01_ellsworth_NP.csv", skip=1)
de_np_DF2 <- read.csv("outputs/EUC_ele_02_ellsworth_NP.csv", skip=1)
de_np_DF3 <- read.csv("outputs/EUC_ele_03_ellsworth_NP.csv", skip=1)

## DE relationship, NP with TPU limitation considered
de_np_tp_DF1 <- read.csv("outputs/EUC_amb_01_ellsworth_NP_tp.csv", skip=1)
de_np_tp_DF2 <- read.csv("outputs/EUC_ele_02_ellsworth_NP_tp.csv", skip=1)
de_np_tp_DF3 <- read.csv("outputs/EUC_ele_03_ellsworth_NP_tp.csv", skip=1)

## Walker relationship, N only
wk_n_DF1 <- read.csv("outputs/EUC_amb_01_walker_N.csv", skip=1)
wk_n_DF2 <- read.csv("outputs/EUC_ele_02_walker_N.csv", skip=1)
wk_n_DF3 <- read.csv("outputs/EUC_ele_03_walker_N.csv", skip=1)

## Walker relationship, N P
wk_np_DF1 <- read.csv("outputs/EUC_amb_01_walker_NP.csv", skip=1)
wk_np_DF2 <- read.csv("outputs/EUC_ele_02_walker_NP.csv", skip=1)
wk_np_DF3 <- read.csv("outputs/EUC_ele_03_walker_NP.csv", skip=1)


### Subset 1st year of data
de_n_DF1_sub <- de_n_DF1[1:366, ]
de_np_DF1_sub <- de_np_DF1[1:366, ]
de_np_tp_DF1_sub <- de_np_tp_DF1[1:366, ]
wk_n_DF1_sub <- wk_n_DF1[1:366, ]
wk_np_DF1_sub <- wk_np_DF1[1:366, ]

### Plotting 1st year GPP
with(de_n_DF1_sub, plot(gpp~doy, ylim=c(0, 0.15), type='l',
                        ylab="GPP [t/ha/yr]"))
with(de_np_DF1_sub, points(gpp~doy, type='l',col="red"))
with(de_np_tp_DF1_sub, points(gpp~doy, type='l',col="blue"))
with(wk_n_DF1_sub, points(gpp~doy, type='l',col="orange"))
with(wk_np_DF1_sub, points(gpp~doy, type='l',col="lightblue"))

### Subset last year of data
s<- nrow(de_n_DF1)-364.0
e<- nrow(de_n_DF1)
de_n_DF1_sub <- de_n_DF1[s:e, ]
de_np_DF1_sub <- de_np_DF1[s:e, ]
de_np_tp_DF1_sub <- de_np_tp_DF1[s:e, ]
wk_n_DF1_sub <- wk_n_DF1[s:e, ]
wk_np_DF1_sub <- wk_np_DF1[s:e, ]

### Plotting 1st year GPP
with(de_n_DF1_sub, plot(gpp~doy, ylim=c(0, 0.15), type='l',
                        ylab="GPP [t/ha/yr]"))
with(de_np_DF1_sub, points(gpp~doy, type='l',col="red"))
with(de_np_tp_DF1_sub, points(gpp~doy, type='l',col="blue"))
with(wk_n_DF1_sub, points(gpp~doy, type='l',col="orange"))
with(wk_np_DF1_sub, points(gpp~doy, type='l',col="lightblue"))


## We don't have P limitation over time. Need to re-parameterize the model - update CO2 result, fixed using 1998 met data
## Update met data so that it reflects Medlyn's paper period (2012-2023) - run using 1998 data both ambCO2 and eleCO2
## Plot CO2 effect at annnual timestep  - add new P limitation GDAY result

## To do list:
# 1. Update met data (fixed 1998 met, variable CO2 for 2012-2023, historic post-industrial transient data)
# 2. Parameterize so that the site is P-limited
# 3. Add plot onto Medlyn 2016
# 4. Generate results for DE
# 5. Make poster using quasi-equil framework (i.e. re-parameterize the simplfied GDAY)
# 6. AmazonFACE talk 
