##############################################################################
#### To analyze GDAY simulation output
#### Author: Mingkai Jiang
#### Created on: 25-May-2017

########################## Main program ######################################
##############################################################################

#### Quality control - check gday output continuity

### Read in files
ambDF <- read.csv("outputs/EUC_amb_equilib.csv", skip=1)
eleDF1 <- read.csv("outputs/EUC_ele_initial.csv", skip=1)
eleDF2 <- read.csv("outputs/EUC_ele_final_equilib.csv", skip=1)

### merge
mgDF <- rbind(ambDF, eleDF1)


unique(ambDF$doy)
