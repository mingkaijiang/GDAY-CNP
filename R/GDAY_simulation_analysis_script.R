##############################################################################
#### To analyze GDAY simulation output
#### Author: Mingkai Jiang
#### Created on: 25-May-2017
#### Modified on: 23-Jun-2017

########################## Main program ######################################
##############################################################################
#### This script analyze newer met results
#### i.e. amb_avg and ele_avg results

### Read in files
## DE relationship, N only
de_n_DF2 <- read.csv("outputs/EUC_amb_avg_02_ellsworth_N.csv", skip=1)
de_n_DF3 <- read.csv("outputs/EUC_ele_avg_03_ellsworth_N.csv", skip=1)

## DE relationship, N P 
de_np_DF2 <- read.csv("outputs/EUC_amb_avg_02_ellsworth_NP.csv", skip=1)
de_np_DF3 <- read.csv("outputs/EUC_ele_avg_03_ellsworth_NP.csv", skip=1)

## DE relationship, NP with TPU limitation considered
de_np_tp_DF2 <- read.csv("outputs/EUC_amb_avg_02_ellsworth_NP_tp.csv", skip=1)
de_np_tp_DF3 <- read.csv("outputs/EUC_ele_avg_03_ellsworth_NP_tp.csv", skip=1)

## Walker relationship, N only
wk_n_DF2 <- read.csv("outputs/EUC_amb_avg_02_walker_N.csv", skip=1)
wk_n_DF3 <- read.csv("outputs/EUC_ele_avg_03_walker_N.csv", skip=1)

## Walker relationship, N P
wk_np_DF2 <- read.csv("outputs/EUC_amb_avg_02_walker_NP.csv", skip=1)
wk_np_DF3 <- read.csv("outputs/EUC_ele_avg_03_walker_NP.csv", skip=1)


### Plot within year result
pdf("R/GDAY_simulation_within_year.pdf")
par(mfrow=c(2,1))
### Subset 1st year of data
de_n_DF2_sub <- de_n_DF2[1:366, ]
de_np_DF2_sub <- de_np_DF2[1:366, ]
de_np_tp_DF2_sub <- de_np_tp_DF2[1:366, ]
wk_n_DF2_sub <- wk_n_DF2[1:366, ]
wk_np_DF2_sub <- wk_np_DF2[1:366, ]

### Plotting 1st year GPP
with(de_n_DF2_sub, plot(gpp~doy, ylim=c(0, 0.15), type='l',
                        ylab="GPP [t/ha/yr]"))
with(de_np_DF2_sub, points(gpp~doy, type='l',col="red"))
with(de_np_tp_DF2_sub, points(gpp~doy, type='l',col="blue"))
with(wk_n_DF2_sub, points(gpp~doy, type='l',col="orange"))
with(wk_np_DF2_sub, points(gpp~doy, type='l',col="lightblue"))
legend("topright", c("DE_N", "DE_NP", "DE_NP_TP", "WK_N", "WK_NP"),
       col=c("black", "red", "blue", "orange", "lightblue"),
       lty = 1, lwd = 0.75, cex=0.5)
title("Year = 2012")

### Subset last year of data
s<- nrow(de_n_DF2)-364.0
e<- nrow(de_n_DF2)
de_n_DF2_sub <- de_n_DF2[s:e, ]
de_np_DF2_sub <- de_np_DF2[s:e, ]
de_np_tp_DF2_sub <- de_np_tp_DF2[s:e, ]
wk_n_DF2_sub <- wk_n_DF2[s:e, ]
wk_np_DF2_sub <- wk_np_DF2[s:e, ]

### Plotting 1st year GPP
with(de_n_DF2_sub, plot(gpp~doy, ylim=c(0, 0.15), type='l',
                        ylab="GPP [t/ha/yr]"))
with(de_np_DF2_sub, points(gpp~doy, type='l',col="red"))
with(de_np_tp_DF2_sub, points(gpp~doy, type='l',col="blue"))
with(wk_n_DF2_sub, points(gpp~doy, type='l',col="orange"))
with(wk_np_DF2_sub, points(gpp~doy, type='l',col="lightblue"))
legend("topright", c("DE_N", "DE_NP", "DE_NP_TP", "WK_N", "WK_NP"),
       col=c("black", "red", "blue", "orange", "lightblue"),
       lty = 1, lwd = 0.75, cex = 0.5)
title("Year = 2023")

dev.off()

#### Process the data to plot annual patterns
## Generate output DF at annual timestep
gppDF <- data.frame(seq(2012, 2023, by=1), NA, NA)
colnames(gppDF) <- c("Year", "gday_n","gday_p")
laiDF <- gppDF
leafDF <- gppDF
npDF <- gppDF

## store data
for (i in 2012:2023) {
    # GPP
    gppDF[gppDF$Year == i, "gday_n"] <- sum(de_n_DF2[de_n_DF2$year == i, "gpp"], na.rm=T)
    gppDF[gppDF$Year == i, "gday_p"] <- sum(de_np_DF2[de_np_DF2$year == i, "gpp"], na.rm=T)
    
    # LAI
    laiDF[laiDF$Year == i, "gday_n"] <- de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "lai"]
    laiDF[laiDF$Year == i, "gday_p"] <- de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "lai"]
    
    # Leaf C
    leafDF[leafDF$Year == i, "gday_n"] <- de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "shoot"]
    leafDF[leafDF$Year == i, "gday_p"] <- de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "shoot"]
    
    # Leaf N:P
    npDF[npDF$Year == i, "gday_n"] <- de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "shootn"]/de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "shootp"]
    npDF[npDF$Year == i, "gday_p"] <- de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "shootn"]/de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "shootp"]
    
}

# GPP
with(gppDF, plot(gday_n*100~Year, ylim=c(0,3000), 
                 ylab = expression(paste("GPP [g C ", m^-2, " ", yr^-1, "]")),
                 type="l", lwd = 2, col = "brown"))
with(gppDF, points(gday_p*100~Year, type="l", col = "red", lwd = 2))

# LAI
with(laiDF, plot(gday_n~Year, ylim=c(1,2), 
                 ylab = "LAI",
                 type="l", lwd = 2, col = "brown"))
with(laiDF, points(gday_p~Year, type="l", col = "red", lwd = 2))

# Leaf C
with(leafDF, plot(gday_n*100~Year, ylim=c(0,300), 
                 ylab = expression(paste("Leaf C [g C ", m^-2, " ", yr^-1, "]")),
                 type="l", lwd = 2, col = "brown"))
with(leafDF, points(gday_p*100~Year, type="l", col = "red", lwd = 2))

# Leaf NP
with(npDF, plot(gday_n~Year, ylim=c(0,100), 
                  ylab = "Leaf N:P ratio",
                  type="l", lwd = 2, col = "brown"))
with(npDF, points(gday_p~Year, type="l", col = "red", lwd = 2))
