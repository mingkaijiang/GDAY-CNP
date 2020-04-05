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
gppDF <- data.frame(seq(2012, 2023, by=1), NA, NA, NA, NA)
colnames(gppDF) <- c("Year", "gday_n_amb","gday_p_amb", "gday_n_ele", "gday_p_ele")
laiDF <- gppDF
leafDF <- gppDF
npDF <- gppDF
woodDF <- gppDF

## store data
for (i in 2012:2023) {
    # GPP
    gppDF[gppDF$Year == i, "gday_n_amb"] <- sum(de_n_DF2[de_n_DF2$year == i, "gpp"], na.rm=T)
    gppDF[gppDF$Year == i, "gday_p_amb"] <- sum(de_np_DF2[de_np_DF2$year == i, "gpp"], na.rm=T)
    gppDF[gppDF$Year == i, "gday_n_ele"] <- sum(de_n_DF3[de_n_DF3$year == i, "gpp"], na.rm=T)
    gppDF[gppDF$Year == i, "gday_p_ele"] <- sum(de_np_DF3[de_np_DF3$year == i, "gpp"], na.rm=T)
    
    # LAI
    laiDF[laiDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "lai"]
    laiDF[laiDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "lai"]
    laiDF[laiDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$year == i & de_n_DF3$doy == 1, "lai"]
    laiDF[laiDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$year == i & de_np_DF3$doy == 1, "lai"]
    
    # Leaf C
    leafDF[leafDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "shoot"]
    leafDF[leafDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "shoot"]
    leafDF[leafDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$year == i & de_n_DF3$doy == 1, "shoot"]
    leafDF[leafDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$year == i & de_np_DF3$doy == 1, "shoot"]

    
    # Leaf N:P
    npDF[npDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "shootn"]/de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "shootp"]
    npDF[npDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "shootn"]/de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "shootp"]
    npDF[npDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$year == i & de_n_DF3$doy == 1, "shootn"]/de_n_DF3[de_n_DF3$year == i & de_n_DF3$doy == 1, "shootp"]
    npDF[npDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$year == i & de_np_DF3$doy == 1, "shootn"]/de_np_DF3[de_np_DF3$year == i & de_np_DF3$doy == 1, "shootp"]
    
    # wood C
    woodDF[woodDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "stem"]
    woodDF[woodDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "stem"]
    woodDF[woodDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$year == i & de_n_DF3$doy == 1, "stem"]
    woodDF[woodDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$year == i & de_np_DF3$doy == 1, "stem"]
}


# Compute CO2 % response
gppDF[,"gday_n"] <- (gppDF[,"gday_n_ele"] - gppDF[,"gday_n_amb"])/gppDF[,"gday_n_amb"] * 100.0
laiDF[,"gday_n"] <- (laiDF[,"gday_n_ele"] - laiDF[,"gday_n_amb"])/laiDF[,"gday_n_amb"] * 100.0
leafDF[,"gday_n"] <- (leafDF[,"gday_n_ele"] - leafDF[,"gday_n_amb"])/leafDF[,"gday_n_amb"] * 100.0
npDF[,"gday_n"] <- (npDF[,"gday_n_ele"] - npDF[,"gday_n_amb"])/npDF[,"gday_n_amb"] * 100.0
woodDF[,"gday_n"] <- (woodDF[,"gday_n_ele"] - woodDF[,"gday_n_amb"])/woodDF[,"gday_n_amb"] * 100.0


gppDF[,"gday_p"] <- (gppDF[,"gday_p_ele"] - gppDF[,"gday_p_amb"])/gppDF[,"gday_p_amb"] * 100.0
laiDF[,"gday_p"] <- (laiDF[,"gday_p_ele"] - laiDF[,"gday_p_amb"])/laiDF[,"gday_p_amb"] * 100.0
leafDF[,"gday_p"] <- (leafDF[,"gday_p_ele"] - leafDF[,"gday_p_amb"])/leafDF[,"gday_p_amb"] * 100.0
npDF[,"gday_p"] <- (npDF[,"gday_p_ele"] - npDF[,"gday_p_amb"])/npDF[,"gday_p_amb"] * 100.0
woodDF[,"gday_p"] <- (woodDF[,"gday_p_ele"] - woodDF[,"gday_p_amb"])/woodDF[,"gday_p_amb"] * 100.0


pdf("R/UK_workshop_parameters.pdf", width = 9, height = 9)
m <- matrix(c(1,2,3,4,5,5),nrow = 3, ncol = 2, byrow = TRUE)

layout(mat = m,heights = c(0.4,0.4,0.2))

par(mar=c(5.1, 6.1, 3.1, 6.1),
    mgp=c(4,1,0))

# GPP
with(gppDF, plot(gday_n~Year, ylim=c(-10,30), 
                 ylab = "GPP % response",cex.axis = 1.5,
                 type="l", lwd = 3, col = "brown", cex.lab = 2.5))
with(gppDF, points(gday_p~Year, type="l", col = "red", lwd = 3))

# LAI
with(laiDF, plot(gday_n~Year, ylim=c(-10,30), 
                 ylab = "LAI % response",cex.axis = 1.5,
                 type="l", lwd = 3, col = "brown", cex.lab = 2.5))
with(laiDF, points(gday_p~Year, type="l", col = "red", lwd = 3))

# Leaf C
with(leafDF, plot(gday_n~Year, ylim=c(-10,30), 
                 ylab = "Leaf C % response",cex.axis = 1.5,
                 type="l", lwd = 3, col = "brown", cex.lab = 2.5))
with(leafDF, points(gday_p~Year, type="l", col = "red", lwd = 3))

# Leaf NP
with(npDF, plot(gday_n~Year, ylim=c(-10,10), 
                  ylab = "Leaf N:P ratio % response",cex.axis = 1.5,
                  type="l", lwd = 3, col = "brown", cex.lab = 2.5))
with(npDF, points(gday_p~Year, type="l", col = "red", lwd = 3))


plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",
       legend = c("N only", "NP model"), 
       col=c("brown", "red"), lwd=5, cex=1.5, horiz = TRUE, pt.cex = 5,
       pt.lwd = 5)

dev.off()



# Wood C
with(woodDF, plot(gday_n~Year, ylim=c(-10,30), 
                  ylab = "Wood C % response",
                  type="l", lwd = 2, col = "brown"))
with(leafDF, points(gday_p~Year, type="l", col = "red", lwd = 2))

