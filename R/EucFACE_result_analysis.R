#### EucFACE CO2 response analysis
#### This script is to accompany the UK workshop EucFACE plot
#### To provide comparison of AmazonFACE CO2 response
#### Although we are not having the same CO2 treatment

### Read in files
## DE relationship, N only
de_n_DF2 <- read.csv("outputs/EUC_amb_avg_02_ellsworth_N.csv", skip=1)
de_n_DF3 <- read.csv("outputs/EUC_ele_avg_03_ellsworth_N.csv", skip=1)

## DE relationship, N P 
de_np_DF2 <- read.csv("outputs/EUC_amb_avg_02_ellsworth_NP.csv", skip=1)
de_np_DF3 <- read.csv("outputs/EUC_ele_avg_03_ellsworth_NP.csv", skip=1)


#### Process the data to plot annual patterns
## Generate output DF at annual timestep
gppDF <- data.frame(seq(2012, 2023, by=1), NA, NA, NA, NA)
colnames(gppDF) <- c("Year", "gday_n_amb","gday_p_amb", "gday_n_ele", "gday_p_ele")
laiDF <- gppDF
nepDF <- gppDF
soilDF <- gppDF

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
    
    # nep 
    nepDF[nepDF$Year == i, "gday_n_amb"] <- sum(de_n_DF2[de_n_DF2$year == i, "nep"], na.rm=T)
    nepDF[nepDF$Year == i, "gday_p_amb"] <- sum(de_np_DF2[de_np_DF2$year == i, "nep"], na.rm=T)
    nepDF[nepDF$Year == i, "gday_n_ele"] <- sum(de_n_DF3[de_n_DF3$year == i, "nep"], na.rm=T)
    nepDF[nepDF$Year == i, "gday_p_ele"] <- sum(de_np_DF3[de_np_DF3$year == i, "nep"], na.rm=T)
    
    
    # Soil C
    soilDF[soilDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$year == i & de_n_DF2$doy == 1, "soilc"]
    soilDF[soilDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$year == i & de_np_DF2$doy == 1, "soilc"]
    soilDF[soilDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$year == i & de_n_DF3$doy == 1, "soilc"]
    soilDF[soilDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$year == i & de_np_DF3$doy == 1, "soilc"]
    
}


# Compute CO2 % response
gppDF[,"gday_n"] <- (gppDF[,"gday_n_ele"] - gppDF[,"gday_n_amb"])/gppDF[,"gday_n_amb"] * 100.0
laiDF[,"gday_n"] <- (laiDF[,"gday_n_ele"] - laiDF[,"gday_n_amb"])/laiDF[,"gday_n_amb"] * 100.0
nepDF[,"gday_n"] <- (nepDF[,"gday_n_ele"] - nepDF[,"gday_n_amb"])
soilDF[,"gday_n"] <- (soilDF[,"gday_n_ele"] - soilDF[,"gday_n_amb"])/soilDF[,"gday_n_amb"] * 100.0


gppDF[,"gday_p"] <- (gppDF[,"gday_p_ele"] - gppDF[,"gday_p_amb"])/gppDF[,"gday_p_amb"] * 100.0
laiDF[,"gday_p"] <- (laiDF[,"gday_p_ele"] - laiDF[,"gday_p_amb"])/laiDF[,"gday_p_amb"] * 100.0
nepDF[,"gday_p"] <- (nepDF[,"gday_p_ele"] - nepDF[,"gday_p_amb"])
soilDF[,"gday_p"] <- (soilDF[,"gday_p_ele"] - soilDF[,"gday_p_amb"])/soilDF[,"gday_p_amb"] * 100.0


pdf("R/UK_workshop_co2_effect_eucface.pdf", width = 9, height = 9)
m <- matrix(c(1,2,3,4,5,5),nrow = 3, ncol = 2, byrow = TRUE)

layout(mat = m,heights = c(0.4,0.4,0.2))

par(mar=c(5.1, 6.1, 3.1, 6.1),
    mgp=c(4,1,0))

# GPP
with(gppDF, plot(gday_n~Year, ylim=c(-10,30), 
                 ylab = "GPP % response",cex.axis = 1.5,
                 type="b", lwd = 3, col = "brown", cex.lab = 2.5))
with(gppDF, points(gday_p~Year, type="b", col = "red", lwd = 3))

# LAI
with(laiDF, plot(gday_n~Year, ylim=c(-10,30), 
                 ylab = "LAI % response",cex.axis = 1.5,
                 type="b", lwd = 3, col = "brown", cex.lab = 2.5))
with(laiDF, points(gday_p~Year, type="b", col = "red", lwd = 3))

# NEP
with(nepDF, plot(gday_n*100~Year, ylim=c(-100,500), 
                  ylab = expression(paste("NEP response [g ", m^-2, " ", yr^-1, "]")),cex.axis = 1.5,
                  type="b", lwd = 3, col = "brown", cex.lab = 2))
with(nepDF, points(gday_p*100~Year, type="b", col = "red", lwd = 3))

# Soil C
with(soilDF, plot(gday_n~Year, ylim=c(-2,2), 
                ylab = "Soil C % response",cex.axis = 1.5,
                type="b", lwd = 3, col = "brown", cex.lab = 2.5))
with(soilDF, points(gday_p~Year, type="b", col = "red", lwd = 3))


plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",
       legend = c("N only", "NP model"), 
       col=c("brown", "red"), lwd=5, cex=1.5, horiz = TRUE, pt.cex = 5,
       pt.lwd = 5)

dev.off()


### Time series plot
pdf("R/UK_workshop_temporal_eucface.pdf", width = 9, height = 9)
m <- matrix(c(1,2,3,4,5,5),nrow = 3, ncol = 2, byrow = TRUE)

layout(mat = m,heights = c(0.4,0.4,0.2))

par(mar=c(5.1, 6.1, 3.1, 6.1),
    mgp=c(3,1,0))

# GPP
with(gppDF, plot(gday_n_ele/10~Year, ylim=c(0,5), 
                 ylab = expression(paste("GPP [kg ", m^-2, " ", yr^-1, "]")),cex.axis = 1.5,
                 type="b", lwd = 3, col = "brown", cex.lab = 2.5))
with(gppDF, points(gday_p_ele/10~Year, type="b", col = "red", lwd = 3))

# LAI
with(laiDF, plot(gday_n_ele~Year, ylim=c(0, 8), 
                 ylab = "LAI",cex.axis = 1.5,
                 type="b", lwd = 3, col = "brown", cex.lab = 2.5))
with(laiDF, points(gday_p_ele~Year, type="b", col = "red", lwd = 3))

# NEP
with(nepDF, plot(gday_n_ele/10~Year, ylim=c(-1,1), 
                 ylab = expression(paste("NEP [kg ", m^-2, " ", yr^-1, "]")),cex.axis = 1.5,
                 type="b", lwd = 3, col = "brown", cex.lab = 2))
with(nepDF, points(gday_p_ele/10~Year, type="b", col = "red", lwd = 3))

# Soil C
with(soilDF, plot(gday_n_ele/10~Year, ylim=c(0,12), 
                  ylab = expression(paste("Soil C [kg ", m^-2, "]")),cex.axis = 1.5,
                  type="b", lwd = 3, col = "brown", cex.lab = 2.5))
with(soilDF, points(gday_p_ele/10~Year, type="b", col = "red", lwd = 3))


plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",
       legend = c("N only", "NP model"), 
       col=c("brown", "red"), lwd=5, cex=1.5, horiz = TRUE, pt.cex = 5,
       pt.lwd = 5)

dev.off()

