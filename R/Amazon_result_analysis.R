#### AmazonFACE CO2 response analysis
#### This script is to accompany the UK workshop EucFACE plot
#### To provide comparison of AmazonFACE CO2 response
#### Although we are not having the same CO2 treatment

### Read in files
## DE relationship, N only
de_n_DF2 <- read.csv("~/Documents/Research/Projects/Amazon/AMAZ/drought/outputs/AmaFACE1_D_GDA_AMB_1999_2023.csv", skip=3)
de_n_DF3 <- read.csv("~/Documents/Research/Projects/Amazon/AMAZ/drought/outputs/AmaFACE1_D_GDA_ELE_1999_2023.csv", skip=3)

## DE relationship, N P 
de_np_DF2 <- read.csv("~/Documents/Research/Projects/Amazon/AMAZ/drought_p/outputs/AmaFACE1_D_GDA_AMB_1999_2023.csv", skip=3)
de_np_DF3 <- read.csv("~/Documents/Research/Projects/Amazon/AMAZ/drought_p/outputs/AmaFACE1_D_GDA_ELE_1999_2023.csv", skip=3)


#### Process the data to plot annual patterns
## Generate output DF at annual timestep
gppDF <- data.frame(seq(1999, 2023, by=1), NA, NA, NA, NA)
colnames(gppDF) <- c("Year", "gday_n_amb","gday_p_amb", "gday_n_ele", "gday_p_ele")
laiDF <- gppDF
leafDF <- gppDF
npDF <- gppDF
woodDF <- gppDF
nepDF <- gppDF
soilDF <- gppDF

## store data
for (i in 1999:2023) {
    # GPP
    gppDF[gppDF$Year == i, "gday_n_amb"] <- sum(de_n_DF2[de_n_DF2$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "gday_p_amb"] <- sum(de_np_DF2[de_np_DF2$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "gday_n_ele"] <- sum(de_n_DF3[de_n_DF3$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "gday_p_ele"] <- sum(de_np_DF3[de_np_DF3$YEAR == i, "GPP"], na.rm=T)

    # LAI
    laiDF[laiDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$YEAR == i & de_n_DF2$DOY == 1, "LAI"]
    laiDF[laiDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$YEAR == i & de_np_DF2$DOY == 1, "LAI"]
    laiDF[laiDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$YEAR == i & de_n_DF3$DOY == 1, "LAI"]
    laiDF[laiDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$YEAR == i & de_np_DF3$DOY == 1, "LAI"]
    
    # Leaf C
    leafDF[leafDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$YEAR == i & de_n_DF2$DOY == 1, "CL"]
    leafDF[leafDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$YEAR == i & de_np_DF2$DOY == 1, "CL"]
    leafDF[leafDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$YEAR == i & de_n_DF3$DOY == 1, "CL"]
    leafDF[leafDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$YEAR == i & de_np_DF3$DOY == 1, "CL"]
    
    
    # Leaf N:P
    npDF[npDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$YEAR == i & de_n_DF2$DOY == 1, "NL"]/de_n_DF2[de_n_DF2$YEAR == i & de_n_DF2$DOY == 1, "PL"]
    npDF[npDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$YEAR == i & de_np_DF2$DOY == 1, "NL"]/de_np_DF2[de_np_DF2$YEAR == i & de_np_DF2$DOY == 1, "PL"]
    npDF[npDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$YEAR == i & de_n_DF3$DOY == 1, "NL"]/de_n_DF3[de_n_DF3$YEAR == i & de_n_DF3$DOY == 1, "PL"]
    npDF[npDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$YEAR == i & de_np_DF3$DOY == 1, "NL"]/de_np_DF3[de_np_DF3$YEAR == i & de_np_DF3$DOY == 1, "PL"]
    
    # NEP
    nepDF[gppDF$Year == i, "gday_n_amb"] <- sum(de_n_DF2[de_n_DF2$YEAR == i, "NEP"], na.rm=T)
    nepDF[gppDF$Year == i, "gday_p_amb"] <- sum(de_np_DF2[de_np_DF2$YEAR == i, "NEP"], na.rm=T)
    nepDF[gppDF$Year == i, "gday_n_ele"] <- sum(de_n_DF3[de_n_DF3$YEAR == i, "NEP"], na.rm=T)
    nepDF[gppDF$Year == i, "gday_p_ele"] <- sum(de_np_DF3[de_np_DF3$YEAR == i, "NEP"], na.rm=T)
    
    # SOil C
    soilDF[leafDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$YEAR == i & de_n_DF2$DOY == 1, "CL"]
    soilDF[leafDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$YEAR == i & de_np_DF2$DOY == 1, "CL"]
    soilDF[leafDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$YEAR == i & de_n_DF3$DOY == 1, "CL"]
    soilDF[leafDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$YEAR == i & de_np_DF3$DOY == 1, "CL"]
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



pdf("R/UK_workshop_parameters_amazon.pdf", width = 9, height = 9)
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



### Will need to subtract met data for 2012-2023 for AmazonFACE and run model just for these periods using CO2 data from EucFACE
### So that we are consistently comparing the same CO2 effect

### TO do list:
### 1. run AmazonFACE
### 2. run quasi-equil framework with parameters from EucFACE and AmazonFACE, plotting
### 3. Plot inst, L and VL constraint comparisons.