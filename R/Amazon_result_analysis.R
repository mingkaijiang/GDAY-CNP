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
    soilDF[leafDF$Year == i, "gday_n_amb"] <- de_n_DF2[de_n_DF2$YEAR == i & de_n_DF2$DOY == 1, "CSOIL"]
    soilDF[leafDF$Year == i, "gday_p_amb"] <- de_np_DF2[de_np_DF2$YEAR == i & de_np_DF2$DOY == 1, "CSOIL"]
    soilDF[leafDF$Year == i, "gday_n_ele"] <- de_n_DF3[de_n_DF3$YEAR == i & de_n_DF3$DOY == 1, "CSOIL"]
    soilDF[leafDF$Year == i, "gday_p_ele"] <- de_np_DF3[de_np_DF3$YEAR == i & de_np_DF3$DOY == 1, "CSOIL"]
}


# Compute CO2 % response
gppDF[,"gday_n"] <- (gppDF[,"gday_n_ele"] - gppDF[,"gday_n_amb"])/gppDF[,"gday_n_amb"] * 100.0
laiDF[,"gday_n"] <- (laiDF[,"gday_n_ele"] - laiDF[,"gday_n_amb"])/laiDF[,"gday_n_amb"] * 100.0
leafDF[,"gday_n"] <- (leafDF[,"gday_n_ele"] - leafDF[,"gday_n_amb"])/leafDF[,"gday_n_amb"] * 100.0
npDF[,"gday_n"] <- (npDF[,"gday_n_ele"] - npDF[,"gday_n_amb"])/npDF[,"gday_n_amb"] * 100.0
woodDF[,"gday_n"] <- (woodDF[,"gday_n_ele"] - woodDF[,"gday_n_amb"])/woodDF[,"gday_n_amb"] * 100.0
nepDF[,"gday_n"] <- (nepDF[,"gday_n_ele"] - nepDF[,"gday_n_amb"])
soilDF[,"gday_n"] <- (soilDF[,"gday_n_ele"] - soilDF[,"gday_n_amb"])/soilDF[,"gday_n_amb"] * 100.0


gppDF[,"gday_p"] <- (gppDF[,"gday_p_ele"] - gppDF[,"gday_p_amb"])/gppDF[,"gday_p_amb"] * 100.0
laiDF[,"gday_p"] <- (laiDF[,"gday_p_ele"] - laiDF[,"gday_p_amb"])/laiDF[,"gday_p_amb"] * 100.0
leafDF[,"gday_p"] <- (leafDF[,"gday_p_ele"] - leafDF[,"gday_p_amb"])/leafDF[,"gday_p_amb"] * 100.0
npDF[,"gday_p"] <- (npDF[,"gday_p_ele"] - npDF[,"gday_p_amb"])/npDF[,"gday_p_amb"] * 100.0
woodDF[,"gday_p"] <- (woodDF[,"gday_p_ele"] - woodDF[,"gday_p_amb"])/woodDF[,"gday_p_amb"] * 100.0
nepDF[,"gday_p"] <- (nepDF[,"gday_p_ele"] - nepDF[,"gday_p_amb"])
soilDF[,"gday_p"] <- (soilDF[,"gday_p_ele"] - soilDF[,"gday_p_amb"])/soilDF[,"gday_p_amb"] * 100.0



pdf("R/UK_workshop_co2_effect_amazon.pdf", width = 9, height = 9)
m <- matrix(c(1,2,3,4,5,5),nrow = 3, ncol = 2, byrow = TRUE)

layout(mat = m,heights = c(0.4,0.4,0.2))

par(mar=c(5.1, 6.1, 3.1, 6.1),
    mgp=c(4,1,0))

# GPP

with(gppDF, plot(gday_n~Year, ylim=c(-10,30), 
                 ylab = "GPP % response",cex.axis = 1.5,
                 type="b", lwd = 3, col = "green", cex.lab = 2.5))
with(gppDF, points(gday_p~Year, type="b", col = "red", lwd = 3))
x<-par("usr")
rect(x[1],x[3],x[2],x[4],col=adjustcolor("lightgray", 0.2))
grid(lty=6, col="white")

# LAI
with(laiDF, plot(gday_n~Year, ylim=c(-10,30), 
                 ylab = "LAI % response",cex.axis = 1.5,
                 type="b", lwd = 3, col = "green", cex.lab = 2.5))
with(laiDF, points(gday_p~Year, type="b", col = "red", lwd = 3))
x<-par("usr")
rect(x[1],x[3],x[2],x[4],col=adjustcolor("lightgray", 0.2))
grid(lty=6, col="white")

# NEP
with(nepDF, plot(gday_n~Year, ylim=c(-100,500), 
                 ylab = expression(paste("NEP response [g ", m^-2, " ", yr^-1, "]")),cex.axis = 1.5,
                 type="b", lwd = 3, col = "green", cex.lab = 2))
with(nepDF, points(gday_p~Year, type="b", col = "red", lwd = 3))
x<-par("usr")
rect(x[1],x[3],x[2],x[4],col=adjustcolor("lightgray", 0.2))
grid(lty=6, col="white")

# Soil C
with(soilDF, plot(gday_n~Year, ylim=c(-2,2), 
                  ylab = "Soil C % response",cex.axis = 1.5,
                  type="b", lwd = 3, col = "green", cex.lab = 2.5))
with(soilDF, points(gday_p~Year, type="b", col = "red", lwd = 3))
x<-par("usr")
rect(x[1],x[3],x[2],x[4],col=adjustcolor("lightgray", 0.2))
grid(lty=6, col="white")


plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",
       legend = c("N only", "NP model"), 
       col=c("green", "red"), lwd=5, cex=1.5, horiz = TRUE, pt.cex = 5,
       pt.lwd = 5)

dev.off()

### Time series plot
pdf("R/UK_workshop_temporal_amazon.pdf", width = 9, height = 9)
m <- matrix(c(1,2,3,4,5,5),nrow = 3, ncol = 2, byrow = TRUE)

layout(mat = m,heights = c(0.4,0.4,0.2))

par(mar=c(5.1, 6.1, 3.1, 6.1),
    mgp=c(3,1,0))

# GPP
with(gppDF, plot(gday_n_ele/1000~Year, ylim=c(0,5), 
                 ylab = expression(paste("GPP [kg ", m^-2, " ", yr^-1, "]")),cex.axis = 1.5,
                 type="b", lwd = 3, col = "green", cex.lab = 2.5))
with(gppDF, points(gday_p_ele/1000~Year, type="b", col = "red", lwd = 3))
x<-par("usr")
rect(x[1],x[3],x[2],x[4],col=adjustcolor("lightgray", 0.2))
grid(lty=6, col="white")

# LAI
with(laiDF, plot(gday_n_ele~Year, ylim=c(0, 8), 
                 ylab = "LAI",cex.axis = 1.5,
                 type="b", lwd = 3, col = "green", cex.lab = 2.5))
with(laiDF, points(gday_p_ele~Year, type="b", col = "red", lwd = 3))
x<-par("usr")
rect(x[1],x[3],x[2],x[4],col=adjustcolor("lightgray", 0.2))
grid(lty=6, col="white")

# NEP
with(nepDF, plot(gday_n_ele/1000~Year, ylim=c(-1,1), 
                 ylab = expression(paste("NEP [kg ", m^-2, " ", yr^-1, "]")),cex.axis = 1.5,
                 type="b", lwd = 3, col = "green", cex.lab = 2))
with(nepDF, points(gday_p_ele/1000~Year, type="b", col = "red", lwd = 3))
x<-par("usr")
rect(x[1],x[3],x[2],x[4],col=adjustcolor("lightgray", 0.2))
grid(lty=6, col="white")

# Soil C
with(soilDF, plot(gday_n_ele/1000~Year, ylim=c(0,12), 
                  ylab = expression(paste("Soil C [kg ", m^-2, "]")),cex.axis = 1.5,
                  type="b", lwd = 3, col = "green", cex.lab = 2.5))
with(soilDF, points(gday_p_ele/1000~Year, type="b", col = "red", lwd = 3))
x<-par("usr")
rect(x[1],x[3],x[2],x[4],col=adjustcolor("lightgray", 0.2))
grid(lty=6, col="white")

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",
       legend = c("N only", "NP model"), 
       col=c("green", "red"), lwd=5, cex=1.5, horiz = TRUE, pt.cex = 5,
       pt.lwd = 5)

dev.off()
