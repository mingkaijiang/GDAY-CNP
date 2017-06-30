#### To reproduce figures in Medlyn et al. 2016 GCB Figure 3 panel a-d

#### Read in all model outputs
#### They are downloaded according to paper supplementary materials
#### Only the fixed met data output was used

### amb CO2 
cabl_amb <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/CABL/outputs/D1CABLEUCAMBAVG.csv",
                 skip=7)
clm4_amb <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/CLM4/outputs/D1CLM4EUCAMBAVG.csv",
                 skip=0)
clmp_amb <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/CLMP/outputs/D1CLMPEUCAMBAVG.csv",
                 skip=0)
gday_amb <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/GDAY/outputs/D1GDAYEUCAMBAVG.csv",
                 skip=3)
#gday_amb <- read.csv("outputs/EUC_amb_avg_02_ellsworth_N.csv", skip=1)
lpjg_amb <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/LPJG/outputs/D1LPJXEUCAMBAVG.csv",
                 skip=2)
ocnx_amb <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/OCNX/outputs/D1OCNXEUCAMBAVG.csv",
                 skip=0)
sdvm_amb <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/SDVM/outputs/D1SDVMEUCAMBAVG.csv",
                 skip=0)
gdap_amb <- read.csv("outputs/EUC_amb_avg_02_ellsworth_NP.csv", skip=1)



## elev CO2
cabl_ele <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/CABL/outputs/D1CABLEUCELEAVG.csv",
                     skip=7)
clm4_ele <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/CLM4/outputs/D1CLM4EUCELEAVG.csv",
                     skip=0)
clmp_ele <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/CLMP/outputs/D1CLMPEUCELEAVG.csv",
                     skip=0)
gday_ele <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/GDAY/outputs/D1GDAYEUCELEAVG.csv",
                     skip=3)
#gday_ele <- read.csv("outputs/EUC_ele_avg_03_ellsworth_N.csv", skip=1)
lpjg_ele <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/LPJG/outputs/D1LPJXEUCELEAVG.csv",
                     skip=2)
ocnx_ele <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/OCNX/outputs/D1OCNXEUCELEAVG.csv",
                     skip=0)
sdvm_ele <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/SDVM/outputs/D1SDVMEUCELEAVG.csv",
                     skip=0)
gdap_ele <- read.csv("outputs/EUC_ele_avg_03_ellsworth_NP.csv", skip=1)

## exclude day 366 from all data
cabl_amb <- subset(cabl_amb, DOY<=365)
clm4_amb <- subset(clm4_amb, DOY<=365)
clmp_amb <- subset(clmp_amb, DOY<=365)
gday_amb <- subset(gday_amb, DOY<=365)
lpjg_amb <- subset(lpjg_amb, DOY<=365)
ocnx_amb <- subset(ocnx_amb, DOY<=365)
sdvm_amb <- subset(sdvm_amb, DOY<=365)
gdap_amb <- subset(gdap_amb, doy<=365)
cabl_ele <- subset(cabl_ele, DOY<=365)
clm4_ele <- subset(clm4_ele, DOY<=365)
clmp_ele <- subset(clmp_ele, DOY<=365)
gday_ele <- subset(gday_ele, DOY<=365)
lpjg_ele <- subset(lpjg_ele, DOY<=365)
ocnx_ele <- subset(ocnx_ele, DOY<=365)
sdvm_ele <- subset(sdvm_ele, DOY<=365)
gdap_ele <- subset(gdap_ele, doy<=365)

## convert -9999 to missing values
cabl_amb[cabl_amb == -9999] <- NA
clm4_amb[clm4_amb == -9999] <- NA
clmp_amb[clmp_amb == -9999] <- NA
gday_amb[gday_amb == -9999] <- NA
lpjg_amb[lpjg_amb == -9999] <- NA
ocnx_amb[ocnx_amb == -9999] <- NA
sdvm_amb[sdvm_amb == -9999] <- NA
gdap_amb[gdap_amb == -9999] <- NA

cabl_ele[cabl_ele == -9999] <- NA
clm4_ele[clm4_ele == -9999] <- NA
clmp_ele[clmp_ele == -9999] <- NA
gday_ele[gday_ele == -9999] <- NA
lpjg_ele[lpjg_ele == -9999] <- NA
ocnx_ele[ocnx_ele == -9999] <- NA
sdvm_ele[sdvm_ele == -9999] <- NA
gdap_ele[gdap_ele == -9999] <- NA


## Generate output DF at annual timestep
gppDF <- data.frame(seq(2012, 2023, by=1), NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA)
colnames(gppDF) <- c("Year", "cabl_amb", "clm4_amb", "clmp_amb",
                     "gday_amb", "lpjg_amb", "ocnx_amb", "sdvm_amb", "gdap_amb",
                     "cabl_ele", "clm4_ele", "clmp_ele",
                     "gday_ele", "lpjg_ele", "ocnx_ele", "sdvm_ele", "gdap_ele")
nppDF <- gppDF
cueDF <- gppDF
laiDF <- gppDF

## store data
for (i in 2012:2023) {
    # GPP
    gppDF[gppDF$Year == i, "cabl_amb"] <- sum(cabl_amb[cabl_amb$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "clm4_amb"] <- sum(clm4_amb[clm4_amb$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "clmp_amb"] <- sum(clmp_amb[clmp_amb$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "gday_amb"] <- sum(gday_amb[gday_amb$YEAR == i, "GPP"], na.rm=T)
    #gppDF[gppDF$Year == i, "gday_amb"] <- sum(gday_amb[gday_amb$year == i, "gpp"], na.rm=T)
    gppDF[gppDF$Year == i, "lpjg_amb"] <- sum(lpjg_amb[lpjg_amb$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "ocnx_amb"] <- sum(ocnx_amb[ocnx_amb$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "sdvm_amb"] <- sum(sdvm_amb[sdvm_amb$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "gdap_amb"] <- sum(gdap_amb[gdap_amb$year == i, "gpp"], na.rm=T)
    
    
    gppDF[gppDF$Year == i, "cabl_ele"] <- sum(cabl_ele[cabl_ele$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "clm4_ele"] <- sum(clm4_ele[clm4_ele$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "clmp_ele"] <- sum(clmp_ele[clmp_ele$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "gday_ele"] <- sum(gday_ele[gday_ele$YEAR == i, "GPP"], na.rm=T)
    #gppDF[gppDF$Year == i, "gday_ele"] <- sum(gday_ele[gday_ele$year == i, "gpp"], na.rm=T)
    gppDF[gppDF$Year == i, "lpjg_ele"] <- sum(lpjg_ele[lpjg_ele$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "ocnx_ele"] <- sum(ocnx_ele[ocnx_ele$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "sdvm_ele"] <- sum(sdvm_ele[sdvm_ele$YEAR == i, "GPP"], na.rm=T)
    gppDF[gppDF$Year == i, "gdap_ele"] <- sum(gdap_ele[gdap_ele$year == i, "gpp"], na.rm=T)
    
    
    # NPP
    nppDF[nppDF$Year == i, "cabl_amb"] <- sum(cabl_amb[cabl_amb$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "clm4_amb"] <- sum(clm4_amb[clm4_amb$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "clmp_amb"] <- sum(clmp_amb[clmp_amb$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "gday_amb"] <- sum(gday_amb[gday_amb$YEAR == i, "NPP"], na.rm=T)
    #nppDF[nppDF$Year == i, "gday_amb"] <- sum(gday_amb[gday_amb$year == i, "npp"], na.rm=T)
    
    nppDF[nppDF$Year == i, "lpjg_amb"] <- sum(lpjg_amb[lpjg_amb$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "ocnx_amb"] <- sum(ocnx_amb[ocnx_amb$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "sdvm_amb"] <- sum(sdvm_amb[sdvm_amb$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "gdap_amb"] <- sum(gdap_amb[gdap_amb$year == i, "npp"], na.rm=T)
    
    
    nppDF[nppDF$Year == i, "cabl_ele"] <- sum(cabl_ele[cabl_ele$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "clm4_ele"] <- sum(clm4_ele[clm4_ele$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "clmp_ele"] <- sum(clmp_ele[clmp_ele$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "gday_ele"] <- sum(gday_ele[gday_ele$YEAR == i, "NPP"], na.rm=T)
    #nppDF[nppDF$Year == i, "gday_ele"] <- sum(gday_ele[gday_ele$year == i, "npp"], na.rm=T)
    
    nppDF[nppDF$Year == i, "lpjg_ele"] <- sum(lpjg_ele[lpjg_ele$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "ocnx_ele"] <- sum(ocnx_ele[ocnx_ele$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "sdvm_ele"] <- sum(sdvm_ele[sdvm_ele$YEAR == i, "NPP"], na.rm=T)
    nppDF[nppDF$Year == i, "gdap_ele"] <- sum(gdap_ele[gdap_ele$year == i, "npp"], na.rm=T)
    
    
    # LAI
    laiDF[laiDF$Year == i, "cabl_amb"] <- cabl_amb[cabl_amb$YEAR == i & cabl_amb$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "clm4_amb"] <- clm4_amb[clm4_amb$YEAR == i & clm4_amb$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "clmp_amb"] <- clmp_amb[clmp_amb$YEAR == i & clmp_amb$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "gday_amb"] <- gday_amb[gday_amb$YEAR == i & gday_amb$DOY == 365, "LAI"]
    #laiDF[laiDF$Year == i, "gday_amb"] <- gday_amb[gday_amb$year == i & gday_amb$doy == 365, "lai"]
    
    laiDF[laiDF$Year == i, "lpjg_amb"] <- lpjg_amb[lpjg_amb$YEAR == i & lpjg_amb$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "ocnx_amb"] <- ocnx_amb[ocnx_amb$YEAR == i & ocnx_amb$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "sdvm_amb"] <- sdvm_amb[sdvm_amb$YEAR == i & sdvm_amb$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "gdap_amb"] <- gdap_amb[gdap_amb$year == i & gdap_amb$doy == 365, "lai"]
    
    
    laiDF[laiDF$Year == i, "cabl_ele"] <- cabl_ele[cabl_ele$YEAR == i & cabl_ele$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "clm4_ele"] <- clm4_ele[clm4_ele$YEAR == i & clm4_ele$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "clmp_ele"] <- clmp_ele[clmp_ele$YEAR == i & clmp_ele$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "gday_ele"] <- gday_ele[gday_ele$YEAR == i & gday_ele$DOY == 365, "LAI"]
    #laiDF[laiDF$Year == i, "gday_ele"] <- gday_ele[gday_ele$year == i & gday_ele$doy == 365, "lai"]
    
    laiDF[laiDF$Year == i, "lpjg_ele"] <- lpjg_ele[lpjg_ele$YEAR == i & lpjg_ele$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "ocnx_ele"] <- ocnx_ele[ocnx_ele$YEAR == i & ocnx_ele$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "sdvm_ele"] <- sdvm_ele[sdvm_ele$YEAR == i & sdvm_ele$DOY == 365, "LAI"]
    laiDF[laiDF$Year == i, "gdap_ele"] <- gdap_ele[gdap_ele$year == i & gdap_ele$doy == 365, "lai"]
    
                                                                                        
}

## Get cue
cueDF[,2:ncol(cueDF)] <- nppDF[,2:ncol(cueDF)]/gppDF[,2:ncol(cueDF)] 

## compute % CO2 response (ele at each year - amb in 2012): anomaly 1
gppRSP1 <- data.frame(seq(2012, 2023, by=1), NA, NA, NA, NA, NA, NA, NA, NA)
colnames(gppRSP1) <- c("Year", "cabl", "clm4", "clmp",
                     "gday", "lpjg", "ocnx", "sdvm", "gdap")
nppRSP1 <- gppRSP1
cueRSP1 <- gppRSP1
laiRSP1 <- gppRSP1

for (i in 1:nrow(gppRSP1)) {
    gppRSP1[i,2:9] <- (gppDF[i,10:17] - gppDF[1,2:9])/gppDF[1,2:9] * 100.0
    nppRSP1[i,2:9] <- (nppDF[i,10:17] - nppDF[1,2:9])/nppDF[1,2:9] * 100.0
    cueRSP1[i,2:9] <- (cueDF[i,10:17] - cueDF[1,2:9])/cueDF[1,2:9] * 100.0
    laiRSP1[i,2:9] <- (laiDF[i,10:17] - laiDF[1,2:9])/laiDF[1,2:9] * 100.0
}

## compute % CO2 response (ele at each year - amb at each year): anomaly 2
gppRSP2 <- gppRSP1
nppRSP2 <- gppRSP1
cueRSP2 <- gppRSP1
laiRSP2 <- gppRSP1

gppRSP2[,2:9] <- (gppDF[,10:17] - gppDF[,2:9])/gppDF[,2:9] * 100.0
nppRSP2[,2:9] <- (nppDF[,10:17] - nppDF[,2:9])/nppDF[,2:9] * 100.0
cueRSP2[,2:9] <- (cueDF[,10:17] - cueDF[,2:9])/cueDF[,2:9] * 100.0
laiRSP2[,2:9] <- (laiDF[,10:17] - laiDF[,2:9])/laiDF[,2:9] * 100.0

## Plot anomaly 1
with(gppRSP1, plot(cabl~Year, ylim = c(0, 30), type="l", lwd = 2,
                   ylab = "Response of GPP (%)", lty = 2, col = "darkgreen"))
with(gppRSP1, points(clm4~Year, type="l", lty = 1, col = "blue", lwd = 2))
with(gppRSP1, points(clmp~Year, type="l", lty = 1, col = "darkgreen", lwd = 2))
with(gppRSP1, points(gday~Year, type="l", lty = 2, col = "blue", lwd = 2))
with(gppRSP1, points(lpjg~Year, type="l", lty = 4, col = "blue", lwd = 2))
with(gppRSP1, points(ocnx~Year, type="l", lty = 3, col = "blue", lwd = 2))
with(gppRSP1, points(sdvm~Year, type="l", lty = 1, col = "red", lwd = 2))
with(gppRSP1, points(gdap~Year, type="l", lty = 1, col = "orange", lwd = 2.5))

# Plot graphics
col.list <- c('#4DAF4A','#377EB8','#4DAF4A', '#377EB8', '#377EB8',
              '#377EB8','#fc8d62')


## Plot anomaly 2
pdf('R/Medlyn_2016_Figure3_reproduce.pdf', width=8, height=6)
par(mfrow=c(2,2), mar=c(2.1, 4.1, 2.1, 4.1),
    mgp=c(2,1,0))

# GPP
with(gppRSP2, plot(cabl~Year, ylim = c(-10, 30), type="l", lwd = 1,
                   ylab = "Response of GPP (%)", lty = 2, col = "darkgreen"))
with(gppRSP2, points(clm4~Year, type="l", lty = 2, col = "blue", lwd = 1))
with(gppRSP2, points(clmp~Year, type="l", lty = 2, col = "orange", lwd = 1))
with(gppRSP2, points(gday~Year, type="l", lty = 1, col = "black", lwd = 2))
with(gppRSP2, points(lpjg~Year, type="l", lty = 2, col = "lightblue", lwd = 1))
with(gppRSP2, points(ocnx~Year, type="l", lty = 2, col = "green", lwd = 1))
with(gppRSP2, points(sdvm~Year, type="l", lty = 2, col = "lightgreen", lwd = 1))
with(gppRSP2, points(gdap~Year, type="l", lty = 1, col = "red", lwd = 2.5))



# NPP
with(nppRSP2, plot(cabl~Year, ylim = c(-10, 40), type="l", lwd = 1,
                   ylab = "Response of NPP (%)", lty = 2, col = "darkgreen"))
with(nppRSP2, points(clm4~Year, type="l", lty = 2, col = "blue", lwd = 1))
with(nppRSP2, points(clmp~Year, type="l", lty = 2, col = "orange", lwd = 1))
with(nppRSP2, points(gday~Year, type="l", lty = 1, col = "black", lwd = 2))
with(nppRSP2, points(lpjg~Year, type="l", lty = 2, col = "lightblue", lwd = 1))
with(nppRSP2, points(ocnx~Year, type="l", lty = 2, col = "green", lwd = 1))
with(nppRSP2, points(sdvm~Year, type="l", lty = 2, col = "lightgreen", lwd = 1))
with(nppRSP2, points(gdap~Year, type="l", lty = 1, col = "red", lwd = 2.5))

# CUE
with(cueRSP2, plot(cabl~Year, ylim = c(-20, 50), type="l", lwd = 1,
                   ylab = "Response of CUE (%)", lty = 2, col = "darkgreen"))
with(cueRSP2, points(clm4~Year, type="l", lty = 2, col = "blue", lwd = 1))
with(cueRSP2, points(clmp~Year, type="l", lty = 2, col = "orange", lwd = 1))
with(cueRSP2, points(gday~Year, type="l", lty = 1, col = "black", lwd = 2))
with(cueRSP2, points(lpjg~Year, type="l", lty = 2, col = "lightblue", lwd = 1))
with(cueRSP2, points(ocnx~Year, type="l", lty = 2, col = "green", lwd = 1))
with(cueRSP2, points(sdvm~Year, type="l", lty = 2, col = "lightgreen", lwd = 1))
with(cueRSP2, points(gdap~Year, type="l", lty = 1, col = "red", lwd = 2.5))
#legend("topright", c("CABLE", "CLM4", "CLM-P", "GDAY", "LPJ", "OCN", "SDVM", "GDAY-P"),
#       lty = c(2, 2, 2, 1, 2, 2, 2, 1), col = c("darkgreen", "blue", "orange", "black",
#                                                "lightblue", "green", "lightgreen", "red"),
#        cex = 0.5, bg = "white")

# LAI
with(laiRSP2, plot(cabl~Year, ylim = c(-10, 20), type="l", lwd = 1,
                   ylab = "Response of LAI (%)", lty = 2, col = "darkgreen"))
with(laiRSP2, points(clm4~Year, type="l", lty = 2, col = "blue", lwd = 1))
with(laiRSP2, points(clmp~Year, type="l", lty = 2, col = "orange", lwd = 1))
with(laiRSP2, points(gday~Year, type="l", lty = 1, col = "black", lwd = 2))
with(laiRSP2, points(lpjg~Year, type="l", lty = 2, col = "lightblue", lwd = 1))
with(laiRSP2, points(ocnx~Year, type="l", lty = 2, col = "green", lwd = 1))
with(laiRSP2, points(sdvm~Year, type="l", lty = 2, col = "lightgreen", lwd = 1))
with(laiRSP2, points(gdap~Year, type="l", lty = 1, col = "red", lwd = 2.5))

dev.off()

## Plot annual pattern
pdf('R/annual_pattern.pdf')
par(mfrow=c(2,1), mar=c(2.1, 3.1, 2.1, 3.1),
    mgp=c(2,1,0))

# GPP
with(gppDF, plot(cabl_amb~Year, ylim = c(0, 3000), type="l", lwd = 1,
                 ylab = "GPP [g C m-2 yr-1]", lty = 2, col = "darkgreen"))
with(gppDF, points(clm4_amb~Year, type="l", lty = 2, col = "blue", lwd = 1))
with(gppDF, points(clmp_amb~Year, type="l", lty = 2, col = "orange", lwd = 1))
with(gppDF, points(gday_amb~Year, type="l", lty = 1, col = "black", lwd = 2))
with(gppDF, points(lpjg_amb~Year, type="l", lty = 2, col = "lightblue", lwd = 1))
with(gppDF, points(ocnx_amb~Year, type="l", lty = 2, col = "green", lwd = 1))
with(gppDF, points(sdvm_amb~Year, type="l", lty = 2, col = "lightgreen", lwd = 1))
with(gppDF, points(gdap_amb*100~Year, type="l", lty = 1, col = "red", lwd = 2.5))
#legend("bottomright", c("CABLE", "CLM4", "CLM-P", "GDAY", "LPJ", "OCN", "SDVM", "GDAY-P"),
#       lty = c(2, 1, 1, 2, 4, 3, 1, 1), col = c("darkgreen", "blue", "darkgreen", "blue",
#                                                "blue", "blue", "red", "orange"),
#       cex = 0.7, bg = "white")

# LAI
with(laiDF, plot(cabl_amb~Year, ylim = c(1, 5), type="l", lwd = 1,
                 ylab = "LAI", lty = 2, col = "darkgreen"))
with(laiDF, points(clm4_amb~Year, type="l", lty = 2, col = "blue", lwd = 1))
with(laiDF, points(clmp_amb~Year, type="l", lty = 2, col = "orange", lwd = 1))
with(laiDF, points(gday_amb~Year, type="l", lty = 1, col = "black", lwd = 2))
with(laiDF, points(lpjg_amb~Year, type="l", lty = 2, col = "lightblue", lwd = 1))
with(laiDF, points(ocnx_amb~Year, type="l", lty = 2, col = "green", lwd = 1))
with(laiDF, points(sdvm_amb~Year, type="l", lty = 2, col = "lightgreen", lwd = 1))
with(laiDF, points(gdap_amb~Year, type="l", lty = 1, col = "red", lwd = 2.5))


dev.off()


## Plot annual pattern
pdf('R/annual_pattern_amb_vs_ele.pdf')
par(mfrow=c(2,1), mar=c(2.1, 3.1, 2.1, 3.1),
    mgp=c(2,1,0))

# GPP
with(gppDF, plot(gday_amb~Year, ylim = c(1000, 1800), type="l", lwd = 1.5,
                 ylab = "GPP [g C m-2 yr-1]", lty = 2, col = "darkgreen"))
with(gppDF, points(gday_ele~Year, type="l", lty = 1, col = "darkgreen", lwd = 1.5))
with(gppDF, points(gdap_amb*100~Year, type="l", lty = 2, col = "red", lwd = 1.5))
with(gppDF, points(gdap_ele*100~Year, type="l", lty = 1, col = "red", lwd = 1.5))

#legend("bottomright", c("GDAY AMB", "GDAY ELE", "GDAY-P AMB", "GDAY-P ELE"),
#       col=c("darkgreen", "darkgreen", "red", "red"), lty = c(2, 1, 2, 1))

# LAI
with(laiDF, plot(gday_amb~Year, ylim = c(1, 2), type="l", lwd = 1.5,
                 ylab = "LAI", lty = 2, col = "darkgreen"))
with(laiDF, points(gday_ele~Year, type="l", lty = 1, col = "darkgreen", lwd = 1.5))
with(laiDF, points(gdap_amb~Year, type="l", lty = 2, col = "red", lwd = 1.5))
with(laiDF, points(gdap_ele~Year, type="l", lty = 1, col = "red", lwd = 1.5))


dev.off()
