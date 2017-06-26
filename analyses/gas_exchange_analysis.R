#### Script to analyze DE's gas exchange dataset


#### prepare libraries
library(plantecophys)
library(lme4)
library(plotrix)
library(car)

#### read in csv
myDF <- read.csv("~/Documents/Research/Projects/eucface/mesophyll_conductance/data/MEANS_species_gas_exchange.csv")

#### Check for outliers
myDF$Jmax25_TP[myDF$Jmax25_TP >= 1000] <- NA

#### visional check if the data matches with figures from DE
with(myDF, plot(Vcmax25_TP~Pmolarea, ylim=c(0, 150)))
with(myDF, plot(Jmax25_TP~Pmolarea, ylim=c(0,300)))
with(myDF, plot(TPU25~Pmolarea, ylim=c(0,15)))
### all seems fine

#### Compute linear fit coefficients and check if they match with DE results
fit_vcmax_TP <- with(myDF, lm(Vcmax25_TP~Pmolarea))
fit_jmax_TP <- with(myDF, lm(Jmax25_TP~Pmolarea))
fit_vcmax <- with(myDF, lm(Vcmax25~Pmolarea))
fit_jmax <- with(myDF, lm(Jmax25~Pmolarea))

#### DE's equation Vcmax25 = 49.958 + 6.4412*Parea_mmol, r^2 = 0.605
fit_vcmax_TP

#### DE's equation Jmax25 = 79.87 + 14.2347*Parea_mmol, r^2 = 0.723
fit_jmax_TP

### note: not the same but considering se, the difference is probably acceptable. 

#### Checking statistical significance
summary(fit_vcmax_TP)
summary(fit_jmax_TP) 
summary(fit_vcmax)
summary(fit_jmax) 

#### Compare with and without TP limitation
with(myDF, plot(Vcmax25_TP~Vcmax25))
abline(a=0, b=1)

with(myDF, plot(Jmax25_TP~Jmax25))
abline(a=0, b=1)

#### Linear models
lm1 <- with(myDF, lm(Vcmax25~-1+Narea))
summary(lm1)

lm1.2 <- with(myDF, lm(Vcmax25_TP~-1+Narea))
summary(lm1.2)

lm2 <- with(myDF, lm(Vcmax25~-1+Parea))
summary(lm2)

lm2.2 <- with(myDF, lm(Vcmax25_TP~-1+Parea))
summary(lm2.2)

lm3 <- with(myDF, lm(Jmax25~-1+Narea))
summary(lm3)

lm3.2 <- with(myDF, lm(Jmax25_TP~-1+Narea))
summary(lm3.2)

lm4 <- with(myDF, lm(Jmax25~-1+Parea))
summary(lm4)

lm4.2 <- with(myDF, lm(Jmax25_TP~-1+Parea))
summary(lm4.2)



### Note: species-mean data is not helpful for computing mixed effect linear model
###       so next to read in all species data
###       1. need to exclude certain species - cross check with the species mean data to identify which
###       2. the all species data does not exclude mesophyll conductance limitation, but
###          because TP effect doesn't appear to be big, currently I am only fitting with TP embedded.

#### read in all species data
allDF <- read.csv("~/Documents/Research/Projects/eucface/mesophyll_conductance/data/Allspecies_gas_exchange.csv")

#### Exclude species not in species mean data
allDF2 <- allDF[allDF$Species %in% myDF$Species, ]

#### Convert N and P into the correct unit
#### Current in %
#### Convert into mass/area then molarea
# unit: g / area ( I think area is m2)
allDF2$Narea <- allDF2$X.N/100 * allDF2$LMA
allDF2$Parea <- allDF2$X.P/100 * allDF2$LMA

# unit: Nmolarea = mmol / area ( I think area is m2)
allDF2$Nmolarea <- allDF2$Narea/14/0.001
allDF2$Pmolarea <- allDF2$Parea/31/0.001

allDF2$N_g_m2 <- allDF2$Narea
allDF2$P_g_m2 <- allDF2$Parea

allDF2$LMA_g_m2 <- allDF2$LMA

#### Linear model for N and P only relationships
### Step 1: N only relationship
lm1 <- with(allDF2, lm(Vcmax25~N_g_m2))
summary(lm1)

lm1.2 <- lm(Vcmax25~N_g_m2+LMA_g_m2, data=allDF2)
summary(lm1.2)

lm1.3 <- lm(Vcmax25~-1+N_g_m2, data=allDF2)
summary(lm1.3)

lm2 <- with(allDF2, lm(Vcmax25~P_g_m2))
summary(lm2)

lm2.2 <- with(allDF2, lm(Vcmax25~-1+P_g_m2))
summary(lm2.2)

lm3 <- with(allDF2, lm(Jmax25~N_g_m2))
summary(lm3)

lm3.2 <- with(allDF2, lm(Jmax25~-1+N_g_m2))
summary(lm3.2)

lm4 <- with(allDF2, lm(Jmax25~P_g_m2))
summary(lm4)

lm4.2 <- with(allDF2, lm(Jmax25~-1+P_g_m2))
summary(lm4.2)

#### Mixed effect linear model step by step
### Step 1: linear model
fit1 <- with(allDF2, lm(Vcmax25~Species + Pmolarea+ Species:Pmolarea))
summary(fit1)

allDF2$Vcmax25_pred <- predict(fit1, allDF2)
with(allDF2, plot(Pmolarea, Vcmax25, pch=1, col=Species))
with(allDF2, points(Pmolarea, Vcmax25_pred, pch=19, col=Species))

### Step 2: Linear mixed model
fit2 <- lmList(Vcmax25 ~ Pmolarea | Species, data = allDF2)
coef2 <- coef(fit2)
fitsplit <- split(allDF2, allDF2$Species)
with(allDF2, plot(Pmolarea, Vcmax25, col=Species))
for (i in 1:length(fitsplit)) {
    xmin <- min(fitsplit[[i]]$Pmolarea)
    xmax <- max(fitsplit[[i]]$Pmolarea)
    
    ablineclip(coef2[i,1], coef2[i,2], x1=xmin, x2=xmax)
}

# Random intercept only
fit3 <- lmer(Vcmax25~ Nmolarea + Pmolarea + Nmolarea:Pmolarea + (1|Species), data = allDF2)
fit4 <- lmer(Vcmax25~ Nmolarea + Pmolarea + Nmolarea:Pmolarea + (Pmolarea|Species), data=allDF2)
AIC(fit3, fit4)
anova(fit3, fit4)

# output from randome intercept model
summary(fit3)

# Using Anova from car, to get p-values for the main effects
Anova(fit3)

# Standard deviation of the intercept and slope between Species, as estiamted from the ramdom effects
VarCorr(fit4)
VarCorr(fit3)
