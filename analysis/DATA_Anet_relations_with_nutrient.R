###########################################################################
#### Script to analyze data from EucFACE
#### Purpose: finding relationship between Jmax, Vcmax and Leaf N and Leaf P
#### Created by: Mingkai Jiang

###########################################################################
#### read in csv
myDF <- read.csv("data/FACE_P0020_R0_FERTILISATION-GASEXCHANGE-FITS_20120416 - 20140409_L2.csv")

rlt <- lm(Jmax~Leaf.N, data=myDF)
summary(rlt)

rlt <- lm(Jmax~Leaf.P, data=myDF)
summary(rlt)

rlt <- lm(Vcmax~Leaf.N, data=myDF)
summary(rlt)

rlt <- lm(Vcmax~Leaf.P, data=myDF)
summary(rlt)

with(myDF, plot(Jmax25~Leaf.N))
with(myDF, plot(Jmax25~Leaf.P))


library(MASS)
myDF2 <- na.omit(myDF)
fit1 <- lm(Jmax25~Leaf.N+Leaf.Pi+LMA+Tleaf+Amax+Leaf.P,data=myDF2)
step1 <- stepAIC(fit1, direction="both")
step1$anova # display results

fit2 <- lm(Vcmax25~Leaf.N+Leaf.Pi+LMA+Tleaf+Amax+Leaf.P,data=myDF2)
step2 <- stepAIC(fit2, direction="both")
step2$anova # display results


fit3 <- lm(Rd25~Leaf.N+Leaf.Pi+LMA+Tleaf+Amax+Leaf.P,data=myDF2)
step3 <- stepAIC(fit3, direction="both")
step3$anova # display results

fit4 <- lm(Anet~Leaf.N+Leaf.Pi+LMA+Tleaf+Amax+Leaf.P,data=myDF2)
step4 <- stepAIC(fit4, direction="both")
step4$anova # display results

fit5 <- lm(Amax~Leaf.N+Leaf.Pi+LMA+Tleaf+Leaf.P,data=myDF2)
step5 <- stepAIC(fit5, direction="both")
step5$anova # display results

library(visreg)
par(mfrow=c(2,2))
visreg(fit1, "Leaf.N")
visreg(fit1, "Leaf.P")
visreg(fit1, "LMA")
visreg(fit1, "Tleaf")
visreg(fit1, "Amax")
visreg(fit1,c("Leaf.N","Leaf.P","LMA"))
visreg2d(fit1, x="Leaf.N", y="Leaf.P")
visreg2d(fit2, x="Leaf.N", y="Leaf.P")
visreg2d(fit4, x="Leaf.N", y="Leaf.P")
visreg2d(fit5, x="Leaf.N", y="Leaf.P")

# convert
myDF2$leafNC <- myDF2$Leaf.N / 100.0 * 2
myDF2$leafPC <- myDF2$Leaf.P / 100.0 * 2
rlt <- lm(Amax~leafNC+leafPC, data=myDF2)
summary(rlt)
