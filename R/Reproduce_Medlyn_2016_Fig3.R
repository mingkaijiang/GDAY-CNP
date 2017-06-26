#### To reproduce figures in Medlyn et al. 2016 GCB Figure 3 panel a-d

#### Read in all model outputs
#### They are downloaded according to paper supplementary materials
#### Only the fixed met data output was used

cabl <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/CABL/outputs/D1CABLEUCAMBAVG.csv",
                 skip=7)
clm4 <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/CLM4/outputs/D1CLM4EUCAMBAVG.csv",
                 skip=0)
clmp <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/CLMP/outputs/D1CLMPEUCAMBAVG.csv",
                 skip=0)
gday <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/GDAY/outputs/D1GDAYEUCAMBAVG.csv",
                 skip=4)
lpjg <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/LPJG/outputs/D1LPJXEUCAMBAVG.csv",
                 skip=2)
ocnx <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/OCNX/outputs/D1OCNXEUCAMBAVG.csv",
                 skip=0)
sdvm <- read.csv("~/Documents/Research/Projects/eucface/Cross_model_comparison/Medlyn_raw_data/mdekauwe-EucFACE_paper-50fb228/models/SDVM/outputs/D1SDVMEUCAMBAVG.csv",
                 skip=0)