### Create met data for Amazon site
###
### Simply repeat year 1 climate data, except CO2
###
### Author: Mingkai Jiang
### Created on: 27-Jun-2017

#### Function

Process_Amazon_Data <- function(inName, outName) {
    #### read in file
    myDF <- read.csv(inName, 
                     skip=5, header=F)
    
    colnames(myDF) <- c("year", "doy", "tair", "rain", "tsoil", "tam", "tpm", "tmin", "tmax", "tday",
                        "vpd_am", "vpd_pm", "co2", "ndep", "nfix", "pdep", "wind", "pres", "wind_am",
                        "wind_pm", "par_am", "par_pm")
    
    yr.range <- 2023 - 1999 + 1
    
    d <- Sys.Date()
    
    myDF <- subset(myDF, doy <= 365)
    myDF <- subset(myDF, year <= 2023)
    
    for (i in c(3:12, 14:22)) {
        myDF[,i] <- rep(myDF[1:365, i], times=yr.range)
    }
    
    myDF[(nrow(myDF)+1), 3:22] <- myDF[myDF$year == 2000 & myDF$doy == 365, 3:22]
    myDF[nrow(myDF), "year"] <- 2000
    myDF[nrow(myDF), "doy"] <- 366
    
    myDF[(nrow(myDF)+1), 3:22] <- myDF[myDF$year == 2004 & myDF$doy == 365, 3:22]
    myDF[nrow(myDF), "year"] <- 2004
    myDF[nrow(myDF), "doy"] <- 366
    
    myDF[(nrow(myDF)+1), 3:22] <- myDF[myDF$year == 2008 & myDF$doy == 365, 3:22]
    myDF[nrow(myDF), "year"] <- 2008
    myDF[nrow(myDF), "doy"] <- 366
    
    myDF[(nrow(myDF)+1), 3:22] <- myDF[myDF$year == 2012 & myDF$doy == 365, 3:22]
    myDF[nrow(myDF), "year"] <- 2012
    myDF[nrow(myDF), "doy"] <- 366
    
    myDF[(nrow(myDF)+1), 3:22] <- myDF[myDF$year == 2016 & myDF$doy == 365, 3:22]
    myDF[nrow(myDF), "year"] <- 2016
    myDF[nrow(myDF), "doy"] <- 366
    
    myDF[(nrow(myDF)+1), 3:22] <- myDF[myDF$year == 2020 & myDF$doy == 365, 3:22]
    myDF[nrow(myDF), "year"] <- 2020
    myDF[nrow(myDF), "doy"] <- 366
    
    outDF <- myDF[order(myDF$year,myDF$doy),]
    
    ### rows
    row1 <- "# simplified gday met forcing transient"
    row2 <- "# Data from 1 - 500 years"
    row3 <- paste("# Created by Mingkai Jiang: ", d, sep="")
    
    row4 <- as.list(c("#--", "--", "degC", "mm/d", "degC", "degC", "degC", "degC", "degC", "degC", 
                      "kPa", "kPa", "ppm", "t/ha/d", "t/ha/d", "t/ha/d", "m/s", "kPa", "m/s", "m/s", 
                      "mj/m2/d", "mj/m2/d"))
    row5 <- as.list(as.character(c("#year", "doy", "tair", "rain", "tsoil", "tam", "tpm", "tmin", "tmax", "tday",
                                   "vpd_am", "vpd_pm", "co2", "ndep", "nfix", "pdep", "wind", "pres", "wind_am",
                                   "wind_pm", "par_am", "par_pm")))
    
    # write into folder
    write.table(row1, outName,
                col.names=F, row.names=F, sep=",", append=F, quote = F)
    
    write.table(row2, outName,
                col.names=F, row.names=F, sep=",", append=T, quote=F)
    
    write.table(row3, outName,
                col.names=F, row.names=F, sep=",", append=T, quote=F)
    
    write.table(row4, outName,
                col.names=F, row.names=F, sep=",", append=T, quote=F)
    
    write.table(row5, outName,
                col.names=F, row.names=F, sep=",", append=T, quote=F)
    
    write.table(outDF, outName,
                col.names=F, row.names=F, sep=",", append=T, quote=F)
    
}



######### Scripts
amb_obs_n <- "~/Documents/Research/Projects/Amazon/AMAZ/drought/met_data/AmaFACE_met_data_amb_obs_co2.csv"
ele_obs_n <- "~/Documents/Research/Projects/Amazon/AMAZ/drought/met_data/AmaFACE_met_data_ele_obs_co2.csv"
amb_obs_p <- "~/Documents/Research/Projects/Amazon/AMAZ/drought_p/met_data/AmaFACE_met_data_amb_obs_co2.csv"
ele_obs_p <- "~/Documents/Research/Projects/Amazon/AMAZ/drought_p/met_data/AmaFACE_met_data_ele_obs_co2.csv"

amb_out_n <- "~/Documents/Research/Projects/Amazon/AMAZ/drought/met_data/AmaFACE_met_data_amb_1999_2023.csv"
ele_out_n <- "~/Documents/Research/Projects/Amazon/AMAZ/drought/met_data/AmaFACE_met_data_ele_1999_2023.csv"
amb_out_p <- "~/Documents/Research/Projects/Amazon/AMAZ/drought_p/met_data/AmaFACE_met_data_amb_1999_2023.csv"
ele_out_p <- "~/Documents/Research/Projects/Amazon/AMAZ/drought_p/met_data/AmaFACE_met_data_ele_1999_2023.csv"

Process_Amazon_Data(amb_obs_n, amb_out_n)
Process_Amazon_Data(ele_obs_n, ele_out_n)
Process_Amazon_Data(amb_obs_p, amb_out_p)
Process_Amazon_Data(ele_obs_p, ele_out_p)



