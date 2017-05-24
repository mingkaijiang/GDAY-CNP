######################################################################################
# Download met data from HIEv and process them
# extraction code provided by Jim Yang



#### Incomplete!!!!!!

######################################################################################
### packages
if(!dir.exists("download"))dir.create("download")

if(!require(HIEv)){
    stop("Install the HIEv package first from bitbucket.org/remkoduursma/HIEv")
}

setToken(tokenfile="tokenfile.txt", quiet=TRUE)
setToPath("download")

if(!require(pacman))install.packages("pacman")
pacman::p_load(dplyr, doBy, readxl, lubridate) # add other packages needed to this list

######################################################################################
### download met data

### set start and end date
startDate <- as.Date("2013-01-01", format = "%Y-%m-%d")
endDate <- as.Date("2016-12-31", format = "%Y-%m-%d")


ros15 <- downloadTOA5("ROS_WS_Table15",
                      startDate = startDate,
                      endDate = endDate)
ros05 <- downloadTOA5("ROS_WS_Table05",
                      startDate = startDate,
                      endDate = endDate)


ros05_30 <- as.data.frame(dplyr::summarize(group_by(ros05,DateTime=nearestTimeStep(DateTime,30)),
                                           PPFD=mean(PPFD_Avg, na.rm=TRUE),
                                           Tair=mean(AirTC_Avg, na.rm=TRUE),
                                           RH=mean(RH, na.rm=TRUE)))
ros15_30 <- as.data.frame(dplyr::summarize(group_by(ros15,DateTime=nearestTimeStep(DateTime,30)),
                                           Rain=sum(Rain_mm_Tot, na.rm=TRUE)))

######################################################################################
# get real co2 data from hiev
fn <- sprintf("R%s_FCPLOGG_R",Ring)
Rawdata <- downloadTOA5(fn,
                        maxnfiles = 600,
                        rowbind=FALSE,
                        startDate = startDate,
                        endDate = endDate)     