##########################################
##script to look at snow dynamics ########
##########################################

install.packages(c("raster","ncdf4","gdalUtils"))

library(raster)
library(ncdf4)
library(gdalUtils)

#read in test data file

#list the monthly files

swe <- nc_open("c:\\Users\\hkropp\\Google Drive\\swe_test\\swe\\SWE_obsMEAN_M2BGS_ease2.2010.nc")
attributes(swe$dim)$names[1]


lat <- ncvar_get(swe, "y")
nlat <- dim(lat)


depletion<- function(b,day,mid,M){
	(M/(1+exp(b*(day-mid))))+2
}

plot(seq(50,100),depletion(1,seq(50,100),75,10), type="l")