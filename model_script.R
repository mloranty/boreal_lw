##########################
#
# examine snow dynamics 
# for picker project
# preliminary agu analyses
# for 2016 fall mtg
#
# MML 12/01/16
##########################

rm(list=ls())

require(raster)
require(ncdf4)
require(xlsx)
require(gdalUtils)


setwd("C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker")
load("C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/monthly_swe_analysis_29Nov.RData")