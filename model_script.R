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
require(plyr)
require(raster)
require(ncdf4)
require(xlsx)
require(gdalUtils)
require()

setwd("C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker")
load("C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/monthly_swe_analysis_29Nov.RData")

##extract MODIS tree cover data
vcf <- as.data.frame(getValues(vcf.2014.ease))
colnames(vcf) <- 'vcf'
vcf$cell.id <- 1:nrow(vcf)

# extract GLC200 land cover- note these are only cells with at least 50% of dominant LC class
glc <- as.data.frame(getValues(glc.50pct.ease))
colnames(glc) <- 'lc'
glc$cell.id <- 1:nrow(glc)

#extract yearley swe/temp values, each layer is a column
snow <- getValues(swe.cru)
# convert to vector - data extracted by column, we checked
snow <- as.data.frame(as.vector(snow))
colnames(snow) <- 'swe'



snow$year <- rep(1980:2014, each=nrow(glc))
snow$cell.id <- rep(1:nrow(glc),35)