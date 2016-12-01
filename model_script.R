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
vcf <- na.omit(vcf)

# extract GLC200 land cover- note these are only cells with at least 50% of dominant LC class
glc <- as.data.frame(getValues(glc.50pct.ease))
colnames(glc) <- 'lc'
glc$cell.id <- 1:nrow(glc)
glc <- na.omit(glc)

#extract yearley swe/temp values, each layer is a column
snow <- getValues(swe.cru)
# convert to vector - data extracted by column, we checked
snow <- as.data.frame(as.vector(snow))
colnames(snow) <- 'swe'
snow <- na.omit(snow)
snow$year <- rep(1980:2014, each=nrow(glc))
snow$cell.id <- rep(1:nrow(glc),35)

#combine data frames for ensuing wizardry
covar <- join(vcf,glc,type='inner',by='cell.id')
all.dat <- join(covar,snow,type='inner',by='cell.id')

#make unique ids for year and LC

year <- data.frame(year=unique(all.dat$year),year.id=seq(1,length(unique(all.dat$year))))

all.dat <- join(all.dat,year,type='left',by='year')

lc <- data.frame(lc=unique(all.dat$lc),lc.id=seq(1,length(unique(all.dat$lc))))
all.dat <- join(all.dat,lc,type='left',by='lc')


