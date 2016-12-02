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
require(rjags)
require(xtable)

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
snow$year <- rep(1980:2014, each=nrow(snow)/35)
snow$cell.id <- rep(1:(nrow(snow)/35),35)
snow <- na.omit(snow)

#combine data frames for ensuing wizardry
covar <- join(vcf,glc,type='inner',by='cell.id')
all.dat <- join(covar,snow,type='inner',by='cell.id')

#make unique ids for year and LC

year <- data.frame(year=unique(all.dat$year),year.id=seq(1,length(unique(all.dat$year))))

all.dat <- join(all.dat,year,type='left',by='year')

lc <- data.frame(lc=unique(all.dat$lc),lc.id=seq(1,length(unique(all.dat$lc))))
all.dat <- join(all.dat,lc,type='left',by='lc')


## perpare data for the model

data.model.list <- list(Nobs=dim(all.dat)[1],SWE=all.dat$swe,Tree.cov=all.dat$vcf,LandID=all.dat$lc.id,NyearS=dim(year)[1],
                        yearID=all.dat$year.id,Nland=dim(lc)[1],Nyear=dim(year)[1],xS=rep(1,dim(year)[1]),yS=year$year)

sample.list <- c("Beta1star","Beta2","eps.star","sig.eps","rho.eps","sig.SWE")

inits <- list(list(t.eps=1,rho.eps=0.9),list(t.eps=1.5,rho.eps=0.8),list(t.eps=0.5,rho.eps=0.7))

##
model.init <- jags.model(file="C:\\Users\\mloranty\\Documents\\GitHub\\boreal_lw\\model_code.r",
                         data=data.model.list,n.adapt=1000,n.chains=3,inits=inits)

n.iter <- 60000
n.thin <- 20

coda.obj <- coda.samples(model.init,variable.names = sample.list,n.iter=n.iter,thin=n.thin)

##
plot(coda.obj,ask=TRUE)

