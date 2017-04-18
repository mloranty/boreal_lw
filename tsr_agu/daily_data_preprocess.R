##########################
#
# daily temp/snow dynamics
# for picker project
# reanalysis pre-processing
# 
# MML 02/28/17
##########################

rm(list=ls())
require(caTools)
require(raster)
require(ncdf4)
require(xlsx)
require(gdalUtils)

## link to notes on era time-steps - 
## https://software.ecmwf.int/wiki/pages/viewpage.action?pageId=56658233 
## https://software.ecmwf.int/wiki/display/CKB/How+to+get+daily+temperature+max,+min+and+mean+from+ERA-Interim


## get example file with ease2 projection
ease2 <- raster("C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/nhtsd25e2_20121220_v01r01.nc")
projection(ease2) <- "+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m"
ease <- glc.50 <- raster("C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_50_mask.tif")

####################################################
## see monthly_snow_temp.R for veg map processing ##
####################################################
## read in and reproject veg data ##
glc.50 <- raster("C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_50_mask.tif")
glc.50 <- projectRaster(glc.50,ease2,method='ngb',progress='TEXT',overwrite=T,
                        filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE2_from_0.25_50_mask.tif')

glc.75 <- raster("C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_75_mask.tif")
glc.75 <- projectRaster(glc.75,ease2,method='ngb',progress='TEXT',overwrite=T,
                        filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE2_from_0.25_75_mask.tif')

vcf <- raster("C:/Users/mloranty/Google Drive/GIS_Data/MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_ease.tif")
vcf <- projectRaster(vcf,ease2,progress='text',overwrite=T,  
                     filename='C:/Users/mloranty/Google Drive/GIS_Data/MODIS/MOD44B/hdf_tiles/MOD44B_2014_ease2.tif')
###################################################
### load daily SWE files from Lawrence Mudryk    ##
### these are a combo of GlobSnow, MERRA & Brown ##
###################################################
setwd("C:/Users/mloranty/Google Drive/GIS_Data/swe_mudryk/")
swe <- brick('SWE_obsMEAN_M2BGS_ease2.2010.nc')
extent(swe) <- extent(ease2)

## note this is a modified PROJ.4 description for EASE2, but not sure this works with Lawrence's files ##
projection(swe) <- "+proj=laea +lat_0=90 +lon_0=-90 +ellps=WGS84 +datum=WGS84 +units=m"

swe <- subset(swe, which(substr(swe@z$Date,6,7) =="03" | substr(swe@z$Date,6,7) =="04" | substr(swe@z$Date,6,7) =="05"))
writeRaster(swe,filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/SWE_obsMEAN_M2BGS_ease2_MAM_2010.tif', overwrite=T)

swe <- projectRaster(swe,ease2, overwrite=T,
                     filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/SWE_obsMEAN_M2BGS_ease2_MAM_2010.tif')


####################################################
##      ERA Interim Reanalysis Data               ##
## see notes and links for processing             ##
## Temp is Analsyses data - daily means           ##
## radiation vars or 12hr forecasts - daily sums  ##
####################################################

## note that the ERA radiation data are in a 12hr timestep and units are joules/m2 - so divide by 12hrs for W/m2
## http://apps.ecmwf.int/codes/grib/param-db?id=169
## http://apps.ecmwf.int/codes/grib/param-db?id=176

########################################## Downwelling Longwave Radiation ##########################################

era.lw <- stack('C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_2014_Jan_May_strd_0.125deg.nc')

#subset  MAM 2010 for now #
era.lw <-  subset(era.lw,which(substr(names(era.lw),2,8)=="2010.03" 
							| substr(names(era.lw),2,8)=="2010.04"
							| substr(names(era.lw),2,8)=="2010.05"))

## data are 12hr forecasts beginning at 0 and 1200, so need to combine to get daily flux
d <- unique(substr(names(era.lw),2,11))
era.lw.d <- subset(era.lw, which(substr(names(era.lw),2,11)==d[1]))
era.lw.d <- calc(era.lw.d, fun=sum)

for(i in 2:length(d))
{
  n <- subset(era.lw, which(substr(names(era.lw),2,11)==d[i]))
  n <- calc(n, fun=sum)
  era.lw.d <- stack(era.lw.d,n)
}
rm(n)  
names(era.lw.d) <- d

era.lw.d <- crop(era.lw.d,c(0,360,45,90),overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_Mar_Apr_LW_daily_NH_WGS84.tif')
#convert from J/s m2 to W/m2
era.lw.d <- era.lw.d/(24*60*60)

era.lw.d <- rotate(era.lw.d,progress='TEXT',overwrite=T,
              filename="C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_MAM_LW_daily_NH_0.125_WGS84.tif")
era.lw.d.w <- crop(era.lw.d,c(-180,0,45,90),overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_MAM_LW_daily_NH_West_0.125_WGS84.tif')
era.lw.d.e <- crop(era.lw.d,c(-180,0,45,90),overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_MAM_LW_daily_NH_East_0.125_WGS84.tif')
			
############################################# Downwelling Shortwave Radiation ##########################################

era.sw <- stack('C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010__MAM_ssrd_boreal_0.125deg.nc')

era.sw <-  subset(era.sw,which(substr(names(era.sw),2,8)=="2010.03" 
							| substr(names(era.sw),2,8)=="2010.04"
							| substr(names(era.sw),2,8)=="2010.05"))

## data are 12hr forecasts beginning at 0 and 1200, so need to combine to get daily flux
d <- unique(substr(names(era.sw),2,11))
era.sw.d <- subset(era.sw, which(substr(names(era.sw),2,11)==d[1]))
era.sw.d <- calc(era.sw.d, fun=sum)

for(i in 2:length(d))
{
  n <- subset(era.sw, which(substr(names(era.sw),2,11)==d[i]))
  n <- calc(n, fun=sum)
  era.sw.d <- stack(era.sw.d,n)
}
rm(n)
names(era.sw.d) <- d

#era.sw.d <- crop(era.sw.d,c(0,360,45,90),overwrite=T,
#                 filename='C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_Mar_Apr_SW_daily_NH_0.125_WGS84.tif')
#convert from J/s m2 to W/m2
era.sw.d <- era.sw.d/(24*60*60)
writeRaster(era.sw.d,overwrite=T,
                 filename='C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_MAM_SW_daily_NH_0.125_WGS84.tif')

# rotate not needed when file are clipped before donwload from EWMCF
#era.sw.d <- rotate(era.sw.d,progress='TEXT',overwrite=T,
 #                  filename="C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_MAM_SW_daily_NH_0.125_WGS84.tif")
era.sw.d.w <- crop(era.sw.d,c(-180,0,45,90),overwrite=T,snap='near',
                 filename='C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_MAM_SW_daily_NH_West_0.125_WGS84.tif')
era.sw.d.e <- crop(era.sw.d,c(0,180,45,90),overwrite=T,snap='near',
                 filename='C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_MAM_SW_daily_NH_East_0.125_WGS84.tif')				 
########################################## 2m Air Temperature ##########################################

era.t <- stack('C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_MAM_6hr_2m_temp_0.125deg.nc')

era.t <-  subset(era.t,which(substr(names(era.t),2,8)=="2010.03" 
							| substr(names(era.t),2,8)=="2010.04"
							| substr(names(era.t),2,8)=="2010.05"))

## data are 3hr analyses, so need to calculate daily mean
d <- unique(substr(names(era.t),2,11))
era.t.d <- subset(era.t, which(substr(names(era.t),2,11)==d[1]))
era.t.d <- calc(era.t.d, fun=mean)

for(i in 2:length(d))
{
  n <- subset(era.t, which(substr(names(era.t),2,11)==d[i]))
  n <- calc(n, fun=mean)
  era.t.d <- stack(era.t.d,n)
}
rm(n)
names(era.t.d) <- d

era.t.d <- crop(era.t.d,c(0,360,45,90),overwrite=T,
                 filename='C:/Users/mloranty/Google Drive/GIS_Data/era_interim/era_interim_2010_MAM_2mTemp_daily_NH_0.125_WGS84.tif')


era.t.d <- rotate(era.t.d,progress='TEXT',overwrite=T,
                   filename="C:\\Users\\mloranty\\Google Drive\\GIS_Data\\era_interim\\era_interim_2010_MAM_2mTemp_daily_NH_0.125_WGS84.tif")
era.t.d.w <- crop(era.t.d, c(-180,0,45,90),overwrite=T,
					filename="C:\\Users\\mloranty\\Google Drive\\GIS_Data\\era_interim\\era_interim_2010_MAM_2mTemp_daily_NH_West_0.125_WGS84.tif")
era.t.d.e <- crop(era.t.d, c(0,180,45,90),overwrite=T,
					filename="C:\\Users\\mloranty\\Google Drive\\GIS_Data\\era_interim\\era_interim_2010_MAM_2mTemp_daily_NH_East_0.125_WGS84.tif")	
					
					
########################################## Reproject for Analysis ##########################################
## writing all data to this project folder now ##
setwd('C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data')

era.t.d <- projectRaster(era.t.d, ease2,filename="era_interim_2010_MAM_2mTemp_daily_NH_EASE2.tif",overwrite=T) 
era.t.d.w <- projectRaster(era.t.d.w, ease2,filename="era_interim_2010_MAM_2mTemp_daily_NH_west_EASE2.tif",overwrite=T)
era.t.d.e <- projectRaster(era.t.d.e, ease2,filename="era_interim_2010_MAM_2mTemp_daily_NH_east_EASE2.tif",overwrite=T)
era.sw.d <- projectRaster(era.sw.d, ease2,filename="era_interim_2010_MAM_SW_daily_NH_EASE2.tif",overwrite=T) 
era.sw.d.w <- projectRaster(era.sw.d.w, ease2,filename="era_interim_2010_MAM_SW_daily_NH_west_EASE2.tif",overwrite=T)
era.sw.d.e <- projectRaster(era.sw.d.e, ease2,filename="era_interim_2010_MAM_SW_daily_NH_east_EASE2.tif",overwrite=T)
era.lw.d <- projectRaster(era.lw.d, ease2,filename="era_interim_2010_MAM_LW_daily_NH_EASE2.tif",overwrite=T) 
era.lw.d.w <- projectRaster(era.lw.d.w, ease2,filename="era_interim_2010_MAM_LW_daily_NH_west_EASE2.tif",overwrite=T)
era.lw.d.e <- projectRaster(era.lw.d.e, ease2,filename="era_interim_2010_MAM_LW_daily_NH_east_EASE2.tif",overwrite=T)

#####################################################################################################################
#######################         SEPARATE VEG DATA AND CREATE BINS OF TREE COVER FOR ANALYSIS	#####################
#####################################################################################################################
setwd('C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data')

glc.50 <- raster('GLC2000_EASE2_from_0.25_50_mask.tif')				   			   
glc.75 <- raster('GLC2000_EASE2_from_0.25_75_mask.tif')
vcf <- raster('MOD44B_2014_ease2.tif')	

#####################################################################################################################
## make multiband rasters with each GLC2000 forest class in separate layer 
#deciduous broadleaf
lc <- cbind(unique(glc.50), NA)
lc[1:2,2] <- 3
glc <- reclassify(glc.50, lc)
glc75 <- reclassify(glc.75, lc)
#deciduous needleleaf
lc <- cbind(unique(glc.50), NA)
lc[3,2] <- 4
glc <- stack(glc,reclassify(glc.50, lc))
glc75 <- stack(glc75,reclassify(glc.75, lc))
#evergreen needleleaf
lc <- cbind(unique(glc.50), NA)
lc[4,2] <- 5
glc <- stack(glc,reclassify(glc.50, lc))
glc75 <- stack(glc75,reclassify(glc.75, lc))
#mixed leaf
lc <- cbind(unique(glc.50), NA)
lc[5,2] <- 6
glc <- stack(glc,reclassify(glc.50, lc))
glc75 <- stack(glc75,reclassify(glc.75, lc))
#mosaic forest
lc <- cbind(unique(glc.50), NA)
lc[6,2] <- 9
glc <- stack(glc,reclassify(glc.50, lc))
glc75 <- stack(glc75,reclassify(glc.75, lc))

writeRaster(glc, filename='GLC2000_EASE2_from_0.25_50_mask_forest_layers.tif',overwrite=T)
writeRaster(glc75,filename='GLC2000_EASE2_from_0.25_75_mask_forest_layers.tif',overwrite=T)

## make multiband rasters with each MODIS VCF for each GLC2000 forest class 
vcf50 <- mask(vcf, glc,filename="MOD44B_2014_ease2_GLC2000_50mask_forest_layers.tif",overwrite=T)
vcf75 <- mask(vcf, glc75,filename="MOD44B_2014_ease2_GLC2000_75mask_forest_layers.tif",overwrite=T)

## make multiband rasters with bins of MODIS VCF for each GLC2000 forest class 
bin <- matrix(c(0,20,40,60,80,20,40,60,80,100,10,30,50,70,90),ncol=3)
vcf50bin <- reclassify(vcf50,bin,filename="MOD44B_2014_ease2_GLC2000_50mask_forest_layers_vcf_bins.tif",overwrite=T)
vcf75bin <- reclassify(vcf75,bin,filename="MOD44B_2014_ease2_GLC2000_75mask_forest_layers_vcf_bins.tif",overwrite=T)

#####################################################################################################################
#####################################################################################################################
#######################         READ IN CLIPPED AND REPROJECTED PROJECT DATA	#####################################
##########			SEE ABOVE FOR INFORMATION ON PRE-PROCESSING AND ORIGINAL DATA FILES				     ############
#####################################################################################################################
#####################################################################################################################
setwd('C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data')

glc <- stack('GLC2000_EASE2_from_0.25_50_mask_forest_layers.tif')
glc75 <- stack('GLC2000_EASE2_from_0.25_50_mask_forest_layers.tif')
vcf <- stack("MOD44B_2014_ease2_GLC2000_50mask_forest_layers.tif")
vcf75 <- stack("MOD44B_2014_ease2_GLC2000_75mask_forest_layers.tif")
vcf.bin <- stack("MOD44B_2014_ease2_GLC2000_50mask_forest_layers_vcf_bins.tif")
vcf75bin <- stack("MOD44B_2014_ease2_GLC2000_75mask_forest_layers_vcf_bins.tif")

swe <- stack('SWE_obsMEAN_M2BGS_ease2_MAM_2010.tif')
era.t <- stack("era_interim_2010_MAM_2mTemp_daily_NH_EASE2.tif")
#era.tw <- stack("era_interim_2010_MAM_2mTemp_daily_NH_west_EASE2.tif")
#era.te <- stack("era_interim_2010_MAM_2mTemp_daily_NH_east_EASE2.tif")
era.sw <- stack("era_interim_2010_MAM_SW_daily_NH_EASE2.tif") 
#era.sw.w <- stack("era_interim_2010_MAM_SW_daily_NH_west_EASE2.tif")
#era.sw.e <- stack("era_interim_2010_MAM_SW_daily_NH_east_EASE2.tif")
era.lw <- stack("era_interim_2010_MAM_LW_daily_NH_EASE2.tif") 
#era.lw.w <- stack("era_interim_2010_MAM_LW_daily_NH_west_EASE2.tif")
#era.lw.e <- stack("era_interim_2010_MAM_LW_daily_NH_east_EASE2.tif")

######### claculate TSR and SSR ############
swe.dif <- swe[[1:91]]-swe[[2:92]]
writeRaster(swe.dif,filename='swedif_SWE_obsMEAN_M2BGS_ease2_MAM_2010.tif')
swe.dif <- reclassify(swe.dif,matrix(c(-Inf,0,NA),ncol=3),overwrite=T,
					filename='swedif_rcl_SWE_obsMEAN_M2BGS_ease2_MAM_2010.tif')


sw.dif <- era.sw[[2:92]]-era.sw[[1:91]]
writeRaster(t.dif,filename='swdif_era_interim_2010_MAM_2mTemp_daily_NH_EASE2.tif')

era.t.rcl <- reclassify(era.t,matrix(c(-Inf,273.15,NA),ncol=3),overwrite=T,
						filename='rcl_era_interim_2010_MAM_2mTemp_daily_NH_EASE2.tif')
t.dif <- era.t.rcl[[2:92]]-era.t.rcl[[1:91]]
writeRaster(t.dif,filename='tdif_era_interim_2010_MAM_2mTemp_daily_NH_EASE2.tif')
t.dif <- reclassify(t.dif,matrix(c(-Inf,0,NA),ncol=3),overwrite=T,
					filename='tdif_rcl_era_interim_2010_MAM_2mTemp_daily_NH_EASE2.tif')

tsr <- swe.dif/t.dif
tsr <- reclassify(tsr,rbind(c(-Inf,0,NA),c(20,Inf,NA)),overwrite=T,
                  filename='tsr_era_interim_SWE_M2BGS.tif')
tsr.mean <- calc(tsr,fun=mean, na.rm=T)
writeRaster(tsr.mean,filename='ssr_era_interim_SWE_M2BGS.tif',overwrite=T)

ssr <- swe.dif/sw.dif
writeRaster(ssr,filename='ssr_era_interim_SWE_M2BGS.tif',overwrite=T)






##################################################
## summarize daily met data by veg type/VCF bin ##
##################################################
dbf <- list(zonal(era.t,vcf.bin[[1]],fun=mean))
dbf[[2]] <- zonal(swe,vcf.bin[[1]],fun=mean)
dbf[[3]] <- zonal(era.sw,vcf.bin[[1]],fun=mean)
dbf[[4]] <- zonal(era.lw,vcf.bin[[1]],fun=mean)

enf <- list(zonal(era.t,vcf.bin[[2]],fun=mean))
enf[[2]] <- zonal(swe,vcf.bin[[2]],fun=mean)
enf[[3]] <- zonal(era.sw,vcf.bin[[2]],fun=mean)
enf[[4]] <- zonal(era.lw,vcf.bin[[2]],fun=mean)

dnf <- list(zonal(era.t,vcf.bin[[3]],fun=mean))
dnf[[2]] <- zonal(era.t,vcf.bin[[3]],fun=sd)
dnf[[3]] <- zonal(swe,vcf.bin[[3]],fun=mean)
dnf[[4]] <- zonal(swe,vcf.bin[[3]],fun=sd)
dnf[[5]] <- zonal(era.sw,vcf.bin[[3]],fun=mean)
dnf[[6]] <- zonal(era.sw,vcf.bin[[3]],fun=sd)
dnf[[7]] <- zonal(era.lw,vcf.bin[[3]],fun=mean)	
dnf[[8]] <- zonal(era.lw,vcf.bin[[3]],fun=sd)	
			   
##############################################				   
## plot MAM timeseries of daily MET drivers ##
##############################################

pdf('dnf_MAM_met.pdf',5,11)		
par(mfcol=c(4,1))		   
plot(61:152,dnf[[1]][1,2:93],type='l',
		ylim=c(0,300),xlab="DOY",ylab="Temp/SW/LW/SWE")				   
lines(61:152,dnf[[3]][1,2:93],type='l',col="blue")				   
lines(61:152,dnf[[5]][1,2:93],type='l',col='orange')					   
lines(61:152,dnf[[7]][1,2:93],type='l',col='gray')	
abline(h=273.15,lty='dotted')
			   
plot(61:152,dnf[[1]][2,2:93],type='l',
		ylim=c(0,300),xlab="DOY",ylab="Temp/SW/LW/SWE")				   
lines(61:152,dnf[[3]][2,2:93],type='l',col="blue")				   
lines(61:152,dnf[[5]][2,2:93],type='l',col='orange')					   
lines(61:152,dnf[[7]][2,2:93],type='l',col='gray')
abline(h=273.15,lty='dotted')
				   
plot(61:152,dnf[[1]][3,2:93],type='l',ylim=c(0,300),ylab="Temp/SW/LW/SWE",xlab="DOY")				   
lines(61:152,dnf[[3]][3,2:93],type='l',col="blue")				   
lines(61:152,dnf[[5]][3,2:93],type='l',col='orange')					   
lines(61:152,dnf[[7]][3,2:93],type='l',col='gray')					   
abline(h=273.15,lty='dotted')
	
plot(61:152,dnf[[1]][4,2:93],type='l',ylim=c(0,300),ylab="Temp/SW/LW/SWE",xlab="DOY")				   
lines(61:152,dnf[[3]][4,2:93],type='l',col="blue")				   
lines(61:152,dnf[[5]][4,2:93],type='l',col='orange')					   
lines(61:152,dnf[[7]][4,2:93],type='l',col='gray')
abline(h=273.15,lty='dotted')
				   
legend('bottomleft',c('Tair','SWE','SW','LW'),ncol=4,cex=0.8,
		fill=c('black', 'blue', 'orange','gray'),bty='n')			   
dev.off()				   
				   
				   
##############################################
# get some TSR values
dnf.tsr <- cbind(getValues(vcf[[3]]),getValues(tsr[[1]]))
enf.tsr <- cbind(getValues(vcf[[2]]),getValues(tsr[[1]]))
for(i in 2:nlayers(tsr))
{
 enf.tsr <- rbind(enf.tsr,cbind(getValues(vcf[[3]]),getValues(tsr[[i]])))
 dnf.tsr <- rbind(dnf.tsr,cbind(getValues(vcf[[3]]),getValues(tsr[[i]])))
}				   
				   
plot(dnf.tsr[,1],dnf.tsr[,2])
				   

pdf('TSR_daily_mean_MAM_2010.pdf',5,11)		
par(mfcol=c(4,1))				   			   
plot(getValues(vcf[[1]]),getValues(tsr.mean),xlab='VCF',ylab='TSR',col='lightgreen')
legend('topright','DBF',bty='n')				   
plot(getValues(vcf[[2]]),getValues(tsr.mean),xlab='VCF',ylab='TSR',col='darkgreen')
legend('topright','ENF',bty='n')				   
plot(getValues(vcf[[3]]),getValues(tsr.mean),xlab='VCF',ylab='TSR',col='orange')
legend('topright','DNF',bty='n')				   
plot(getValues(vcf[[4]]),getValues(tsr.mean),xlab='VCF',ylab='TSR',col='maroon')
legend('topright','MIX',bty='n')					   
dev.off()				   
				   
				   
				   
				   