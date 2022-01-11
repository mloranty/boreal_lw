############################
# 
# download MODIS VCF data
# for the year 2000 and 
# pre-process
#
# MML 07/07/21
############################

### set up working environment ----
# install luna and MODIStsp
library(remotes)
#remotes::install_github("rspatial/luna")
#remotes::install_github("ropensci/MODIStsp")
library(luna)
library(MODIStsp)
library(raster)
library(ncdf4)


### calculat standard deviatin of vcf within 50km grid cells - added 01/11/22
setwd("L:/data_repo/gis_data/")

# read a topo raster for reprojecting
topo <- rast("L:/projects/boreal_swe_depletion/data/ETOPO1_Ice_c_50km_ease_mtn_mask.tif")
# define the projection - EASE2.0 equal area grid - will use 50km resolution
# https://epsg.io/6931

# read original MODIS sinusoidal mosaic
vcf <- rast("MODIS/MOD44B/MOD44B_2014_500m_mosaic_sin.tif")

# reproject to EASR grid, retaining 500m resolution
vcf2 <- project(vcf,"epsg:6931", filename = "MODIS/MOD44B/MOD44B_2014_500m_mosaic_EASE.tif")

# mask non-tree cover pixels
vcf3 <- classify(vcf2,cbind(101,Inf,NA), "MODIS/MOD44B/MOD44B_2014_500m_mosaic_EASE_mask.tif" )

# aggregate to 50km resolution
vcf4 <- aggregate(vcf3,fact = 50000/res(vcf3), fun = sd, na.rm = T)

# resample to align with other data sets 
vcf50e <- resample(vcf4,topo,method = "near", 
                   filename = "L:/projects/boreal_swe_depletion/data/MOD44B_2014_stdev_mosaic_50km_ease.tif", overwrite = T)


# note - aggregating before reprojecting leads to large messy errors near the international date line. 

#-----------------------------------------------------------------------------------#
# section to try downloading and processing newer mosaics from new collection. 
#specify download directory and set as working directory
dd <- "L:/data_repo/gis_data/MODIS/MOD44B/hdf_tiles/V6_2000_boreal/"
setwd(dd)

### download MODIS data using luna package ----
# download year 2000 MODIS vcf for boreal region
getModis(product = "MOD44B",
         start_date = "2000-01-01",
         end_date = "2000-12-31",
         aoi = c(-180,180,50,90),
         download = TRUE,
         path = dd,
         username = un, #specify username and password outside of script
         password = pw)

### mask and mosaic MODIS data using system calls to LDOPE tools ----
# see the LDOPE users manual here https://lpdaac.usgs.gov/documents/308/LDOPE_Users_Manual.pdf

# list hdf MODIS tiles
f <- list.files(path = dd, pattern = ".hdf")

# MODIS sinusoidal projection
sinus = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
# create two mosaics and write to tif 
# a single mosaic at 250m resolution is too large. 

# mosaic the western hemisphere
system(paste("mosaic_sds -sds=Percent_Tree_Cover -of=MOD44B.A2000065.west.tif",
             paste(c(f[1:37]), collapse = " "), sep = " "))


# mosaic the eastern hemisphere
system(paste("mosaic_sds -sds=Percent_Tree_Cover -of=MOD44B.A2000065.east.tif",
             paste(c(f[38:78]), collapse = " "), sep = " "))

# maybe try MODIStsp
MODIStsp(
  gui = FALSE,
  out_folder = dd,
  selprod = "Veg_Cont_Fields_Yearly_250m (MOD44B)",
  user = un, 
  password = pw,
  start_date = "2000.01.01", 
  end_date = "2000.12.31",
  spatmeth = "bbox", 
  
  
)


### work with data from NASA HEG tools ----

vcf <- raster("MOD44B_50km_lambert_GRID_v2.tif")
vcf.rcl <- reclassify(vcf,matrix(c(101,255,NA), ncol = 3))

