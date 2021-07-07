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

