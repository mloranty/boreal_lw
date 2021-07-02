library(raster)
require(ncdf4)
require(rgdal)
require(gdalUtils)
require(lubridate)


# define the projection - EASE2.0 equal area grid - will use 50km resolution
# https://epsg.io/6931
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
#EPSG: 6931


###########################################
# read in data


########## SWE data from Mudryk ##########
# swe data in m from Mudryk
# daily mean SWE (GS2,MERRA2,Brown,Crocus)
# list the monthly files
swe.files <- list.files(pattern =".nc",path ="E:/Google Drive/GIS/swe_mudryk_blended",full.names=TRUE)

# read one file in to use for reprojecting
pr <- raster(swe.files[1])
# crop to boreal region
pr <- crop(pr,c(-180,180,50,90))
# reproject to 50km EASE2.0 gird
pr <- projectRaster(pr,res=50000,crs=laea,progress='text')

########## ETOPO1 Elevation ##########
#https://www.ngdc.noaa.gov/mgg/global/global.html
# info on grid vs. cell registration
topo <- raster("E:/Google Drive/GIS/boreal_swe_all_data/ETOPO1_Ice_c_geotiff.tif")
topo@crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
topoN <- crop(topo, extent(-180,180,45,90))


topo.ease <- projectRaster(topoN, crs=laea)

topo.sd <- aggregate(topo.ease,fact=50000/res(topo.ease),fun=sd)


topo.mask <- reclassify(topo.sd, matrix(c(0,200,1,
                                        200.00001,3000,0),byrow=TRUE))


topo.maskR <- resample(topo.mask, pr, method="ngb")




pr.m <- mask(pr, topo.maskR,maskvalue=0)


plot(pr.m)
