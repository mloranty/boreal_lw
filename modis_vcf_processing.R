list.files(path = "https://e4ftl01.cr.usgs.gov/MOLT/MOD44B.006/2000.03.05/")

library("HelpersMG")
setwd("L:/data_repo/gis_data/MODIS/MOD44B/hdf_tiles/V6_2000_boreal/")


xy <- coordinates(pr)

sf_transform

st_as_sf(coordinates(pr), as_points = TRUE, merge = FALSE)

# install luna package1
install.packages("remotes")
remotes::install_github("rspatial/luna")