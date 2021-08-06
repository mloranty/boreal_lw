###########################################################
###########################################################
###########################################################
############ SWE stats model v2.0              ############  
############ modified from v1.0                ############
############ H. Kropp                          ############
###########################################################
###########################################################
############ Data outputs:                     ############
############  A: Raster stacks:                ############
############    1: dailySwe                    ############
############      yearly list stack daily swe  ############
############    2: dailySwe.mask               ############
############      masked daily swe for analysis############
############    3: melt.mm.day                 ############
############      annual spring melt rate      ############
############    4: meltDuration                ############
############      annual duration of melt      ############
############    5: doyStart                    ############
############      annual start of melt         ############
############    6: glc2000                     ############
############      land cover of interest       ############
############    7: maxSwe                      ############
############      annual maximum               ############
############      Swe at start melt            ############
############    8: vcf.mask                    ############
############      % tree cover                 ############
############  B: Dataframes:                   ############
############    1: glcID                       ############
############      table of glcID codes         ############
############    2: analysisDF                  ############
############      table with no missing data   ############
############      table with no missing data   ############
###########################################################
###########################################################



###########################################
########## Source data        -----
source("c:/Users/hkropp/Documents/GitHub/boreal_lw/HK_swe_organize.R")

###########################################
########## Libraries       -----
library(sf)
library(maps)
library(rgeos)
library(raster)

###########################################
########## Directories       -----
outDir <- "E:/Google Drive/research/projects/boreal_swe/boreal_2021/figures"


###########################################
########## Additional data      -----


#regression info


#read in model output
datS <- read.csv(paste0(modDir,"\\curve_mod_stats.csv"))
datQ <- read.csv(paste0(modDir,"\\curve_mod_quant.csv"))

#combine data frames
datC <- cbind(datS,datQ)
#pull out parameter names
dexps<-"\\[*[[:digit:]]*\\]"
datC$parm <- gsub(dexps,"",rownames(datC))
unique(datC$parm)
datC$parm2 <- gsub("\\d","",rownames(datC))

#pull out parameters
#transformed intercepts
beta0NL <- datC[datC$parm == "trB0",] 
#nontransformed regression parameters
beta0 <- datC[datC$parm == "betaB0S",] 
beta1 <- datC[datC$parm == "betaB1",] 
beta2 <- datC[datC$parm == "betaB2",] 
beta3 <- datC[datC$parm == "betaB3",] 
beta4 <- datC[datC$parm == "betaB4",] 
#add indicator if parameter is significant
#temp during melt period
beta1$sig <- ifelse(beta1$X2.5.<0&beta1$X97.5.<0,1,
                    ifelse(beta1$X2.5.>0&beta1$X97.5.>0,1,0))
#Canopy cover - 20%						
beta2$sig <- ifelse(beta2$X2.5.<0&beta2$X97.5.<0,1,
                    ifelse(beta2$X2.5.>0&beta2$X97.5.>0,1,0))						
#onset doy - 107 doy (middle of time period)				
beta3$sig <- ifelse(beta3$X2.5.<0&beta3$X97.5.<0,1,
                    ifelse(beta3$X2.5.>0&beta3$X97.5.>0,1,0))
#maximum swe value -0.15m
beta4$sig <- ifelse(beta4$X2.5.<0&beta4$X97.5.<0,1,
                    ifelse(beta4$X2.5.>0&beta4$X97.5.>0,1,0))	

#check if any nonsignificant slopes to account for

length(which(beta1$sig == 0))						

length(which(beta3$sig == 0))	


#pull out regression means for plotting

mu.Temp <- datC[datC$parm2 == "mu.Temp[,]",]
mu.Onset <- datC[datC$parm2 == "mu.Onset[,]",]
mu.Max <- datC[datC$parm2 == "mu.Max[,]",]
mu.Canopy <- datC[datC$parm2 == "mu.Canopy[,]",]	


###########################################
########## Map boundaries      -----


#make a circle for the polar region
#max missing is predominately zero throughout the entire map. Not worth showing
worldmap <- map("world", ylim=c(40,90), fill=TRUE)
#focus on a smaller extent
worldmap2 <- map("world", ylim=c(50,90))

#world map
world <- project(matrix(c(worldmap$x,worldmap$y), ncol=2,byrow=FALSE),laea)
world2 <- project(matrix(c(worldmap2$x,worldmap2$y), ncol=2,byrow=FALSE),laea)

###make polygon to cover up non study area####
#make a point at center of map
pt1 <- SpatialPoints(data.frame(x=0,y=0), CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#make a buffer around center of plot choosing distance that is relevant
ptBuff <- buffer(pt1,4500000)
#set up bounds to extend beyond study area
xcor <- c(-8000000,-8000000,8000000,8000000)
ycor <- c(-8000000,8000000,8000000,-8000000)
#make empty plot polygon
boxC <- cbind(xcor,ycor)
p <- Polygon(boxC)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
proj4string(sps) = CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" )
#remove study area from empty plot
PolyBlock <- gDifference(sps,ptBuff, byid=TRUE)


###########################################
########## Color palette       -----

glcID$gcID <- seq(1,5)

vegePallete <- c(rgb(50/255,80/255,10/255), #evergreen needleleaf,
                 rgb(130/255,160/255,190/255),# deciduous needleleaf,
                 rgb(250/255,120/255,80/255), #mixed boreal,
                 rgb(60/255,60/255,110/255), #deciduous shrub
                  rgb(170/255,190/255,140/255)) #herbaceous)
vegePallete2 <-	c(rgb(50/255,80/255,10/255,.1),	
                  rgb(130/255,160/255,190/255,.1),
                  rgb(250/255,120/255,80/255,.1),
                  rgb(60/255,60/255,110/255,.1),
                  rgb(170/255,190/255,140/255,.1))		
vegePallete3 <-	c(rgb(50/255,80/255,10/255,.5),	
                  rgb(130/255,160/255,190/255,.5),
                  rgb(250/255,120/255,80/255,.5),
                  rgb(60/255,60/255,110/255,.5),
                      rgb(170/255,190/255,140/255,.5))		

treePallete <- c(rgb(229,245,224,max=255),
                 rgb(199,233,192,max=255),
                 rgb(161,217,155,max=255),
                 rgb(116,196,118,max=255),
                 rgb(65,171,93,max=255),
                 rgb(35,139,69,max=255),
                 rgb(0,109,44,max=255),
                 rgb(0,68,27,max=255))

#add better names
glcID$name2 <- c("Evergreen needleleaf", "Deciduous needleleaf","Mixed boreal","Deciduous shrub","Herbaceous")




###########################################
########## Figure 1: map of data ##########
########## inputs including %    ##########
########## tree cover & glc      ##########
########## Figure 1       -----

hd <- 12
wd1 <- 12
wd2 <- 8
water <- rgb(149/255,218/255,255/255,.3)
land <- rgb(250,230,190, max=255)
#size of panel label
mx <- 2
#line for panel label
pll <- .5
#vege type breaks
vegeBr <- c(0,1,2,3,4,5)
canopyBr <- c(0,10,20,30,40,50,60,70,80)
sweBr <-c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5)
#size of axis
cxa <- 1.75
