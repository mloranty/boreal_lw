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
############    9: meltMeanT                   ############
############     ave temp melt period (c)      ############
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
library(BAMMtools)

###########################################
########## Directories       -----
plotDI <- "E:/Google Drive/research/projects/boreal_swe/boreal_2021/figures"
modDir <- "E:/Google Drive/research/projects/boreal_swe/boreal_2021/model/run3"

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
mu.Max[400:500,]
#log transform
analysisDF$log.melt <- log(analysisDF$abs.melt)
#log transform max swe
analysisDF$log.max <- log(analysisDF$maxSwe.m)
#regression mean plot
tempMean <- seq(floor(range(analysisDF$meltTempC)[1]),ceiling(range(analysisDF$meltTempC)[2]), length.out=200)
CanopyMean <- seq(floor(range(analysisDF$vcf)[1]),ceiling(range(analysisDF$vcf)[2]), length.out=200)
SdayMean <- seq(floor(range(analysisDF$doyStart)[1]),ceiling(range(analysisDF$doyStart)[2]), length.out=200)
MaxMean <- seq(floor(range(analysisDF$log.max)[1]*10)/10,ceiling(range(analysisDF$log.max)[2]*10)/10, length.out=200)

###########################################
########## Map boundaries      -----
laea <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

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


swePallete <- rev(c(rgb(178,24,43,max=255),
                    rgb(214,96,77,max=255),
                    rgb(244,165,130,max=255),
                    rgb(253,219,199,max=255),
                    rgb(209,229,240,max=255),
                    rgb(146,197,222,max=255),
                    rgb(67,147,195,max=255),
                    rgb(33,102,172,max=255)))
#alice blue: rgb(149/255,218/255,255/255)
#option 2: "#c0d6df"
#option 3: "#d6e2e9"
#option 4: water <- "#f0f8ff"
#option 5: water <- "#80C5DE85"
water <- "#80C5DE20"

#option 1: land <-"#f1e7dd"
#option 2: land <-"#d9d9d9"
#option 4: land <- "#fff7f0"
land <-"#ccc5b9"







###########################################
########## Figure 1: map of data ##########
########## inputs including %    ##########
########## tree cover & glc      ##########
########## Figure 1       -----

hd <- 12
wd1 <- 12
wd2 <- 8

#size of panel label
mx <- 2
#line for panel label
pll <- .5
#vege type breaks
vegeBr <- c(0,4.5,5.5,6.5,12.5,13.5)
vegeBri <- seq(0,5)
canopyBr <- c(0,10,20,30,40,50,60,70,80)

#size of axis
cxa <- 1.75


plot(vcf.mask)
plot(glc2000)


png(paste0(plotDI,"\\figure1_data_maps.png"), width = 17, height = 7, units = "in", res=300)
layout(matrix(seq(1,4),ncol=4), width=c(lcm(wd1),lcm(wd2),lcm(wd1),lcm(wd2)),height=lcm(hd))
#set up empty plot
### plot 1 vegetation type ###
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
#plot points
image(glc2000,breaks=vegeBr,col=vegePallete,add=TRUE )
mtext("A",at=4100000,side=2,line=pll, las=2,cex=mx)
plot(PolyBlock, col="white",border="white", add=TRUE)
### plot 1 legend ###
par(mai=c(0,.25,0,2))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
for(i in 1:(length(vegeBr)-1)){
  polygon(c(0,0,1,1), 
          c(vegeBri[i]/(length(vegeBri)-1),vegeBri[i+1]/(length(vegeBri)-1),vegeBri[i+1]/(length(vegeBri)-1),vegeBri[i]/(length(vegeBri)-1)),
          col=vegePallete[i],border=NA)
}
axis(4,(vegeBri[1:5]/5)+.1,glcID$name2,cex.axis=cxa,las=2)
mtext("Vegetation cover", side=3, line=1, cex=2)

### plot 2 canopy cover ###
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
#plot points
image(vcf.mask, breaks=canopyBr, col=treePallete, add=TRUE)
mtext("B",at=4100000,side=2,line=pll, las=2,cex=mx)
plot(PolyBlock, col="white",border="white", add=TRUE)
### plot 2 legend ###
par(mai=c(0,.25,0,2))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
for(i in 1:(length(canopyBr)-1)){
  polygon(c(0,0,1,1), 
          c(canopyBr[i]/canopyBr[length(canopyBr)],canopyBr[i+1]/canopyBr[length(canopyBr)],canopyBr[i+1]/canopyBr[length(canopyBr)],canopyBr[i]/canopyBr[length(canopyBr)]),
          col=treePallete[i],border=NA)
}
axis(4,canopyBr/canopyBr[length(canopyBr)],canopyBr,cex.axis=cxa,las=2)	
mtext("Canopy cover (%)", side=3, line=1, cex=2)
dev.off()



###########################################
########## Figure 2&3: map of     #########
########## melt rate & swe max    #########
########## temp & doy start       #########
########## with violins           #########
########## Figure 2 & 3      -----
#break up names 
nameSplit1 <- character(0)
nameSplit2 <- character(0)
for(i in 1:5){
  nameSplit1[i] <- strsplit(glcID$name2[i], " ")[[1]][1]
  nameSplit2[i] <- strsplit(glcID$name2[i], " ")[[1]][2]
}
nameSplit2 <- ifelse(is.na(nameSplit2), " ", nameSplit2)


#average melt over 10 year time period
melt.ave <- calc(abs(melt.mm.day), fun=mean,na.rm=TRUE)
melt.sd <- calc(abs(melt.mm.day), fun=sd,na.rm=TRUE)
melt.count <- calc(melt.mm.day, fun=function(x){length(x[!is.na(x)])})


#average swe max over 10 year time period
max.ave <- calc(maxSwe, fun=mean,na.rm=TRUE)
max.mm.ave <- max.ave * 1000


#average doy start
doyS.ave <- calc(doyStart, fun=mean, na.rm=TRUE)


#average melt temp
temp.ave <- calc(meltMeanT, fun=mean,na.rm=TRUE)

#set up violin plots
#look at absoulute melt
analysisDF$absRate <- abs(analysisDF$melt.mm.day)
#convert max to mm
analysisDF$maxSwe.mm <- analysisDF$maxSwe.m * 1000

histL <- list()
densityH <- numeric(0)
maxS <- numeric(0)
minS <- numeric(0)
quant <- list()
for(i in 1:nrow(glcID)){
  
  histL[[i]] <- hist(analysisDF$absRate[analysisDF$glc == glcID$glc[i]], breaks=seq(0,32, by=2))
  #get max and min
  
  densityH[i] <- max(histL[[i]]$density)
  histL[[i]]$densityScale <-histL[[i]]$density*(0.5/ densityH[i])
  maxS[i] <- round(max(analysisDF$absRate[analysisDF$glc == glcID$glc[i]]),1)
  minS[i] <- floor(min(analysisDF$absRate[analysisDF$glc == glcID$glc[i]]))
  #create box plot quantiles

  quant[[i]] <- quantile(analysisDF$absRate[analysisDF$glc == glcID$glc[i]], probs=c(0.025,0.25,0.50,0.75,0.975))
}




#get quantiles for max


histLM <- list()
densityHM <- numeric(0)
maxSM <- numeric(0)
minSM <- numeric(0)
quantM <- list()
for(i in 1:nrow(glcID)){
  
  histLM[[i]] <- hist(analysisDF$maxSwe.mm[analysisDF$glc == glcID$glc[i]], breaks=seq(0,500, by=10))
  #get max and min
  
  densityHM[i] <- max(histLM[[i]]$density)
  histLM[[i]]$densityScale <-histLM[[i]]$density*(0.5/ densityHM[i])
  maxSM[i] <- round(max(analysisDF$maxSwe.mm[analysisDF$glc == glcID$glc[i]]),1)
  minSM[i] <- floor(min(analysisDF$maxSwe.mm[analysisDF$glc == glcID$glc[i]]))
  #create box plot quantiles
  
  quantM[[i]] <- quantile(analysisDF$maxSwe.mm[analysisDF$glc == glcID$glc[i]], probs=c(0.025,0.25,0.50,0.75,0.975))
}

quantT <- list()
for(i in 1:5){
  quantT[[i]] <- quantile(analysisDF$meltTempC[analysisDF$glc == glcID$glc[i]], probs=c(0.025,0.25,0.50,0.75,0.975),na.rm=TRUE)
}

histLT <- list()
densityHT <- numeric(0)
maxST <- numeric(0)
minST <- numeric(0)

for(i in 1:5){
  
  histLT[[i]] <- hist(analysisDF$meltTempC[analysisDF$glc == glcID$glc[i]], breaks=seq(-9,15, by=0.5))
  #get max and min
  
  densityHT[i] <- max(histLT[[i]]$density)
  histLT[[i]]$densityScale <-histLT[[i]]$density*(0.5/ densityHT[i])
  maxST[i] <- round(max(analysisDF$meltTempC[analysisDF$glc == glcID$glc[i]], na.rm=TRUE),1)
  minST[i] <- floor(min(analysisDF$meltTempC[analysisDF$glc == glcID$glc[i]]*10, na.rm=TRUE))/10
}


#get quantiles for onset
quantO <- list()
for(i in 1:5){
  quantO[[i]] <- quantile(analysisDF$doyStart[analysisDF$glc == glcID$glc[i]], probs=c(0.025,0.25,0.50,0.75,0.975),na.rm=TRUE)
}
histLO <- list()
densityHO <- numeric(0)
maxSO <- numeric(0)
minSO <- numeric()

for(i in 1:5){
  
  histLO[[i]] <- hist(analysisDF$doyStart[analysisDF$glc == glcID$glc[i]], breaks=seq(30,170, by=4))
  #get max and min
  
  densityHO[i] <- max(histLO[[i]]$density)
  histLO[[i]]$densityScale <-histLO[[i]]$density*(0.5/ densityHO[i])
  maxSO[i] <- round(max(analysisDF$doyStart[analysisDF$glc == glcID$glc[i]], na.rm=TRUE),1)
  minSO[i] <- floor(min(analysisDF$doyStart[analysisDF$glc == glcID$glc[i]], na.rm=TRUE))
}





#set up figures
sweBr <- round(getJenksBreaks(getValues(melt.ave),9),2)	
sweMaxBr <- round(getJenksBreaks(getValues(max.mm.ave),9),2)
meltTBr <- round(getJenksBreaks(getValues(temp.ave),9),2)
OnsetBr <- round(getJenksBreaks(getValues(doyS.ave),9),0)

hd <- 22
wd1 <- 22
wd2 <- 8

#size of panel label
mx <- 4
#line for panel label
pll <- 2.5
#size of axis
cxa <- 3	
#size of ticks
tlw <- 3

#plotting
xseqV <- seq(1,10, by=2)

wd1V <-30

ylV <- 0
yhV <- 31	
ylT <- -9
yhT <- 15
ylVC <- 0
yhVC <- 500
ylO <- 30
yhO <- 170
cax <- 1.75	
lcx <- 3
lll <- 6
lllx <- 7
al1 <- 2
al2 <- 4

##########Map plots part 1

png(paste0(plotDI,"\\figure2_maps_swe_p1.png"), width = 25, height = 20, units = "in", res=300)
layout(matrix(seq(1,6),ncol=3, byrow=TRUE), width=c(lcm(wd1),lcm(wd2),lcm(wd1V)),height=c(lcm(hd),lcm(hd)))	

### plot 1 swe ave ###
par(mai=c(.5,.5,.5,.5))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
#plot points
image(melt.ave,breaks=sweBr, col=swePallete, add=TRUE)

plot(PolyBlock, col="white",border="white", add=TRUE)
text(-4150000,4150000,"A",cex=mx)
### legent plot 1 swe ave ###
par(mai=c(0.5,0.5,0.5,2))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
for(i in 1:(length(sweBr)-1)){
  polygon(c(0,0,1,1), 
          c(i/length(sweBr),
            (i+1)/length(sweBr),
            (i+1)/length(sweBr),
            i/length(sweBr)),
          col=swePallete[i],border=NA)
}
axis(4,seq(1,length(sweBr))/length(sweBr),round(sweBr,0),cex.axis=cxa,las=2)	

plot(c(0,1),c(0,1), xlim=c(0,12), ylim=c(ylV,yhV), axes=FALSE, type="n", xlab = " ", ylab= " ",
     xaxs="i", yaxs="i")
for(j in 1:5){
  #if plot order needs to be changed do so here
  i <- j
  polygon(c(xseqV[j]+(0-histL[[i]]$densityScale[histL[[i]]$mids<=maxS[i] & histL[[i]]$mids >= minS[i]]), 
            rev(xseqV[j]+(histL[[i]]$densityScale[histL[[i]]$mids<=maxS[i] & histL[[i]]$mids >= minS[i]]))),
          c(histL[[i]]$mids[histL[[i]]$mids<=maxS[i] & histL[[i]]$mids >= minS[i]],
            rev(histL[[i]]$mids[histL[[i]]$mids<=maxS[i] & histL[[i]]$mids >= minS[i]])), 
          lwd=0.75,  col=vegePallete3[i])
  arrows(	xseqV[j],quant[[i]][1], xseqV[j],quant[[i]][5], code=0, lwd=1)
  polygon(c(xseqV[j]-0.15,xseqV[j]-0.15,xseqV[j]+0.15,xseqV[j]+0.15),
          c(quant[[i]][2],quant[[i]][4],quant[[i]][4],quant[[i]][2]),
          border=NA, col=rgb(0.25,0.25,0.25,0.5))
  
  arrows( xseqV[j]-0.15,quant[[i]][3], xseqV[j]+0.15,quant[[i]][3],code=0, lwd=4, col=vegePallete[i])	
}	

axis(1, xseqV, rep(" ",length(xseqV)),lwd.ticks=tlw)
axis(2, seq(0,31, by=2),rep(" ",length(seq(0,31, by=2))),lwd.ticks=tlw)
mtext(paste(nameSplit1),at=xseqV,side=1,line=al1,cex=cax)
mtext(paste(nameSplit2),at=xseqV,side=1,line=al2,cex=cax)
mtext(seq(0,31, by=2), at=seq(0,31, by=2), side=2, las=2, line=al1, cex=cax)
mtext(expression(paste("Melt rate (mm day"^"-1",")")), side=2, line=lll, cex=lcx)
mtext("Landcover type", side=1, line=lllx, cex=lcx)
text(0.5,30,"B",cex=mx)


### plot 2 swe max ###
par(mai=c(.5,.5,.5,.5))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
#plot points
image(max.mm.ave,breaks=sweMaxBr, col=swePallete, add=TRUE)

plot(PolyBlock, col="white",border="white", add=TRUE)
text(-4150000,4150000,"C",cex=mx)
### legent plot 1 swe ave ###
par(mai=c(0.5,0.5,0.5,2))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
for(i in 1:(length(sweMaxBr)-1)){
  polygon(c(0,0,1,1), 
          c(i/length(sweMaxBr),
            (i+1)/length(sweMaxBr),
            (i+1)/length(sweMaxBr),
            i/length(sweMaxBr)),
          col=swePallete[i],border=NA)
}
axis(4,seq(1,length(sweMaxBr))/length(sweMaxBr),round(sweMaxBr,0),cex.axis=cxa,las=2)		
par(mai=c(0.5,0.5,0.5,0.5))
plot(c(0,1),c(0,1), xlim=c(0,12), ylim=c(ylVC,yhVC), axes=FALSE, type="n", xlab = " ", ylab= " ",
     xaxs="i", yaxs="i")
for(j in 1:5){
  i <- j
  polygon(c(xseqV[j]+(0-histLM[[i]]$densityScale[histLM[[i]]$mids<=maxSM[i] & histLM[[i]]$mids>=minSM[i] ]), 
            rev(xseqV[j]+(histLM[[i]]$densityScale[histLM[[i]]$mids<=maxSM[i] & histLM[[i]]$mids>=minSM[i]]))),
          c(histLM[[i]]$mids[histLM[[i]]$mids<=maxSM[i] & histLM[[i]]$mids>=minSM[i]],
            rev(histLM[[i]]$mids[histLM[[i]]$mids<=maxSM[i] & histLM[[i]]$mids>=minSM[i]])), 
          lwd=0.75,  col=vegePallete3[i])
  arrows(	xseqV[j],quantM[[i]][1], xseqV[j],quantM[[i]][5], code=0, lwd=1)
  polygon(c(xseqV[j]-0.15,xseqV[j]-0.15,xseqV[j]+0.15,xseqV[j]+0.15),
          c(quantM[[i]][2],quantM[[i]][4],quantM[[i]][4],quantM[[i]][2]),
          border=NA, col=rgb(0.25,0.25,0.25,0.5))
  
  arrows( xseqV[j]-0.15,quantM[[i]][3], xseqV[j]+0.15,quantM[[i]][3],code=0, lwd=4, col=vegePallete[i])		
}	

axis(1, xseqV, rep(" ",length(xseqV)),lwd.ticks=tlw)
axis(2, seq(0,500, by=100),rep(" ",length(seq(0,500, by=100))),lwd.ticks=tlw)
mtext(paste(nameSplit1),at=xseqV,side=1,line=al1,cex=cax)
mtext(paste(nameSplit2),at=xseqV,side=1,line=al2,cex=cax)
mtext(seq(0,500, by=100), at=seq(0,500, by=100), side=2, las=2, line=al1, cex=cax)
mtext(expression(paste("Maximum SWE (mm)")), side=2, line=lll, cex=lcx)
mtext("Landcover type", side=1, line=lllx, cex=lcx)
text(0.5,470,"D",cex=mx)
dev.off()	

##########Map plots part 2

png(paste0(plotDI,"\\figure3_maps_swe_p2.png"), width = 25, height = 20, units = "in", res=300)
layout(matrix(seq(1,6),ncol=3, byrow=TRUE), width=c(lcm(wd1),lcm(wd2),lcm(wd1V)),height=c(lcm(hd),lcm(hd)))	

### plot 1 swe ave ###
par(mai=c(.5,.5,.5,.5))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
#plot points
image(temp.ave,breaks=meltTBr, col=swePallete, add=TRUE)

plot(PolyBlock, col="white",border="white", add=TRUE)
text(-4150000,4150000,"A",cex=mx)
### legent plot 1 swe ave ###
par(mai=c(0.5,0.5,0.5,2))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
for(i in 1:(length(meltTBr)-1)){
  polygon(c(0,0,1,1), 
          c(i/length(meltTBr),
            (i+1)/length(meltTBr),
            (i+1)/length(meltTBr),
            i/length(meltTBr)),
          col=swePallete[i],border=NA)
}
axis(4,seq(1,length(meltTBr))/length(meltTBr),round(meltTBr,1),cex.axis=cxa,las=2)	

plot(c(0,1),c(0,1), xlim=c(0,12), ylim=c(ylT,yhT), axes=FALSE, type="n", xlab = " ", ylab= " ",
     xaxs="i", yaxs="i")
for(j in 1:5){
  i <- j
  polygon(c(xseqV[j]+(0-histLT[[i]]$densityScale[histLT[[i]]$mids<=maxST[i] & histLT[[i]]$mids>=minST[i]]), 
            rev(xseqV[j]+(histLT[[i]]$densityScale[histLT[[i]]$mids<=maxST[i] & histLT[[i]]$mids>=minST[i]]))),
          c(histLT[[i]]$mids[histLT[[i]]$mids<=maxST[i] & histLT[[i]]$mids>=minST[i]],
            rev(histLT[[i]]$mids[histLT[[i]]$mids<=maxST[i] & histLT[[i]]$mids>=minST[i]])), 
          lwd=0.75,  col=vegePallete3[i])
  arrows(	xseqV[j],quantT[[i]][1], xseqV[j],quantT[[i]][5], code=0, lwd=1)
  polygon(c(xseqV[j]-0.15,xseqV[j]-0.15,xseqV[j]+0.15,xseqV[j]+0.15),
          c(quantT[[i]][2],quantT[[i]][4],quantT[[i]][4],quantT[[i]][2]),
          border=NA, col=rgb(0.25,0.25,0.25,0.5))
  
  arrows( xseqV[j]-0.15,quantT[[i]][3], xseqV[j]+0.15,quantT[[i]][3],code=0, lwd=4, col=vegePallete[i])	
}	

axis(1, xseqV, rep(" ",length(xseqV)),lwd.ticks=tlw)
axis(2, seq(-9,15, by=3),rep(" ",length(seq(-9,15, by=3))),lwd.ticks=tlw)
mtext(paste(nameSplit1),at=xseqV,side=1,line=al1,cex=cax)
mtext(paste(nameSplit2),at=xseqV,side=1,line=al2,cex=cax)
mtext(seq(-9,15, by=3), at=seq(-9,15, by=3), side=2, las=2, line=al1, cex=cax)
mtext(expression(paste("Average temperature (",degree,"C )")), side=2, line=lll, cex=lcx)
mtext("Landcover type", side=1, line=lllx, cex=lcx)
text(0.5,14.5,"B",cex=mx)
### plot 2 swe sd ###
par(mai=c(.5,.5,.5,.5))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ",xlim=c(-4150000,4150000),ylim=c(-4150000,4150000))
#color background
polygon(c(-5000000,-5000000,5000000,5000000),c(-5000000,5000000,5000000,-5000000), border=NA, col=water)
#boundaries
points(world, type="l", lwd=2, col="grey65")
#continent color
polygon(c(world[,1],rev(world[,1])), c(world[,2],rev(world[,2])),col=land,border=NA)
#plot points
image(doyS.ave,breaks=OnsetBr, col=swePallete, add=TRUE)

plot(PolyBlock, col="white",border="white", add=TRUE)
text(-4150000,4150000,"C",cex=mx)
### legent plot 1 swe ave ###
par(mai=c(0.5,0.5,0.5,2))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab=" ", ylab=" ", xlim=c(0,1),ylim=c(0,1)) 
for(i in 1:(length(OnsetBr)-1)){
  polygon(c(0,0,1,1), 
          c(i/length(OnsetBr),
            (i+1)/length(OnsetBr),
            (i+1)/length(OnsetBr),
            i/length(OnsetBr)),
          col=swePallete[i],border=NA)
}
axis(4,seq(1,length(OnsetBr))/length(OnsetBr),floor(OnsetBr),cex.axis=cxa,las=2)		
par(mai=c(0.5,0.5,0.5,0.5))
plot(c(0,1),c(0,1), xlim=c(0,12), ylim=c(ylO,yhO), axes=FALSE, type="n", xlab = " ", ylab= " ",
     xaxs="i", yaxs="i")
for(j in 1:5){
  i <- j
  polygon(c(xseqV[j]+(0-histLO[[i]]$densityScale[histLO[[i]]$mids<=maxSO[i] & histLO[[i]]$mids >= minSO[i]]), 
            rev(xseqV[j]+(histLO[[i]]$densityScale[histLO[[i]]$mids<=maxSO[i]& histLO[[i]]$mids >= minSO[i]]))),
          c(histLO[[i]]$mids[histLO[[i]]$mids<=maxSO[i]& histLO[[i]]$mids >= minSO[i]],
            rev(histLO[[i]]$mids[histLO[[i]]$mids<=maxSO[i]& histLO[[i]]$mids >= minSO[i]])), 
          lwd=0.75,  col=vegePallete3[i])
  arrows(	xseqV[j],quantO[[i]][1], xseqV[j],quantO[[i]][5], code=0, lwd=1)
  polygon(c(xseqV[j]-0.15,xseqV[j]-0.15,xseqV[j]+0.15,xseqV[j]+0.15),
          c(quantO[[i]][2],quantO[[i]][4],quantO[[i]][4],quantO[[i]][2]),
          border=NA, col=rgb(0.25,0.25,0.25,0.5))
  
  arrows( xseqV[j]-0.15,quantO[[i]][3], xseqV[j]+0.15,quantO[[i]][3],code=0, lwd=4, col=vegePallete[i])		
}	

axis(1, xseqV, rep(" ",length(xseqV)),lwd.ticks=tlw)
axis(2, seq(30,160, by=30),rep(" ",length(seq(30,160, by=30))),lwd.ticks=tlw)
mtext(paste(nameSplit1),at=xseqV,side=1,line=al1,cex=cax)
mtext(paste(nameSplit2),at=xseqV,side=1,line=al2,cex=cax)
mtext(seq(30,160, by=30), at=seq(30,160, by=30), side=2, las=2, line=al1, cex=cax)
mtext(expression(paste("Melt Onset Day")), side=2, line=lll, cex=lcx)
mtext("Landcover type", side=1, line=lllx, cex=lcx)
text(0.5,167,"D",cex=mx)
dev.off()



###########################################
########## Figure 4: bivariate   ##########
########## plots of regression   ##########
########## data                  ##########
########## Figure 4: -------

#organize model output
mu.Temp$gcID <- rep(seq(1,5),each=200) 
mu.Onset$gcID <- rep(seq(1,5),each=200) 
mu.Max$gcID <- rep(seq(1,5),each=200) 

#intercepts

plotTree <- c(1,2,3)	
plotTun <- c(4,5)			
analysisDF$logAbsRate <- log(abs(analysisDF$melt.mm.day))
wd1 <- 11
hd1 <- 11	
yl <- - 1
yh <- 3.5	
xl1 <- 	floor(range(analysisDF$meltTempC)[1])
xh1 <-	ceiling(range(analysisDF$meltTempC)[2])
xl2 <- floor(range(analysisDF$vcf)[1])
xh2 <-	ceiling(range(analysisDF$vcf)[2])
xl3 <- range(analysisDF$doyStart)[1] - 1
xh3 <-	range(analysisDF$doyStart)[2] + 1

###check name of swe Max
xl4 <- floor(range(analysisDF$log.max)[1])
xh4 <- round(range(analysisDF$log.max)[2],1)
#axis labels
xs1 <- seq(xl1,xh1-3, by=3)
xs2 <- seq(xl2,xh2, by= 15)
xs3 <- seq(60,xh3, by= 30)
ys <- seq(yl,yh, by = 1)
xs4 <- seq(xl4,xh4,by=0.5)

#width of regression line
mlw <- 4
#width of ticks
tlw <- 4
#axis tick label line
tll <- 2
#axis label size
alc <- 2
#plot label text size
plc <- 3
#x label plot line
xpl <- 6
#size of panel letter
ttx <- 4
#legend size
legcex <- 2.5

dlTemp <- numeric(0)
dhTemp <- numeric(0)
dlvcf<- numeric(0)
dhvcf <- numeric(0)
dlOnset <- numeric(0)
dhOnset <- numeric(0)
dlMax <- numeric(0)
dhMax <- numeric(0)
for(i in 1:5){
  dlTemp[i] <- floor(min(analysisDF$meltTempC[analysisDF$glc == glcID$glc[i]]))
  dhTemp[i] <- ceiling(max(analysisDF$meltTempC[analysisDF$glc == glcID$glc[i]]))
  dlvcf[i] <- floor(min(analysisDF$vcf[analysisDF$glc == glcID$glc[i]]))
  dhvcf[i] <- ceiling(max(analysisDF$vcf[analysisDF$glc == glcID$glc[i]]))
  dlOnset[i] <- floor(min(analysisDF$doyStart[analysisDF$glc == glcID$glc[i]]))
  dhOnset[i] <- ceiling(max(analysisDF$doyStart[analysisDF$glc == glcID$glc[i]]))
  dlMax[i] <- floor(min(analysisDF$log.max[analysisDF$glc == glcID$glc[i]])*10)/10
  dhMax[i] <- ceiling(max(analysisDF$log.max[analysisDF$glc == glcID$glc[i]])*10)/10
  
  
}




png(paste0(plotDI,"\\regression.png"), width = 55, height = 32, units = "cm", res=300)
layout(matrix(seq(1,8),ncol=4, byrow=TRUE), width=rep(lcm(wd1),4),height=rep(lcm(hd1),2))
par(mai=c(0,0,0,0))
#temperature trees
plot(c(0,1),c(0,1), type="n", xlim=c(xl1,xh1), ylim=c(yl,yh), xaxs="i",yaxs="i",
     xlab= " ", ylab=" ", axes=FALSE)
for(i in plotTree){
  points(	analysisDF$meltTempC[analysisDF$glc == glcID$glc[i]],
          analysisDF$logAbsRate[analysisDF$glc == glcID$glc[i]], col=vegePallete2[i], pch=19)
}
for(i in plotTree){	
  polygon(c(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
            rev(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]])),
          c(mu.Temp$X2.5.[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
            rev(mu.Temp$X97.5[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]])),
          border=NA, col=vegePallete3[i])
  
  points(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
         mu.Temp$Mean[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
         type="l", lwd=mlw, col=vegePallete[i])
}
axis(2, ys, rep(" ",length(ys)), lwd.ticks=tlw)
mtext(ys,at=ys, line=tll, cex=alc, side=2,las=2)
box(which="plot")
mtext(expression(paste("log(Melt Rate (cm day"^"-1","))")), side=2, outer=TRUE,line= -5, cex=plc)
text(xl1+(.05*(xh1-xl1)), yh-(.05*(yh-yl)), "a", cex=ttx)
legend("bottomright", paste(glcID$name2[plotTree]), col=vegePallete[plotTree],cex=legcex, lwd=mlw,lty=1, bty="n")
#tree cover trees
par(mai=c(0,0,0,0))	
plot(c(0,1),c(0,1), type="n", xlim=c(xl2,xh2), ylim=c(yl,yh), xaxs="i",yaxs="i",
     xlab= " ", ylab=" ", axes=FALSE)
for(i in plotTree){
  points(	analysisDF$vcf[analysisDF$glc == glcID$glc[i]],
          analysisDF$logAbsRate[analysisDF$glc == glcID$glc[i]], col=vegePallete2[i], pch=19)
}
for(i in plotTree){	
  polygon(c(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]],
            rev(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
          c(rep(beta0$X2.5.[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
            rep(beta0$X97.5[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]]))),
          border=NA, col=vegePallete3[i])
  
  points(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]],
         rep(beta0$Mean[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
         type="l", lwd=mlw, lty=3, col=vegePallete[i])		
  
}
box(which="plot")
text(xl2+(.05*(xh2-xl2)), yh-(.05*(yh-yl)), "c", cex=ttx)
#day of onset trees
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1), type="n", xlim=c(xl3,xh3), ylim=c(yl,yh), xaxs="i",yaxs="i",
     xlab= " ", ylab=" ", axes=FALSE)
for(i in plotTree){
  
  points(analysisDF$doyStart[analysisDF$glc == glcID$glc[i]],
          analysisDF$logAbsRate[analysisDF$glc == glcID$glc[i]], col=vegePallete2[i], pch=19)
}
for(i in plotTree){	
  polygon(c(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
            rev(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]])),
          c(mu.Onset$X2.5.[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
            rev(mu.Onset$X97.5[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]])),
          border=NA, col=vegePallete3[i])
  
  points(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
         mu.Onset$Mean[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
         type="l", lwd=mlw, col=vegePallete[i])
}
box(which="plot")
text(xl3+(.05*(xh3-xl3)), yh-(.05*(yh-yl)), "e", cex=ttx)

#max swe trees
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1), type="n", xlim=c(xl4,xh4), ylim=c(yl,yh), xaxs="i",yaxs="i",
     xlab= " ", ylab=" ", axes=FALSE)
for(i in plotTree){
  points(	analysisDF$log.max[analysisDF$glc == glcID$glc[i]],
          analysisDF$logAbsRate[analysisDF$glc == glcID$glc[i]], col=vegePallete2[i], pch=19)
}
for(i in plotTree){	
  MaxMeans <- MaxMean[MaxMean >= dlMax[i] & MaxMean <= dhMax[i]]
  mu.Maxs <- mu.Max[MaxMean >= dlMax[i] & MaxMean <= dhMax[i],]
  polygon(c(MaxMeans,
            rev(MaxMeans)),
          c(mu.Maxs$X2.5.[mu.Maxs$gcID == i],
            rev(mu.Maxs$X97.5[mu.Maxs$gcID == i])),
          border=NA, col=vegePallete3[i])
  
  points(MaxMeans,
         mu.Maxs$Mean[mu.Maxs$gcID == i],
         type="l", lwd=mlw, col=vegePallete[i])	
  
}
box(which="plot")
text(xl4+(.05*(xh4-xl4)), yh-(.05*(yh-yl)), "g", cex=ttx)	


#temperature tundra
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1), type="n", xlim=c(xl1,xh1), ylim=c(yl,yh), xaxs="i",yaxs="i",
     xlab= " ", ylab=" ", axes=FALSE)
for(i in plotTun){
  points(	analysisDF$meltTempC[analysisDF$glc == glcID$glc[i]],
          analysisDF$logAbsRate[analysisDF$glc == glcID$glc[i]], col=vegePallete2[i], pch=19)
}
for(i in plotTun){
  polygon(c(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
            rev(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]])),
          c(mu.Temp$X2.5.[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
            rev(mu.Temp$X97.5[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]])),
          border=NA, col=vegePallete3[i])
  
  points(tempMean[tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
         mu.Temp$Mean[mu.Temp$gcID == i & tempMean >= dlTemp[i] & tempMean <= dhTemp[i]],
         type="l", lwd=mlw, col=vegePallete[i])		
}	
axis(2, ys, rep(" ",length(ys)), lwd.ticks=tlw)
mtext(ys,at=ys, line=tll, cex=alc, side=2,las=2)
axis(1, xs1, rep(" ",length(xs1)), lwd.ticks=tlw)
mtext(xs1,at=xs1, line=tll, cex=alc, side=1)
box(which="plot")
mtext("Temperature (c)", side=1,line= xpl, cex=plc)
text(xl1+(.05*(xh1-xl1)), yh-(.05*(yh-yl)), "b", cex=ttx)
legend("bottomright", paste(glcID$name2[plotTun]), col=vegePallete[plotTun],cex=legcex, lwd=mlw,lty=1, bty="n")
#tree cover tundra
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1), type="n", xlim=c(xl2,xh2), ylim=c(yl,yh), xaxs="i",yaxs="i",
     xlab= " ", ylab=" ", axes=FALSE)
for(i in plotTun){
  points(		analysisDF$vcf[analysisDF$glc == glcID$glc[i]],
           analysisDF$logAbsRate[analysisDF$glc == glcID$glc[i]], col=vegePallete2[i], pch=19)
}
for(i in plotTun){	
  polygon(c(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]],
            rev(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
          c(rep(beta0$X2.5.[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
            rep(beta0$X97.5[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]]))),
          border=NA, col=vegePallete3[i])
  
  points(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]],
         rep(beta0$Mean[i], length(CanopyMean[CanopyMean >= dlvcf[i] & CanopyMean <= dhvcf[i]])),
         type="l", lwd=mlw, lty=3, col=vegePallete[i])		
  
}
axis(1, xs2, rep(" ",length(xs2)), lwd.ticks=tlw)
mtext(xs2,at=xs2, line=tll, cex=alc, side=1)
box(which="plot")
mtext("Canopy cover (%)", side=1,line= xpl, cex=plc)
text(xl2+(.05*(xh2-xl2)), yh-(.05*(yh-yl)), "d", cex=ttx)
#onset day tundra
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1), type="n", xlim=c(xl3,xh3), ylim=c(yl,yh), xaxs="i",yaxs="i",
     xlab= " ", ylab=" ", axes=FALSE)
for(i in plotTun){
  points(	analysisDF$doyStart[analysisDF$glc == glcID$glc[i]],
          analysisDF$logAbsRate[analysisDF$glc == glcID$glc[i]], col=vegePallete2[i], pch=19)
}
for(i in plotTun){	
  polygon(c(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
            rev(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]])),
          c(mu.Onset$X2.5.[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
            rev(mu.Onset$X97.5[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]])),
          border=NA, col=vegePallete3[i])
  
  points(SdayMean[SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
         mu.Onset$Mean[mu.Onset$gcID == i & SdayMean >= dlOnset[i] & SdayMean <= dhOnset[i]],
         type="l", lwd=mlw, col=vegePallete[i])
}
axis(1, xs3, rep(" ",length(xs3)), lwd.ticks=tlw)
mtext(xs3,at=xs3, line=tll, cex=alc, side=1)
box(which="plot")
mtext("Onset day of year", side=1,line= xpl, cex=plc)
text(xl3+(.05*(xh3-xl3)), yh-(.05*(yh-yl)), "f", cex=ttx)
#latitude tundra
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1), type="n", xlim=c(xl4,xh4), ylim=c(yl,yh), xaxs="i",yaxs="i",
     xlab= " ", ylab=" ", axes=FALSE)
for(i in plotTun){
  points(	analysisDF$log.max[analysisDF$glc == glcID$glc[i]],
          analysisDF$logAbsRate[analysisDF$glc == glcID$glc[i]], col=vegePallete2[i], pch=19)
}


for(i in plotTun){	
  MaxMeans <- MaxMean[MaxMean >= dlMax[i] & MaxMean <= dhMax[i]]
  mu.Maxs <- mu.Max[MaxMean >= dlMax[i] & MaxMean <= dhMax[i],]
  polygon(c(MaxMeans,
            rev(MaxMeans)),
          c(mu.Maxs$X2.5.[mu.Maxs$gcID == i],
            rev(mu.Maxs$X97.5[mu.Maxs$gcID == i])),
          border=NA, col=vegePallete3[i])
  
  points(MaxMeans,
         mu.Maxs$Mean[mu.Maxs$gcID == i],
         type="l", lwd=mlw, col=vegePallete[i])	
  
}
box(which="plot")
text(xl4+(.05*(xh4-xl4)), yh-(.05*(yh-yl)), "h", cex=ttx)	

axis(1, xs4, rep(" ",length(xs4)), lwd.ticks=tlw)
mtext(xs4,at=xs4, line=tll, cex=alc, side=1)
mtext("log(Max SWE )", side=1,line= xpl, cex=plc)
mtext("log(m)", side=1,line= xpl+5, cex=plc)
dev.off()


