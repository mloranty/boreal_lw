################################
#
# pre-AGU data anlaysis for 
# boreal LW project
#
# MML 12/8/16 (YIKES!)
##############################

require(raster)
require(rgdal)
require(plyr)
rm(list=ls())

setwd('/Users/mloranty/Google Drive/Documents/Research/Picker_Boreal/agu_analyses/')
## read all data and make histogram of tsr
all.dat <- read.csv('all.csv',header=T)
tsr <- all.dat$swe[all.dat$swe<50]

pdf(file='tsr_hist.pdf',6,6)
par(cex=1.75)
hist(tsr,col='blue',
     ylab="Frequency",main="",
     xlab=expression(paste('Mar-Apr TSR (mm ',degree*C^-1,")",sep="")))
dev.off()
###########################################

pdf(file='tsr_vcf_all.pdf',6,6)
par(cex=1.5)
plot(all.dat$vcf,all.dat$swe,xlab='Tree Cover (%)',
     ylab=expression(paste('TSR (mm ',degree*C^-1,")",sep="")),
     ylim=c(0,50))
dev.off()

## read in some other dat ##

glc <- raster('GLC2000_EASE_from_0.25_50_mask.tif')
vcf <- raster('MOD44B_2014_mosaic_ease.tif')

leg <- read.csv('initial_model_run/lc.csv',header=T) 
leg <- leg[order(leg$lc),]

yr <- read.csv('initial_model_run/year.csv', header=T)
yr <- yr[order(yr$year),]
mod.res <- read.csv('initial_model_run/glc50_quantile.csv',header=T)
mod.res$zone <-  leg$lc[match(mod.res$id,leg$lc.id)]

veg <- data.frame(cbind(getValues(glc),getValues(vcf)))
names(veg) <- c('glc','vcf')
veg <- na.omit(veg)


veg.sum <- aggregate(veg$vcf,by=list(veg$glc),FUN=mean)
names(veg.sum) <- c('zone','mean')
veg.sum$sd <- aggregate(veg$vcf,by=list(veg$glc),FUN=sd)[,2]

# cell count to see which vegetation zones are represented in the aggregate data set
veg.sum$count <- 0
for(i in 1:length(veg.sum$zone))
{
  veg.sum$count[i] <- length(which(veg$glc==veg.sum$zone[i]))
}

# calculate area in millions of km2
veg.sum$area <- veg.sum$count*(25.06752^2)/10^6
  
### USE THESE CLASSES/ZONES FOR ANALYSIS ###
z <- c(2,4,5)

col <- c('lightgreen','darkgreen','orange')
##################################################################################
##################################################################################
# Histograms of VCF for each vegetation class
# make hisotgrams for each of these jams
for(i in 1:length(z))
{
  
  d <- veg$vcf[which(veg$glc==z[i])]
  pdf(file=paste('zone',veg.sum$zone[which(veg.sum$zone==z[i])],
                 'hist.pdf',sep='.'),6,6)
  par(cex=1.75)
  hist(d,xlab='Tree Cover (%)',ylab="Frequency",col=col[i],
       main=leg$abr[leg$lc==z[i]])
  dev.off()
}
######### END ######

##################################################################################
##################################################################################
# scaterplots by veg type and year
all.2010 <- all.dat[all.dat$year==2010,]

pdf(file="annual_tsf_veg'pdf",8.5,11)
par(mfcol=c(3,2),cex=.9)

plot(all.2010$vcf[all.2010$lc==2],all.2010$swe[all.2010$lc==2],col=col[1],
     xlab= 'Tree Cover (%)',ylim=c(0,20),xlim=c(0,80),
     ylab=expression(paste('TSR (mm ',degree*C^-1,")",sep="")))
l <- lm(all.2010$swe[all.2010$lc==2]~all.2010$vcf[all.2010$lc==2])
abline(l,lty='dashed')

plot(all.2010$vcf[all.2010$lc==4],all.2010$swe[all.2010$lc==4],col=col[2],
     xlab= 'Tree Cover (%)',ylim=c(0,20),xlim=c(0,80),
     ylab=expression(paste('TSR (mm ',degree*C^-1,")",sep="")))
l <- lm(all.2010$swe[all.2010$lc==4]~all.2010$vcf[all.2010$lc==4])
abline(l,lty='dashed')
summary(l)

plot(all.2010$vcf[all.2010$lc==5],all.2010$swe[all.2010$lc==5],col=col[3],
     xlab= 'Tree Cover (%)',ylim=c(0,20),xlim=c(0,80),
     ylab=expression(paste('TSR (mm ',degree*C^-1,")",sep="")))
l <- lm(all.2010$swe[all.2010$lc==5]~all.2010$vcf[all.2010$lc==5])
abline(l,lty='dashed')
summary(l)

# plot(all.2010$vcf[all.2010$lc==6],all.2010$swe[all.2010$lc==6],col=col[4],
#      xlab= 'Tree Cover (%)',ylim=c(0,20),xlim=c(0,80),
#      ylab=expression(paste('TSR (mm ',degree*C^-1,")",sep="")))
# l <- lm(all.2010$swe[all.2010$lc==6]~all.2010$vcf[all.2010$lc==6])
# abline(l,lty='dashed')
# summary(l)

## 2011 ##
all.2011 <- all.dat[all.dat$year==2011,]
plot(all.2011$vcf[all.2011$lc==2],all.2011$swe[all.2011$lc==2],col=col[1],
     xlab= 'Tree Cover (%)',ylim=c(0,20),xlim=c(0,80),
     ylab=expression(paste('TSR (mm ',degree*C^-1,")",sep="")))
l <- lm(all.2011$swe[all.2011$lc==2]~all.2011$vcf[all.2011$lc==2])
abline(l,lty='dashed')
summary(l)

plot(all.2011$vcf[all.2011$lc==4],all.2011$swe[all.2011$lc==4],col=col[2],
     xlab= 'Tree Cover (%)',ylim=c(0,20),xlim=c(0,80),
     ylab=expression(paste('TSR (mm ',degree*C^-1,")",sep="")))
l <- lm(all.2011$swe[all.2011$lc==4]~all.2011$vcf[all.2011$lc==4])
abline(l,lty='dashed')
summary(l)

plot(all.2011$vcf[all.2011$lc==5],all.2011$swe[all.2011$lc==5],col=col[3],
     xlab= 'Tree Cover (%)',ylim=c(0,20),xlim=c(0,80),
     ylab=expression(paste('TSR (mm ',degree*C^-1,")",sep="")))
l <- lm(all.2011$swe[all.2011$lc==5]~all.2011$vcf[all.2011$lc==5])
abline(l,lty='dashed')
summary(l)
dev.off()
# plot(all.2011$vcf[all.2011$lc==6],all.2011$swe[all.2011$lc==6],col=col[4],
#      xlab= 'Tree Cover (%)',ylim=c(0,20),xlim=c(0,80),
#      ylab=expression(paste('TSR (mm ',degree*C^-1,")",sep="")))
# l <- lm(all.2011$swe[all.2011$lc==6]~all.2011$vcf[all.2011$lc==6])
# abline(l,lty='dashed')
# summary(l)
##################################################################################
##################################################################################
# plots of Bayesian Regression Parameters

## intercept barplot
pdf(file='intercept.pdf',6,6)
par(cex=1.5)
int.res <- mod.res[1:17,]  
ord <- match(z,int.res$zone)
int.bar <- barplot(int.res$Mean[ord],horiz=F,names=leg$abr[match(z,leg$lc)],col=col,
                   ylab=expression(paste("Intercept (mm",degree*C^-1,")",sep="")))
arrows(int.bar,int.res$X2.50.[ord],int.bar,int.res$X97.50.[ord],length=0,lwd=3)
dev.off()

## intercept vs. tree cover
pdf(file='intercept_vcf.pdf',6,6)
plot(veg.sum$mean[match(z,veg.sum$zone)],int.res$Mean[ord],col=col,pch=16,
     xlim=c(0,65),ylim=c(0,10),
     xlab="Tree Cover (%)",
     ylab="Intercept")

l <- lm(int.res$Mean[ord]~veg.sum$mean[match(z,veg.sum$zone)])
abline(l,lty='dashed')
dev.off()

## year barplot
pdf(file='year.pdf',10,5)
par(cex=1.5)
yr.res <- mod.res[35:69,]
yr.bar <- barplot(yr.res$Mean[yr$year.id],ylim=c(-1.75,2),
                  names=yr$year)
arrows(yr.bar,yr.res$X2.50.[yr$year.id],yr.bar,yr.res$X97.50.[yr$year.id],lwd=2,length=0)
dev.off()

## slope barplot
pdf(file='slope.pdf',6,6)
par(cex=1.5)
slp.res <- mod.res[18:34,]
ord <- match(z,slp.res$zone)
slp.bar <- barplot(slp.res$Mean[ord],col=col,names=leg$abr[match(z,leg$lc)],
                   ylab=expression(paste("Slope (mm",degree*C^-1,")",sep="")))
arrows(slp.bar,slp.res$X2.50.[ord],slp.bar,slp.res$X97.50.[ord],lwd=2,length=0)
dev.off()
######### END ######


