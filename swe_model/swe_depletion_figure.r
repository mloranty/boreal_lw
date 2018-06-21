#######################################################
# Script for making figures                           #
# a model on a blended swe                            #
# product to look at differences in timing of         #
# swe depletion                                       #
#######################################################

library(plyr)
library(coda)
library(mcmcplots)

############################################
###  model input directory               ###
############################################
modDI <- "z:\\Projects\\boreal_swe_depletion\\model\\run5"
plotDI <- "z:\\Projects\\boreal_swe_depletion\\figures\\model\\run5"

############################################
###  read in swe data  and organize      ###
############################################
dat.swe <- read.csv("z:\\Projects\\boreal_swe_depletion\\data\\swe_depletion_model_data_vcf_no_topo.csv")
dat.glc <- read.csv("z:\\Projects\\boreal_swe_depletion\\data\\glc50_table.csv")

#calculate proportion of land cover
dat.swe$glc1.p <- dat.swe$glc1f/3136
hist(dat.swe$glc1.p )
dat.swe$glc2.p <- dat.swe$glc2f/3136
hist(dat.swe$glc2.p )

#just use first glc class 				
dat.swe$zone <- dat.swe$glc1


##########################
##### Subset point 1 #####
##########################
#only focus on 2000-2009 for now
dat.swe <- dat.swe[dat.swe$year==2009,]

##########################
##### Filter point 1 #####
##########################
#start with the simplest classification of glc 
#>50% is dominated by the landcover class
dat.swe1 <- dat.swe[dat.swe$glc1.p>=.5,]


##########################
##### Filter point 2 #####
##########################
#filter out land cover types not of interest
#don't include glc 3,9,10,14,15,16,18,19,20,21,22
dat.swe2 <- dat.swe1[dat.swe1$zone!=2&dat.swe1$zone!=3&dat.swe1$zone!=15&dat.swe1$zone!=16&dat.swe1$zone!=10&dat.swe1$zone!=14&
				dat.swe1$zone!=9&dat.swe1$zone<18,]

##########################
##### Filter point 3 #####
##########################				
				
#filter out areas with no snow extent in spring
sweMax <- aggregate(dat.swe2$swe,by=list(dat.swe2$cell,dat.swe2$year),FUN="max")
colnames(sweMax) <- c("cell","year","sweMax")
sweMin <- aggregate(dat.swe2$swe,by=list(dat.swe2$cell,dat.swe2$year),FUN="min")
colnames(sweMin) <- c("cell","year","sweMin")
sweMax$Diff <- sweMax$sweMax-sweMin$sweMin
te <- hist(sweMin$sweMin, breaks=seq(0,0.5, by=0.01))
te2 <- hist(sweMax$Diff, breaks=seq(0,0.55, by=0.01))
te3 <- hist(sweMax$sweMax, breaks=seq(0,0.6, by=0.01))
#vast majority of minimum swe 97% is below 0.02 in a year
#exclude any sites and years where the maximum swe does not exceed 0.04 in the year
sweMaxF <- sweMax[sweMax$sweMax >=0.01,]
#join back to swe to filter
dat.swe3 <- join(dat.swe2, sweMaxF, by=c("cell","year"), type="inner")


#######################################################
# organize data for model run                         #
#######################################################
#create a table of identifiers
#unique glc

IDSglc <- unique(data.frame(zone=dat.swe3$zone))
IDSglc <- join(IDSglc, dat.glc[,1:2], by="zone", type="left")
IDSglc$gcID <- seq(1,dim(IDSglc)[1])


#join glc ID into dataframe
dat.swe4 <- join(dat.swe3,IDSglc, by="zone", type="left")
dat.swe4$t.airC <- dat.swe4$t.air-273.15

#get average temperature over a pixel, year

temp.py <- aggregate(dat.swe4$t.airC, by=list(dat.swe4$cell,dat.swe4$year,dat.swe4$zone), FUN="mean")
colnames(temp.py) <- c("cell","year","zone","temp")

#calculate average to center for each zone
temp.z <- aggregate(temp.py$temp,by=list(temp.py$zone),FUN="mean")
colnames(temp.z) <- c("zone","temp.zoneM")
#join zone mean back into temp.py
temp.py <- join(temp.py, temp.z, by="zone", type="left")
temp.py$tempCent <- temp.py$temp-temp.py$temp.zoneM

#join back into swe data
dat.swe5 <- join(dat.swe4,temp.py, by=c("cell","year","zone"), type="left")

#normalize swe
#calculate percent
dat.swe5$sweP <- dat.swe5$swe/dat.swe5$sweMax

#round swe for 20% of peak to 1 and 20% of low to zero
dat.swe5$sweN <- ifelse(dat.swe5$sweP>=0.8,1,
					ifelse(dat.swe5$sweP<=0.2,0,dat.swe5$sweP))



#get unique pixel id in each glc for parameter id

pixID <- unique(data.frame(cell=dat.swe5$cell, gcID=dat.swe5$gcID))

#make sure pix doesn't change class
length(unique(dat.swe5$cell))
 dim(pixID)

#subset into each glc
pixList <- list()
pixGLC <- numeric(0)
for(i in 1:dim(IDSglc)[1]){
	pixList[[i]] <- pixID[pixID$gcID==i,]
	pixList[[i]]$pixID <- seq(1,dim(pixList[[i]])[1])
	pixGLC[i] <- dim(pixList[[i]])[1]
}

pixJ <- ldply(pixList,data.frame)

#join back into swe
dat.swe6 <- join(dat.swe5,pixJ, by=c("cell","gcID"), type="left")


print("finish data organize")
############################################
###  read in model input                 ###
############################################
#list all output files

modFiles <- list.files(modDI)

#pull out parameter name
modParam <- character(0)
modGLC <- numeric(0)
modChainT <- character(0)
modChain <- numeric(0)
for(i in 1:length(modFiles)){
	modParam[i] <- strsplit(modFiles, "\\_")[[i]][1]
	modGLC[i] <- as.numeric(gsub("\\D","",strsplit(modFiles, "\\gc")[[i]][2]))
	modChainT[i] <- strsplit(modFiles, "\\gc")[[i]][1]
	modChain[i] <- as.numeric(gsub("\\D","",strsplit(modChainT,"\\chain")[[i]][2]))
	
}	

modelOut <- data.frame(parms=modParam,glc=modGLC,chain=modChain)
modelOut <- modelOut[modelOut$parms!="sig",]
#subset to only look at glc 1 for now

#pull out each model chain

chain1 <- which(modelOut$chain==1)
chain2 <- which(modelOut$chain==2)
chain3 <- which(modelOut$chain==3)

#join number of parms to model info
pixGLCDF <- data.frame(glc=seq(1,6),nparms=pixGLC)

modelOut <- join(modelOut,pixGLCDF, by="glc",type="left")
#read in coda
chain1Out <- list()
chain2Out <- list()
chain3Out <- list()
for(i in 1:length(chain1)){
	chain1Out[[i]] <- read.csv(paste0(modDI,"\\",modFiles[chain1[i]]))
	chain2Out[[i]] <- read.csv(paste0(modDI,"\\",modFiles[chain2[i]]))
	chain3Out[[i]] <- read.csv(paste0(modDI,"\\",modFiles[chain3[i]]))
	colnames(chain1Out[[i]]) <- paste0(modelOut$parms[chain1[i]],"_",modelOut$glc[chain1[i]],"_",seq(1,modelOut$nparms[chain1[i]]))
	colnames(chain2Out[[i]]) <-  paste0(modelOut$parms[chain2[i]],"_",modelOut$glc[chain2[i]],"_",seq(1,modelOut$nparms[chain2[i]]))
	colnames(chain3Out[[i]]) <-  paste0(modelOut$parms[chain3[i]],"_",modelOut$glc[chain3[i]],"_",seq(1,modelOut$nparms[chain3[i]]))
}
############################################
###  turn into a coda object             ###
############################################
#organize coda

chain1all <- list()
chain2all <- list()
chain3all <- list()
codaAll <- list()
modSum <- list()

############################################
###  check history                       ###
###  and model output                    ###
############################################
for(i in 1:length(chain1)){
	chain1all[[i]] <- as.mcmc(chain1Out[[i]])
	chain2all[[i]]  <- as.mcmc(chain2Out[[i]] )
	chain3all[[i]]  <- as.mcmc(chain3Out[[i]] )
	codaAll <- mcmc.list(chain1all[[i]],chain2all[[i]],chain3all[[i]])
	modSum[[i]] <- summary(codaAll)
	#mcmcplot(codaAll,dir=paste0(plotDI,"\\history") )
}

datC <- list()

for(i in 1:length(chain1)){
	datC[[i]] <- data.frame(cbind(modSum[[i]]$statistics,modSum[[i]]$quantiles))
	datC[[i]]$glcID <- rep(modelOut$glc[chain1[i]],dim(datC[[i]])[1])
	datC[[i]]$parms <- rep(modelOut$parms[chain1[i]],dim(datC[[i]])[1])
	
}

IDSglc$name <- c("needleleaf deciduous","deciduous shrub","evergreen shrub","herbaceous",
					"needleleaf evergreen","mixed tree")

	
#join pixel info to data
#first add pixel id to the data frame
b0glc <- list()
for(i in 1:6){
	b0glc[[i]] <- cbind(datC[[i]],pixList[[i]])

}
mid0glc <- list()
for(i in 1:6){
	mid0glc[[i]] <- cbind(datC[[i+6]],pixList[[i]])

}
#organize into data frame
#get tree cover data
treeDF <- unique(data.frame(cell=dat.swe6$cell,vcf=dat.swe6$vcf))


#calculate max swe for each cell
maxSWE <- aggregate(dat.swe6$swe, by=list(dat.swe6$cell),FUN="max")
colnames(maxSWE) <- c("cell","max.swe")

#covatiate data
covDF <- join(temp.py,treeDF, by="cell",type="left")
covDF <- join(covDF,maxSWE, by="cell",type="left")

for(i in 1:6){
	b0glc[[i]] <- join(b0glc[[i]],covDF,by="cell",type="left")
	mid0glc[[i]] <- join(mid0glc[[i]],covDF,by="cell",type="left")
}
	
############################################
###  plot relationships                  ###
############################################				
					
#plotting info					
wd <- 35
hd <- 35

xl1 <- -30
xh1 <- 5
yl1 <- 0
yh1 <- 100
xl2 <- 0
xh2 <- 80
yl2 <- 0
yh2 <- 1
xl3 <- 0
xh3 <- 0.6

xS1 <- seq(-30,0, by=10)
yS1 <- seq(0,100,by=10)

xS2 <- seq(0,80, by=10)
yS2 <- seq(0,.9,by=.1)
xS3 <- seq(.1,.6,by=.1)
px <- 8
llw <- 15
ltw <- 6
lmx <- 6
yll <- 6
xll <- 6
tx <- 13
#dim(IDSglc)[1]
for(i in 1:6){
	fit1 <- lm(b0glc[[i]]$Mean~b0glc[[i]]$temp)
	fit2 <- lm(b0glc[[i]]$Mean~b0glc[[i]]$vcf)
	fit3 <- lm(mid0glc[[i]]$Mean~mid0glc[[i]]$temp)
	fit4 <- lm(mid0glc[[i]]$Mean~mid0glc[[i]]$vcf)
	fit5 <- lm(b0glc[[i]]$Mean~b0glc[[i]]$max.swe)
	fit6 <- lm(mid0glc[[i]]$Mean~mid0glc[[i]]$max.swe)
	
	jpeg(paste0(plotDI,"\\relationships_",IDSglc$name[i],".jpeg"), width=4500,height=2500,quality=100)
	layout(matrix(seq(1,6),ncol=3, byrow=TRUE),width=rep(lcm(wd),6),height=rep(lcm(hd),6))
	par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl1,xh1), ylim=c(yl1,yh1), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")				
		points(b0glc[[i]]$temp,b0glc[[i]]$Mean, pch=19, cex=px,col=rgb(0,0,0,.15))
		abline(fit1, col="firebrick2", lwd=llw)
	
	
	mtext("swe melt slope",side=2,line=40,cex=10)
	mtext("(b)",side=2,line=25,cex=10)
	
	text(xl1+11, yl1+10,paste("R2=",round(summary(fit1)$r.squared,4)),cex=tx,col="firebrick2")
	text(xl1+11, yl1+25,paste("p=",round(summary(fit1)$coefficients[2,4],4)),cex=tx,col="firebrick2")
	
	mtext(IDSglc$name[i],side=3, outer=TRUE,line=-15,cex=10)
	axis(2, yS1, rep("", length(yS1)), lwd.ticks=ltw)
	mtext(yS1, at=yS1, side=2, cex=lmx,las=2, line=yll)

	
	box(which="plot")				
	par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl2,xh2), ylim=c(yl1,yh1), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")				
		points(b0glc[[i]]$vcf,b0glc[[i]]$Mean, pch=19, cex=px,col=rgb(0,0,0,.15))				
		abline(fit2, col="firebrick2", lwd=llw)			
	
	text(xl2+31, yl1+10,paste("R2=",round(summary(fit2)$r.squared,4)),cex=tx,col="firebrick2")
	text(xl2+31, yl1+25,paste("p=",round(summary(fit2)$coefficients[2,4],4)),cex=tx,col="firebrick2")
	
	box(which="plot")
	
	par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl3,xh3), ylim=c(yl1,yh1), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")				
		points(b0glc[[i]]$max.swe,b0glc[[i]]$Mean, pch=19, cex=px,col=rgb(0,0,0,.15))			
	abline(fit5, col="firebrick2", lwd=llw)	
	axis(4, yS1, rep("", length(yS1)), lwd.ticks=ltw)
	mtext(yS1, at=yS1, side=4, cex=lmx,las=2, line=yll)
	text(xl3+.21, yl1+10,paste("R2=",round(summary(fit5)$r.squared,4)),cex=tx,col="firebrick2")
	text(xl3+.21, yl1+25,paste("p=",round(summary(fit5)$coefficients[2,4],4)),cex=tx,col="firebrick2")
	mtext("swe melt slope",side=4,line=25,cex=10)
	mtext("(b)",side=4,line=40,cex=10)
	box(which="plot")					

	box(which="plot")				
	par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl1,xh1), ylim=c(yl2,yh2), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")				
					
		points(mid0glc[[i]]$temp,mid0glc[[i]]$Mean, pch=19, cex=px,col=rgb(0,0,0,.15))	
		abline(fit3, col="firebrick2", lwd=llw)
	
	
	text(xl1+11, yl2+.1,paste("R2=",round(summary(fit3)$r.squared,4)),cex=tx,col="firebrick2")
	text(xl1+11, yl2+.25,paste("p=",round(summary(fit3)$coefficients[2,4],4)),cex=tx,col="firebrick2")
	
	
	axis(1, xS1, rep("", length(xS1)), lwd.ticks=ltw)
	mtext(xS1, at=xS1, side=1, cex=lmx, line=xll)	
	
	axis(2, yS2, rep("", length(yS2)), lwd.ticks=ltw)
	mtext(yS2, at=yS2, side=2, cex=lmx,las=2, line=yll)
	mtext("swe melt midpoint",side=2,line=40,cex=10)
	mtext("(mid)",side=2,line=25,cex=10)
	mtext("Ave melt temp (C)",side=1,line=20,cex=10)
	box(which="plot")	
	
						
	par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl2,xh2), ylim=c(yl2,yh2), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")				
					
		points(mid0glc[[i]]$vcf,mid0glc[[i]]$Mean, pch=19, cex=px,col=rgb(0,0,0,.15))	
		abline(fit4, col="firebrick2", lwd=llw)	
		
	text(xl2+31, yl2+.1,paste("R2=",round(summary(fit4)$r.squared,4)),cex=tx,col="firebrick2")
	text(xl2+31, yl2+.25,paste("p=",round(summary(fit4)$coefficients[2,4],4)),cex=tx,col="firebrick2")
	axis(1, xS2, rep("", length(xS2)), lwd.ticks=ltw)
	mtext(xS2, at=xS2, side=1, cex=lmx, line=xll)
	mtext("Canopy cover (%)",side=1,line=20,cex=10)	
	box(which="plot")
	
	par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl3,xh3), ylim=c(yl2,yh2), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")				
					
		points(mid0glc[[i]]$max.swe,mid0glc[[i]]$Mean, pch=19, cex=px,col=rgb(0,0,0,.15))		
	abline(fit6, col="firebrick2", lwd=llw)	
	text(xl3+.21, yl2+.1,paste("R2=",round(summary(fit6)$r.squared,4)),cex=tx,col="firebrick2")
	text(xl3+.21, yl2+.25,paste("p=",round(summary(fit6)$coefficients[2,4],4)),cex=tx,col="firebrick2")
	
	axis(4, yS2, rep("", length(yS2)), lwd.ticks=ltw)
	mtext(yS2, at=yS2, side=4, cex=lmx,las=2, line=yll)

	axis(1, xS3, rep("", length(xS3)), lwd.ticks=ltw)
	mtext(xS3, at=xS3, side=1, cex=lmx, line=xll)
	box(which="plot")
	
	mtext("swe melt midpoint",side=4,line=25,cex=10)
	mtext("(mid)",side=4,line=40,cex=10)
	mtext("Max swe ",side=1,line=20,cex=10)
	
	
	dev.off()
}					


############################################
###  example curve plot                  ###
############################################
#swe curve function
sweC <- function( b,day,mid){
	1/(1+exp(b*(day-mid)))
}
xseq <- seq(0,1, by=.1)
vegsub <- 1
pixsub <- 1
cellSub <- b0glc[[vegsub]]$cell[which(b0glc[[vegsub]]$pixID==pixsub)]
jpeg(paste0(plotDI,"\\example_curve.jpeg"), width=700,height=700,quality=100)
par(mai=c(1,1,1,1))
	plot(xseq,sweC( b0glc[[vegsub]]$Mean[pixsub]	,xseq, mid0glc[[vegsub]]$Mean[pixsub]),
			type="l", col="firebrick2", lwd=3, xlab="Proportion into melt period", ylab="Proportion of swe peak", cex=2,
			cex.axis=2, cex.lab=2)
dev.off()
vegsub <- 1
pixsub <- 1
cellSub <- b0glc[[vegsub]]$cell[which(b0glc[[vegsub]]$pixID==pixsub)]
jpeg(paste0(plotDI,"\\problem_curve.jpeg"), width=700,height=700,quality=100)
par(mai=c(1,1,1,1))
	plot(dat.swe6$jday[dat.swe6$gcID==vegsub&dat.swe6$cell==cellSub]-32/(182-32),dat.swe6$sweN[dat.swe6$gcID==vegsub&dat.swe6$cell==cellSub],
			pch=19,   xlab="Proportion into melt period", ylab="Proportion of swe peak", 
			cex.axis=2, cex.lab=2)
dev.off()
	

plot(b0glc[[1]]$Mean,mid0glc[[1]]$Mean)	
############################################
###  plot params                         ###
############################################
#make a simplier glc plotting name

					
		
					
					
#normalized day
dat.swe6$dayN <- (dat.swe5$jday-32)/(182-32)
#swe curve function
sweC <- function( b,day,mid){
	1/(1+exp(b*(day-mid)))
}
#plotting info					
wd <- 50
hd <- 50

xl <- 0
xh <- 1
yl <- 0
yh <- 0.5
xplot <- seq(0,1,by=.01)

xL <- seq(32,182, by=10)
xS <- seq(0,1, length.out=length(xL))
yS <- seq(0,0.5,by=.05)
#dim(IDSglc)[1]
for(i in 1:1){
	jpeg(paste0(plotDI,"\\curves_",IDSglc$name[i],".jpeg"), width=2000,height=2000,quality=100)
	layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
	par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl,xh), ylim=c(yl,yh), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")
				
		points(dat.swe5$dayN,dat.swe5$swe,pch=19, col="grey75")

		points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i],
							datC$Mean[rownames(datC)=="b01"][i],
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i],
							datC$Mean[rownames(datC)=="base01"][i]),type="l", lwd=6,
							col="goldenrod3")
							
							
		points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i]+(datC$Mean[rownames(datC)=="M11"][i]*-1),
							datC$Mean[rownames(datC)=="b01"][i]+(datC$Mean[rownames(datC)=="b11"][i]*-1),
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i]+(datC$Mean[rownames(datC)=="mid11"][i]*-1),
							datC$Mean[rownames(datC)=="base01"][i]+(datC$Mean[rownames(datC)=="base11"][i]*-1)),type="l", lwd=6,
							col="royalblue")	

		points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i]+(datC$Mean[rownames(datC)=="M11"][i]*-12),
							datC$Mean[rownames(datC)=="b01"][i]+(datC$Mean[rownames(datC)=="b11"][i]*-12),
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i]+(datC$Mean[rownames(datC)=="mid11"][i]*-12),
							datC$Mean[rownames(datC)=="base01"][i]+(datC$Mean[rownames(datC)=="base11"][i]*-12)),type="l", lwd=6,
							col="royalblue4")	
		points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i]+(datC$Mean[rownames(datC)=="M11"][i]*1),
							datC$Mean[rownames(datC)=="b01"][i]+(datC$Mean[rownames(datC)=="b11"][i]*1),
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i]+(datC$Mean[rownames(datC)=="mid11"][i]*1),
							datC$Mean[rownames(datC)=="base01"][i]+(datC$Mean[rownames(datC)=="base11"][i]*1)),type="l", lwd=6,
							col="tomato3")	

			points(xplot, sweC(datC$Mean[rownames(datC)=="M01"][i]+(datC$Mean[rownames(datC)=="M11"][i]*14),
							datC$Mean[rownames(datC)=="b01"][i]+(datC$Mean[rownames(datC)=="b11"][i]*14),
							xplot,
							datC$Mean[rownames(datC)=="mid01"][i]+(datC$Mean[rownames(datC)=="mid11"][i]*14),
							datC$Mean[rownames(datC)=="base01"][i]+(datC$Mean[rownames(datC)=="base11"][i]*14)),type="l", lwd=6,
							col="tomato4")							
		axis(1, xS,rep("", length(xS)), lwd.ticks=5)
		mtext(xL,at=xS, side=1, line=5, cex=2)
		axis(2,yS, rep("", length(yS)), lwd.ticks=5)
		mtext(yS,at=yS,side=2, line=5, cex=2, las=2)
	dev.off()
}					
					

hist(temp.py$tempCent[temp.py$zone==5])


############################################
###  plot data                           ###
############################################

##plot maximum swe value at beginning of season before melt
#get the maximum swe value by averaging well before melt

dat.sweMax <- dat.swe5[dat.swe5$jday<90,]
#aggregate by cell, year
sweMax <- aggregate(dat.sweMax$swe, by=list(dat.sweMax$cell,dat.sweMax$year,dat.sweMax$zone), FUN="mean")
colnames(sweMax) <- c("cell","year","zone","sweMax")

#join with temperature data
sweMax2 <- join(sweMax, temp.py, by=c("cell","year","zone"), type="left")
#join glc data
sweMax3 <- join(sweMax2, IDSglc, by=c("zone"),type="left")
# calculate min swe

dat.sweMin <- dat.swe5[dat.swe5$jday>175,]

#aggregate by cell, year
sweMin <- aggregate(dat.sweMin$swe, by=list(dat.sweMin$cell,dat.sweMin$year,dat.sweMin$zone), FUN="mean")
colnames(sweMin) <- c("cell","year","zone","sweMin")

#join with temperature data
sweMin2 <- join(sweMin, temp.py, by=c("cell","year","zone"), type="left")
#join glc data
sweMin3 <- join(sweMax2, IDSglc, by=c("zone"),type="left")
# calculate min swe


#get unique tree cover
treeCover <- unique(data.frame(vcf=dat.swe5$vcf,cell=dat.swe5$cell))
length(unique(dat.swe5$cell))
dim(treeCover)[1]

#join unique tree cover into min and max
sweMin3 <- join(sweMin3,treeCover, by="cell",type="left")
sweMax3 <- join(sweMax3,treeCover, by="cell",type="left")

###### plot swe max ###### 

#plotting info					
wd <- 50
hd <- 50

xl <- -20
xh <- 20
yl <- 0
yh <- 0.5
yearU <- list()
yearN <- numeric(0)
xS <- seq(xl,xh,by=10)
yS <- seq(0,.5, by=.1)



for(i in 1:dim(IDSglc)[1]){
	yearU[[i]] <- unique(sweMax3$year[sweMax3$gcID==i])
	yearN[i] <- length(yearU[[i]])
	jpeg(paste0(plotDI,"\\new_data\\maxSWE\\maxSWE_",IDSglc$name[i],".jpeg"), width=2000,height=2000,quality=100)
	layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
		par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl,xh), ylim=c(yl,yh), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")
		for(j in 1:yearN[i]){		
			points(sweMax3$tempCent[sweMax3$gcID==i&sweMax3$year==yearU[[i]][j]],sweMax3$sweMax[sweMax3$gcID==i&sweMax3$year==yearU[[i]][j]],pch=19, col=j, cex=3)
		}
		axis(1, xS, rep("",length(xS)), lwd.ticks=4)
		mtext(xS, at=xS, cex=4,side=1,line=4)
		axis(2, yS, rep("",length(yS)), lwd.ticks=4)
		mtext(yS, at=yS, cex=4,side=2,line=4,las=2)
		mtext("Average spring temperature - long spring temperature average in glc (C)", side=1, cex=5, line=8)
		mtext("SWE", side=2, cex=5, line=8)
	dev.off()
}	





#plotting info					
wd <- 50
hd <- 50

xl <- -20
xh <- 20
yl <- 0
yh <- 0.5
yearU <- list()
yearN <- numeric(0)
xS <- seq(xl,xh,by=10)
yS <- seq(0,.5, by=.1)



for(i in 1:dim(IDSglc)[1]){
	yearU[[i]] <- unique(sweMin3$year[sweMin3$gcID==i])
	yearN[i] <- length(yearU[[i]])
	jpeg(paste0(plotDI,"\\new_data\\minSWE\\minSWE_",IDSglc$name[i],".jpeg"), width=2000,height=2000,quality=100)
	layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
		par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl,xh), ylim=c(yl,yh), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")
		for(j in 1:yearN[i]){		
			points(sweMin3$tempCent[sweMin3$gcID==i&sweMax3$year==yearU[[i]][j]],sweMin3$sweMax[sweMin3$gcID==i&sweMin3$year==yearU[[i]][j]],pch=19, col=j, cex=3)
		}
		axis(1, xS, rep("",length(xS)), lwd.ticks=4)
		mtext(xS, at=xS, cex=4,side=1,line=4)
		axis(2, yS, rep("",length(yS)), lwd.ticks=4)
		mtext(yS, at=yS, cex=4,side=2,line=4,las=2)
		mtext("Average spring temperature - long spring temperature average in glc (C)", side=1, cex=5, line=8)
		mtext("end of spring SWE", side=2, cex=5, line=10)
	dev.off()
}	
	
###### plot swe max and tree cover ###### 

#plotting info					
wd <- 50
hd <- 50

xl <- 0
xh <- 100
yl <- 0
yh <- 0.5
yearU <- list()
yearN <- numeric(0)
xS <- seq(xl,xh,by=10)
yS <- seq(0,.5, by=.1)



for(i in 1:dim(IDSglc)[1]){
	yearU[[i]] <- unique(sweMax3$year[sweMax3$gcID==i])
	yearN[i] <- length(yearU[[i]])
	jpeg(paste0(plotDI,"\\new_data\\maxSWEc\\maxSWE_",IDSglc$name[i],".jpeg"), width=2000,height=2000,quality=100)
	layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
		par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl,xh), ylim=c(yl,yh), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")
		for(j in 1:yearN[i]){		
			points(sweMax3$vcf[sweMax3$gcID==i&sweMax3$year==yearU[[i]][j]],sweMax3$sweMax[sweMax3$gcID==i&sweMax3$year==yearU[[i]][j]],pch=19, col=j, cex=3)
		}
		axis(1, xS, rep("",length(xS)), lwd.ticks=4)
		mtext(xS, at=xS, cex=4,side=1,line=4)
		axis(2, yS, rep("",length(yS)), lwd.ticks=4)
		mtext(yS, at=yS, cex=4,side=2,line=4,las=2)
		mtext("Tree cover", side=1, cex=5, line=8)
		mtext("SWE", side=2, cex=5, line=8)
	dev.off()
}	





#plotting info					
wd <- 50
hd <- 50

xl <-0
xh <- 100
yl <- 0
yh <- 0.5
yearU <- list()
yearN <- numeric(0)
xS <- seq(xl,xh,by=10)
yS <- seq(0,.5, by=.1)

#plot min and tree cover

for(i in 1:dim(IDSglc)[1]){
	yearU[[i]] <- unique(sweMin3$year[sweMin3$gcID==i])
	yearN[i] <- length(yearU[[i]])
	jpeg(paste0(plotDI,"\\new_data\\minSWEc\\minSWE_",IDSglc$name[i],".jpeg"), width=2000,height=2000,quality=100)
	layout(matrix(c(1),ncol=1),width=lcm(wd),height=lcm(hd))
		par(mai=c(0,0,0,0))
		plot(c(0,1),c(0,1), type="n", xlim=c(xl,xh), ylim=c(yl,yh), axes=FALSE,
				xaxs="i", yaxs="i", xlab=" ", ylab=" ")
		for(j in 1:yearN[i]){		
			points(sweMin3$vcf[sweMin3$gcID==i&sweMax3$year==yearU[[i]][j]],sweMin3$sweMax[sweMin3$gcID==i&sweMin3$year==yearU[[i]][j]],pch=19, col=j, cex=3)
		}
		axis(1, xS, rep("",length(xS)), lwd.ticks=4)
		mtext(xS, at=xS, cex=4,side=1,line=4)
		axis(2, yS, rep("",length(yS)), lwd.ticks=4)
		mtext(yS, at=yS, cex=4,side=2,line=4,las=2)
		mtext("tree cover", side=1, cex=5, line=8)
		mtext("end of spring SWE", side=2, cex=5, line=10)
	dev.off()
}		

dfA <- unique(data.frame(cell=dat.swe4$cell, year=dat.swe4$year))
dfC <- unique(data.frame(cell=dat.swe4$cell))
dfY <- unique(data.frame(year=dat.swe4$year))