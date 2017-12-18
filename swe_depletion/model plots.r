library(coda)
library(plyr)

#start with a test year
dat.swe <- read.csv("c:\\Users\\hkropp\\Google Drive\\swe_test\\glc_SWE_M2BGS_ease2_2010.csv")
dat.gl <- read.csv("c:\\Users\\hkropp\\Google Drive\\swe_test\\glc50_table.csv")

#subset swe data 
#if no glc remove 
#add a cellid
datC <- data.frame(cellID= seq(1, dim(dat.swe)[1]))

dat.swe <- cbind(datC, dat.swe)
sweS <- dat.swe[!is.na(dat.swe$glc),]
#now subset to look at only the 3 months of snow melt
#start Mar 1 end May 31 doy (leap year) 61:152
#+2 for glc column and cellid
sweS2 <- sweS[,63:154]
#join ids back into swe
sweS3 <- cbind(sweS[,1:2],sweS2)
#take out any areas with no snow extent
#over the period
sweSums <- rowSums(sweS3[,3:94], na.rm=TRUE)
sweS4 <- sweS3[sweSums!=0,]

#see how many cells there are in each glc
gridC <- aggregate(sweS4[,1], by=list(sweS4[,2]), FUN="length")
colnames(gridC) <- c("glc", "count")

swe.vector <-as.vector(t(data.matrix(sweS4[,3:94])))
#turn into a data frame
sweDF <- data.frame(glc=rep(sweS4[,2],each=92),
			doy=rep(seq(61,152),times=dim(sweS4)[1]),
			gridID= rep(sweS4[,1],each=92),
			swe=swe.vector)

#plot the mean by vege type to see
sweM <- aggregate(sweDF$swe, by=list(sweDF$doy,sweDF$glc),FUN="mean", na.rm=TRUE)
colnames(sweM) <- c("doy","glc","swe")			
#see how many observations are in
plot(sweM$doy[sweM$glc==2],sweM$swe[sweM$glc==2])
plot(sweM$doy,sweM$swe)

#now subset to only include landcover classes intrested in
datgcI <- dat.gl[,1:2]
colnames(datgcI) <- c("glc","name")
gcInfo <- join(datgcI, gridC, by=c("glc"), type="left")
#don't include glc 3,10,15,16,18,20,21,22
gcUse <- gcInfo[gcInfo$glc!=2&gcInfo$glc!=3&gcInfo$glc!=15&gcInfo$glc!=16&gcInfo$glc!=10&gcInfo$glc!=14&
				gcInfo$glc!=9&gcInfo$glc!=14&gcInfo$glc<18,]
gcUse$gcID <- seq(1,  dim(gcUse)[1])
#subset to include only relevant land cover classes
sweDF2 <- join(sweDF, gcUse, by="glc", type="inner")

plot(sweDF2$doy[sweDF2$gcID==7],sweDF2$swe[sweDF2$gcID==7])
				
#get unique gridcells
gridID <- unique(data.frame(gridID=sweDF2$gridID, gcID=sweDF2$gcID))	
gridID$MgID <- seq(1, dim(gridID)[1])
#join back into
sweDF3 <- join(sweDF2, gridID, by=c("gridID","gcID"), type="left")

####################################
## read in output from the model ###
####################################
modDir <- "c:\\Users\\hkropp\\Google Drive\\swe_test\\stan1\\stan_all1"

MD <- list.files(modDir)
MSP1 <- strsplit(MD, "out")
MSP2 <- strsplit(MD, "gc")
MSP3 <-strsplit(MD, "chain")
#read in output
CD <- list()
Cparm <- character(0)
gN <- character(0)
ch1 <- character(0)
for(i in 1:length(MD)){
	CD[[i]] <- read.csv(paste0(modDir,"\\",MD[i]))
	Cparm[i] <- MSP1[[i]][1]
	gN[i]<- MSP2[[i]][2]
	ch1[i] <- MSP3[[i]][2]
}
Mch2<-strsplit(ch1, "gc")
ch2<- character(0)
for(i in 1:length(MD)){
	ch2[i]<- Mch2[[i]][1]
}


#get glcID
glcM1 <- gsub(".csv","",gN)
glcM <- as.numeric(gsub("_","",glcM1))
#parm
parm <- gsub("_","",Cparm)	
#chain
chain <- as.numeric(gsub("_","",ch2))
#give colnames
CD2 <- list()
for(i in 1:length(MD)){
	colnames(CD[[i]])<- paste0(parm[i],"[",glcM[i],"]")
	CD2[[i]]<- mcmc(CD[[i]])
}




#get unique parm and glc combinations
parmALL <- data.frame(glc=glcM, parm=parm,chain=chain)
parmALL$MID <- seq(1, dim(parmALL)[1])

parmIDS <- unique(data.frame(glc=glcM, parm=parm))
parmSubID <- list()
for(i in 1:dim(parmIDS)[1]){
	parmSubID[[i]] <- parmALL$MID[parmALL$glc==parmIDS$glc[i]&parmALL$parm==parmIDS$parm[i]]
}
#put together into mcmc list
CODA_all<- list()
for(i in 1:dim(parmIDS)[1]){
	CODA_all[[i]]<- mcmc.list(CD2[[parmSubID[[i]][1]]],CD2[[parmSubID[[i]][2]]],CD2[[parmSubID[[i]][3]]])
}


#make history plots
for(i in 1:dim(parmIDS)[1]){
	jpeg(paste0("c:\\Users\\hkropp\\Google Drive\\swe_test\\stan1\\history\\",parmIDS$parm[i],parmIDS$glc[i],".jpg"),
			width=2000, height=1000, units="px", quality=100)
	plot(CODA_all[[i]])		
	dev.off()
}


mod.outS<- numeric(0)
mod.outQ1<- numeric(0)
mod.outQ2<- numeric(0)
for(i in 1:dim(parmIDS)[1]){
	mod.outS[i]<-summary(CODA_all[[i]])$statistics[1]
	mod.outQ1[i]<-summary(CODA_all[[i]])$quantiles[1]
	mod.outQ2[i]<-summary(CODA_all[[i]])$quantiles[5]
}

parmsOut <- data.frame(parmIDS, Mean=mod.outS,pc2.5=mod.outQ1,pc97.5=mod.outQ2 )



#swe depletion curve
depletion<- function(b,day,mid,M,base){
	(M/(1+exp(b*(day-mid))))+base
}


#make a plot of the swe depletion curves
sweMuse <- join(sweM,gcUse, by="glc",type="inner")

colorTable <- data.frame(gcID=gcUse$gcID,
				coli=c("darkgreen","goldenrod3",
						"grey50","royalblue3",
						"plum4","tan4"))
						
sweMuse<- join(sweMuse,colorTable, by="gcID", type="left")						
Mid <-(61+((152-61)/2))
wd <- 45
hd <- 45
yl <- 0
yh <- 150
xlD <- 60
xhD <- 153
xseq <- seq(60,150, by=10)
yseq <- seq(0,150, by=30)
jpeg("c:\\Users\\hkropp\\Google Drive\\swe_test\\stan1\\plots\\swe_deplete.jpg", width=3200, height=1500, units="px", quality=100)
ab <- layout(matrix(c(1), ncol=1, byrow=FALSE), width=rep(lcm(wd),1), height=rep(lcm(hd),1))
	par(mai=c(0,0,0,0))
	plot(c(0,1),c(0,1),type="n", ylim=c(yl,yh), xlim=c(xlD-1,xhD), xlab=" ", ylab=" ", axes=FALSE, xaxs="i",
			yaxs="i")
	points(sweMuse$doy, sweMuse$swe, pch=19, col=paste(sweMuse$coli), cex=4)		
	for(i in 1:dim(gcUse)[1]){
	points(seq(xlD,xhD,by=.1),
		depletion(parmsOut$Mean[parmsOut$parm=="b"&parmsOut$glc==i],
		seq(xlD,xhD,by=.1),Mid,
		parmsOut$Mean[parmsOut$parm=="M"&parmsOut$glc==i],
		parmsOut$Mean[parmsOut$parm=="base"&parmsOut$glc==i]),
		type="l", lwd=5, col=paste(colorTable$coli[i]))
	}
	for(i in 1:dim(gcUse)[1]){
	polygon(c(seq(xlD,xhD,by=.1),rev(seq(xlD,xhD,by=.1))),
	c(depletion(parmsOut$pc2.5[parmsOut$parm=="b"&parmsOut$glc==i],
		seq(xlD,xhD,by=.1),Mid,
		parmsOut$pc2.5[parmsOut$parm=="M"&parmsOut$glc==i],
		parmsOut$pc2.5[parmsOut$parm=="base"&parmsOut$glc==i]),
	rev(depletion(parmsOut$pc97.5[parmsOut$parm=="b"&parmsOut$glc==i],
		seq(xlD,xhD,by=.1),Mid,
		parmsOut$pc97.5[parmsOut$parm=="M"&parmsOut$glc==i],
		parmsOut$pc97.5[parmsOut$parm=="base"&parmsOut$glc==i])))
		
		)	
	}
	axis(1, xseq, cex.axis=3,lwd.ticks=5)
	axis(2, yseq, cex.axis=3,lwd.ticks=5, las=2)
dev.off()	

par(mai=c(1,1,1,1))
plot(sweDF3$doy[sweDF3$gcID==1],sweDF3$swe[sweDF3$gcID==1],pch=19,
		axes=FALSE, xlab="Day of Year",ylab="SWE", main="Evergreen Needleleaf",
		cex.lab=2, ylim=c(0,470))
	points(seq(xlD,xhD,by=.1),
		depletion(parmsOut$Mean[parmsOut$parm=="b"&parmsOut$glc==1],
		seq(xlD,xhD,by=.1),Mid,
		parmsOut$Mean[parmsOut$parm=="M"&parmsOut$glc==1],
		parmsOut$Mean[parmsOut$parm=="base"&parmsOut$glc==1]),
		type="l", lwd=5, col="tomato3")
	points(sweMuse$doy[sweMuse$gcID==1], sweMuse$swe[sweMuse$gcID==1], pch=19, col="royalblue4")	
	axis(1,seq(60,150, by=10), cex.axis=1.25)
	axis(2,seq(0,450, by=50), cex.axis=1.25, las=2)
	legend(70,470, c("SWE observations","observation mean","model mean"),
		col=c("black","royalblue3","tomato3"),pch=c(19,19,NA),lwd=c(NA,NA,4),
		bty="n")
	
par(mai=c(1,1,1,1))
plot(sweDF3$doy[sweDF3$gcID==2],sweDF3$swe[sweDF3$gcID==2],pch=19,
		axes=FALSE, xlab="Day of Year",ylab="SWE", main="Deciduous Needleleaf", 
		cex.lab=2, ylim=c(0,470))
	points(seq(xlD,xhD,by=.1),
		depletion(parmsOut$Mean[parmsOut$parm=="b"&parmsOut$glc==2],
		seq(xlD,xhD,by=.1),Mid,
		parmsOut$Mean[parmsOut$parm=="M"&parmsOut$glc==2],
		parmsOut$Mean[parmsOut$parm=="base"&parmsOut$glc==2]),
		type="l", lwd=5, col="tomato3")
	points(sweMuse$doy[sweMuse$gcID==2], sweMuse$swe[sweMuse$gcID==2], pch=19, col="royalblue4")	
	axis(1,seq(60,150, by=10), cex.axis=1.25)
	axis(2,seq(0,450, by=50), cex.axis=1.25, las=2)	
	legend(70,470, c("SWE observations","observation mean","model mean"),
		col=c("black","royalblue3","tomato3"),pch=c(19,19,NA),lwd=c(NA,NA,4),
		bty="n")	

	
par(mai=c(1,1,1,1))
plot(sweDF3$doy[sweDF3$gcID==6],sweDF3$swe[sweDF3$gcID==6],pch=19,
		axes=FALSE, xlab="Day of Year",ylab="SWE", main="Herbaceous", 
		cex.lab=2, ylim=c(0,470))
	points(seq(xlD,xhD,by=.1),
		depletion(parmsOut$Mean[parmsOut$parm=="b"&parmsOut$glc==6],
		seq(xlD,xhD,by=.1),Mid,
		parmsOut$Mean[parmsOut$parm=="M"&parmsOut$glc==6],
		parmsOut$Mean[parmsOut$parm=="base"&parmsOut$glc==6]),
		type="l", lwd=5, col="tomato3")
	points(sweMuse$doy[sweMuse$gcID==6], sweMuse$swe[sweMuse$gcID==6], pch=19, col="royalblue4")	
	axis(1,seq(60,150, by=10), cex.axis=1.25)
	axis(2,seq(0,450, by=50), cex.axis=1.25, las=2)	
	legend(70,470, c("SWE observations","observation mean","model mean"),
		col=c("black","royalblue3","tomato3"),pch=c(19,19,NA),lwd=c(NA,NA,4),
		bty="n")			
		legend(70,470, c("SWE observations","observation mean","model mean"),
		col=c("black","royalblue3","tomato3"),pch=c(19,19,NA),lwd=c(NA,NA,4),
		bty="n")
	
par(mai=c(1,1,1,1))
plot(sweDF3$doy[sweDF3$gcID==5],sweDF3$swe[sweDF3$gcID==5],pch=19,
		axes=FALSE, xlab="Day of Year",ylab="SWE", main="Deciduous shrub", 
		cex.lab=2, ylim=c(0,300))
	points(seq(xlD,xhD,by=.1),
		depletion(parmsOut$Mean[parmsOut$parm=="b"&parmsOut$glc==5],
		seq(xlD,xhD,by=.1),Mid,
		parmsOut$Mean[parmsOut$parm=="M"&parmsOut$glc==5],
		parmsOut$Mean[parmsOut$parm=="base"&parmsOut$glc==5]),
		type="l", lwd=5, col="tomato3")
	points(sweMuse$doy[sweMuse$gcID==5], sweMuse$swe[sweMuse$gcID==5], pch=19, col="royalblue4")	
	axis(1,seq(60,150, by=10), cex.axis=1.25)
	axis(2,seq(0,450, by=50), cex.axis=1.25, las=2)	
	legend(70,300, c("SWE observations","observation mean","model mean"),
		col=c("black","royalblue3","tomato3"),pch=c(19,19,NA),lwd=c(NA,NA,4),
		bty="n")			
	
#depletion example
#338	
gridID[gridID$gcID==1,]
par(mai=c(1,1,1,1))
plot(sweDF3$doy[sweDF3$MgID==616],sweDF3$swe[sweDF3$MgID==616],pch=19,
		axes=FALSE, xlab="Day of Year",ylab="SWE", 
		main="Evergreen needleleaf gridcell", 
		cex.lab=2, ylim=c(0,150), xlim=c(60,120))	
	axis(1,seq(60,120, by=10), cex.axis=1.25)
	axis(2,seq(0,450, by=50), cex.axis=1.25, las=2)
points(seq(60,120,by=.1),
		depletion(.5,seq(60,120,by=.1),90,110,0),
		type="l", lwd=5, col="tomato3")	
	
	#polygon(c(seq(xlD,xhD,by=.1),rev(seq(xlD,xhD,by=.1))),
	#c(depletion(parmsOut$pc2.5[parmsOut$parm=="b"&parmsOut$glc==1],
	#	seq(xlD,xhD,by=.1),Mid,
	#	parmsOut$pc2.5[parmsOut$parm=="M"&parmsOut$glc==1],
	#	parmsOut$pc2.5[parmsOut$parm=="base"&parmsOut$glc==1]),
	#rev(depletion(parmsOut$pc97.5[parmsOut$parm=="b"&parmsOut$glc==1],
	#	seq(xlD,xhD,by=.1),Mid,
	#	parmsOut$pc97.5[parmsOut$parm=="M"&parmsOut$glc==1],
	#	parmsOut$pc97.5[parmsOut$parm=="base"&parmsOut$glc==1]))),
	#col=rgb(205/255,79/255,57/255,.5),border=NULL)