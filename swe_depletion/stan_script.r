library(coda)
library(plyr)
library(rstan)


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
gcUse <- gcInfo[gcInfo$glc!=3&gcInfo$glc!=15&gcInfo$glc!=16&gcInfo$glc!=10&
				gcInfo$glc!=18&gcInfo$glc<20,]
gcUse$gcID <- seq(1,  dim(gcUse)[1])
#subset to include only relevant land cover classes
sweDF2 <- join(sweDF, gcUse, by="glc", type="inner")

plot(sweDF2$doy[sweDF2$gcID==7],sweDF2$swe[sweDF2$gcID==7])
				
#get unique gridcells
gridID <- unique(data.frame(gridID=sweDF2$gridID, gcID=sweDF2$gcID))	
gridID$MgID <- seq(1, dim(gridID)[1])
#join back into
sweDF3 <- join(sweDF2, gridID, by=c("gridID","gcID"), type="left")

#subset one each gridcells and write fro trial model
#MgID=1,106,58,76,217,82,35,30,14,8
sweS <-sweDF3[sweDF3$MgID==1|sweDF3$MgID==106|sweDF3$MgID==58|sweDF3$MgID==76
		|sweDF3$MgID==217|sweDF3$MgID==82|sweDF3$MgID==35|sweDF3$MgID==30|sweDF3$MgID==14
		|sweDF3$MgID==8,]

gridS <-gridID[gridID$MgID==1|gridID$MgID==106|gridID$MgID==58|gridID$MgID==76
		|gridID$MgID==217|gridID$MgID==82|gridID$MgID==35|gridID$MgID==30|gridID$MgID==14
		|gridID$MgID==8,]
gridS$gID <- seq(1, dim(gridS)[1])		

sweS2 <- join(sweS, gridS, by=c("gridID","gcID","MgID"), type="left")

write.table(sweS2,"c:\\Users\\hkropp\\Google Drive\\swe_test\\test_sub_swe.csv", 
			sep=",", row.names=FALSE)
write.table(gridS,"c:\\Users\\hkropp\\Google Drive\\swe_test\\test_sub_swe.csv", 
			sep=",", row.names=FALSE)

plot(sweS2$doy,sweS2$swe)			
#######################################################
# set up model run                                    #
#######################################################
datalist <- list(Nobs=dim(sweDF3)[1], swe=sweDF3$swe, vegeC=sweDF3$gcID,
			day=sweDF3$doy,  mid=61+((152-61)/2),gridC=sweDF3$MgID,
			NGridC=dim(gridID)[1], vege=gridID$gcID,
			NVeg=dim(gcUse)[1])
			
			
parms <- c("M","base","sig.swe","b","sig.swe","sig.bV",
			"mu.base","mu.b",
			"sig.M","sig.base","sig.b")		