# Delete all objects in the work space
rm(list=ls())

####################################################
##              select clim parameters            ##
####################################################

select <- c("DC7", "Ptot6", "Sftot3", "Tmean-9", "Sftot-9", "GSL-6")
# noms qui apparaitront dans le modèle final
namesclimvar <- c("DC7", "P6", "S3", "Tp9", "Sp9", "GSLp")

####################################################
##                   climatic data                ##
####################################################

# Import Temperature / Precipitation / Snowfall
clim<-read.csv("~/ownCloud/Scan & measures/Climate/PTSFBIC.csv", sep=",")
clim[,"year"]<-as.integer(substr(clim[,1],1,4))
clim[,"month"]<-as.integer(substr(clim[,1],6,7))
clim<-clim[,-1]
colnames(clim)<-c("Tmean","Ptot", "Sftot", "year","month")
clim<-data.frame(clim$year, clim$month, clim$Tmean, clim$Ptot, clim$Sftot)
colnames(clim)<-c("year","month","Tmean","Ptot", "Sftot")

# Import Degree Day
clim1<-read.csv("~/ownCloud/Scan & measures/Climate/DDBIC.csv", sep=",")
clim1[,"year"]<-as.integer(substr(clim1[,1],1,4))
clim1[,"month"]<-as.integer(substr(clim1[,1],6,7))
clim1<-clim1[,-1]
colnames(clim1)<-c("DD","year","month")
clim1<-data.frame(clim1$year, clim1$month, clim1$DD)
colnames(clim1)<-c("year","month","DD")

# Merge P T SF & DD
clim <- cbind(clim, clim1[,3])
colnames(clim)<-c("year","month","Tmean","Ptot", "Sftot", "DD")

# Import Drought Code
clim1<-read.csv("~/ownCloud/Scan & measures/Climate/DCBIC.csv", sep=",")
clim1[,"year"]<-as.integer(substr(clim1[,1],1,4))
clim1[,"month"]<-as.integer(substr(clim1[,1],6,7))
clim1<-clim1[,-1]
clim1<-data.frame(clim1$year, clim1$month, clim1$DC)
colnames(clim1)<-c("year","month","DC")

# Merge DC to the other variables
clim[,"both"] <- paste(clim[,"year"], clim[,"month"], sep="")
clim1[,"both"] <- paste(clim1[,"year"], clim1[,"month"], sep="")
clim <- merge(clim, clim1, by = "both", all = T)
attach(clim)
clim <- clim[order(year.x, month.x),]
detach(clim)
clim <- clim[, -c(1, 8, 9)]
colnames(clim)<-c("year","month","Tmean","Ptot", "Sftot", "DD", "DC")

# Import Growth Season Length
clim1<-read.csv("~/ownCloud/Scan & measures/Climate/GSLBIC.csv", sep=",")
clim1[,"month"] <- 6
clim1 <- clim1[, -c(2,3)]
colnames(clim1)<-c("year","GSL", "month")
clim1<-data.frame(clim1$year, clim1$month, clim1$GSL)
colnames(clim1)<-c("year", "month", "GSL")

# Merge GSL to the other variables
clim[,"both"] <- paste(clim[,"year"], clim[,"month"], sep="")
clim1[,"both"] <- paste(clim1[,"year"], clim1[,"month"], sep="")
clim <- merge(clim, clim1, by = "both", all= T)
attach(clim)
clim <- clim[order(year.x, month.x),]
detach(clim)
clim <- clim[, -c(1, 9, 10)]
colnames(clim)<-c("year","month","Tmean","Ptot", "Sftot", "DD", "DC", "GSL")


# change climat dataset format
nbclimvar <- dim(clim)[2]-2     # number of climatic variables
clim1 <- data.frame(row.names = levels(as.factor(clim$year)))

for (i in 1:nbclimvar){     # i = select the climatic variable in the first data frame
  temp <- data.frame(row.names = levels(as.factor(clim$year)))
  for (j in 1:12){          # j = month
    temp2 <- as.data.frame(clim[clim$month == j,i+2])
    temp <- cbind(temp,temp2[,1])
    colnames(temp)[j] <- paste(colnames(clim)[i+2],j, sep="")  # column title
  }
  clim1 <- cbind(clim1, temp)
}
clim <- clim1
clim <- as.matrix(clim)

# Take into acount month [-6:9] = previous June to current October
temp <- clim[,c(6:12, 18:24, 30:36, 42:48, 54:60, 66:72)]
rownames(temp) <- as.character(as.numeric(rownames(temp))+1)
temp <- temp[-nrow(temp),]
colnames(temp) <- c("Ptot-6","Ptot-7","Ptot-8","Ptot-9","Ptot-10","Ptot-11","Ptot-12","Tmean-6","Tmean-7","Tmean-8","Tmean-9","Tmean-10","Tmean-11","Tmean-12","Sftot-6","Sftot-7","Sftot-8","Sftot-9","Sftot-10","Sftot-11","Sftot-12","DD-6","DD-7","DD-8","DD-9","DD-10","DD-11","DD-12","DC-6","DC-7","DC-8","DC-9","DC-10","DC-11","DC-12","GSL-6","GSL-7","GSL-8","GSL-9","GSL-10","GSL-11","GSL-12")
clim <- clim[-1,]
clim <- cbind(clim, temp)


####################################################
##              select clim parameters            ##
####################################################

clim <- as.data.frame(clim[,select])
colnames(clim) <- namesclimvar

####################################################
##                  save               ##
####################################################

write.csv(clim, "~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/input_clim/climBIC.csv")
