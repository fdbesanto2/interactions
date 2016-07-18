
# Delete all objects in the work space
rm(list=ls(all=TRUE))
library(ggplot2)

####################################################
##                  set radius                    ##
####################################################

## définition du rayon (pas diamètre) du cercle dans lequel la compétition opère autour de l'arbre focal
buffer <- 2000 # (en cm)

####################################################
##                  Data & Packages               ##
####################################################

#setwd("~/Google Drive/Scan & measures")

clim=read.csv("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/input_clim/clim1823.csv",sep=",",dec=".", stringsAsFactors=FALSE)
data <- read.csv("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/input_chron/To1823tree.csv",sep=",",dec=".", stringsAsFactors=FALSE)

####################################################
####################################################

ggplot(data=data, aes(X,Y))+
geom_point(aes(size=DHP11), alpha=0.5)

data$DHP11 <- as.numeric(data$DHP11)
data <- data[!grepl("\\?",data$DHP11),]  ## supprimer les ? dans le fichier
data <- data[data$DHP11!="",]            ## enlever les arbres sans DBH
data <- data[!is.na(data$DHP11),]        ## enlever les arbres sans DBH
data <- data[!is.na(data$X),]        ## enlever les arbres sans coordonnées
data <- data[!is.na(data$Y),]        ## enlever les arbres sans coordonnées
data <- data[data$DHP11!=0,]             ## enlever les arbres DBH = 0

#data <- data[data$DHP11<500,]             ## enlever les DBH trop grands (erreurs???)
data <- data[data$DHP11>=100,]             ## enlever les DBH trop petits (<10)
ggplot(data=data[data$DHP11>0,], aes(DHP11))+
geom_histogram(bins=100)

####################################################
##           set tmporal window for clim          ##
####################################################
# set same temporal windows as the one used for the chronologies
beg <- min(as.numeric(substr(colnames(data)[substr(colnames(data),1,5)== "DHP19"],4,7)))
end <- max(as.numeric(substr(colnames(data)[substr(colnames(data),1,5)== "DHP20"],4,7)))
clim <- as.data.frame(clim[clim$X>=beg & clim$X<=end,])

####################################################
##              Distance Calculation              ##
####################################################
############### delete trees on the edges manually
data[data$X>80 | data$X<20,substr(colnames(data),1,2)=="X1"] <- NA
data[data$X>80 | data$X<20, substr(colnames(data),1,2)=="X2"] <- NA
data[data$Y>80 | data$Y<20, substr(colnames(data),1,2)=="X1"] <- NA
data[data$Y>80 | data$Y<20, substr(colnames(data),1,2)=="X2"] <- NA


tmp <- data[!is.na(data$X1991) | !is.na(data$X2000) | !is.na(data$X2010),]
ggplot(data=data, aes(X,Y))+
geom_point(aes(size=DHP11), alpha=0.5)+
geom_point(aes(size=DHP11), color="red", alpha=0.5, data=tmp)


######## matrice des distance
dist <- as.matrix(dist(data[,c("X","Y")], method="euclidean", diag=F,upper=T)) ## création de la matrice de distance
dist <- dist*100
hist(dist, breaks=1000)
## distance minimale entre les arbres: 10 cm (notamment pour les fourches
## qui ont les même coordonnées)
dist[dist<10] <- 10
diag(dist) <- 0

dist[dist>buffer] <- 0    ## mise a 0 des distance > buffer
dist_inv <- 1/dist        ## calcul de l'inverse de la distance (car NCI = DBH/distance)
dist_inv[dist_inv==Inf] <- 0  ## mise à 0 des distance = Inf (celles précédament mise à 0)
hist(dist_inv[dist_inv>0])
dbh <- as.matrix(data[rownames(dist_inv),"DHP11"])  ## dbh des arbres dans le même ordre que la matrice de distance

stockdbh <- dbh
stockdist_inv <- dist_inv



setwd("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/RUN_MODEL")
save(data, file="chrono.Rdata")
save(clim, file="clim.Rdata")
save(stockdist_inv, file="dist_inv.Rdata")
save(stockdbh, file="dbh.Rdata")
