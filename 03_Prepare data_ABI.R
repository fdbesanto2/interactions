
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

  clim=read.csv("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/input_clim/climABI.csv",sep=",",dec=".", stringsAsFactors=FALSE)
  data <- read.csv("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/input_chron/ArABItree.csv",sep=",",dec=".", stringsAsFactors=FALSE)

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
  # arbres avec données de croissance
  tmp <- data[!is.na(data$X1991) | !is.na(data$X2000) | !is.na(data$X2010),]
  # bande de X mètres
  a <- -0.100
  b <- 5398790
  b2 <- b - 20
  ## arbres qui sont dans la bande
  data$edge <- a*data$X+b2
  data$edge1 <- data$Y - data$edge

  ## supprimer leur données de croissance
  data[data$edge1>0, substr(colnames(data),1,2)=="X1"] <- NA
  data[data$edge1>0, substr(colnames(data),1,2)=="X2"] <- NA
  tmp1 <- data[!is.na(data$X1991) | !is.na(data$X2000) | !is.na(data$X2010),]


  ggplot(data=data, aes(X,Y))+
  geom_point(aes(size=DHP11), alpha=0.5)+
  scale_x_continuous(breaks=seq(614600,614900,20))+
  geom_point(aes(size=DHP11), color="blue", alpha=0.5, data=data[data$edge1>0,])+
  geom_point(aes(size=DHP11), color="red", alpha=0.5, data=tmp)+
  stat_function(fun=function(x)a*x+b, color="red")+
  stat_function(fun=function(x)a*x+b2, color="red")+
  geom_point(color="green", alpha=0.7, data=tmp1)+
  coord_fixed(ratio = 1)+ # même échelle pour les X et les Y
  ylim(5337095,5337330)

data <- data[,-c(ncol(data)-1, ncol(data))] ## supprimer colonnes "edge" & "edge1"


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
