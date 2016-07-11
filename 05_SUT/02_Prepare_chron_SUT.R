# Delete all objects in the work space
rm(list=ls())


####################################################
##          set species and temporal window       ##
####################################################

# Load Package
library(dplR)
library(ggplot2)
# Import growth
rw <- read.tucson("~/ownCloud/Scan & measures/RWL/AbSUTepsme.rwl")
# import diameter & coordinates & IDs
dc <- read.csv("~/ownCloud/Work_directory/Data/Sampling/Sutton/SUT_uptodate_June2016.csv", sep=";")
# replace commas by dots
dc$X_UTM <- as.numeric(gsub(",", ".", gsub("\\.", "", dc$X_UTM)))
dc$Y_UTM <- as.numeric(gsub(",", ".", gsub("\\.", "", dc$Y_UTM)))
dc$DHP_mm_ <- as.numeric(gsub(",", ".", gsub("\\.", "", dc$DHP_mm_)))

# set temporal window
beg <- 1991
end <- 2013
growth <- rw[rownames(rw)>=beg & rownames(rw)<=end ,]

####################################################
##                  growth data                   ##
####################################################
# changes individuals names in rw
#colnames(growth) <- substr(colnames(growth), 4,7)
colnames(growth) <- substr(colnames(growth), 3,7)
colnames(growth)[substr(colnames(growth), 1,1)=="b"] <- substr(colnames(growth), 2,5)

growth <- t(growth)
growth <- as.data.frame(growth)
growth[,"TAG"] <- rownames(growth)

# supprimer les arbres pour lesquels on a pas de données de croissance
# sur la temporal window
growth
growth$somme <- apply(growth[,1:(ncol(growth)-1)],1,sum, na.rm=T)
growth <- growth[growth$somme!=0,]
growth <- growth[,-ncol(growth)]

####################################################
##           diameter & coordinates data          ##
####################################################
# vérifier le nom des espèces
colnames(dc)[colnames(dc)=="Esp"] <- "Sp"
unique(dc$Sp)

colnames(dc)[colnames(dc)=="DHP_mm_"] <- "DHP11"
colnames(dc)[colnames(dc)=="X_UTM"] <- "X"
colnames(dc)[colnames(dc)=="Y_UTM"] <- "Y"
colnames(dc)[colnames(dc)=="Arbre"] <- "TAG"

# only living trees (no DBH measures for dead trees)
#tab <- dc[dc[,"MORT"] == "0",c("TAG", "Sp","DHP11", "X", "Y")]
tab <- dc[dc[,"Etat"] == "V",c("TAG", "Sp","DHP11", "X", "Y")]
tab$TAG <- as.factor(tab$TAG)

# réduire la zone d'étude
tab <- tab[!is.na(tab$X),] # on supprime d'abord les arbres sans coordonnées
tab <- tab[!is.na(tab$Y),]

ggplot(data=tab, aes(X,Y))+
geom_point(aes(size=DHP11), alpha=0.5)
tab <- tab[tab$Y>4997000 & tab$Y<4997500,]
ggplot(data=tab, aes(X,Y))+
geom_point(aes(size=DHP11), alpha=0.5)


####################################################
##                   gather data                  ##
####################################################
# diameter, coordinates, ID & growth for cored trees
treei <- merge(tab, growth, by = "TAG", all.y = T)

# nb d'arbres carottés:
nrow(growth)
# duplicated trees?
sum(duplicated(rownames(growth)))
# nb d'arbres après merge:
nrow(treei)
# arbres en double après merge
sum(duplicated(treei$TAG))
# vérifier le nom des espèces carottés
unique(treei$Sp)

# voire les arbres dupliqués (si même nom et même coordonnées alors
# = fourche, on les exclus par la suite)

pb_TAG <- treei[duplicated(treei$TAG),"TAG"]
treei[treei$TAG %in% pb_TAG,]


#tab[tab$TAG=="3174",]

# changes individuals names in tab
#tab[,"TAG"] <- as.character(tab[,"TAG"])
#for (i in 1:nrow(tab)){
#  if (nchar(as.character(tab[i,"TAG"])) == 2){
#    tab[i,"TAG"] <- paste("00", tab[i,"TAG"], sep="")} else if (nchar(as.character(tab[i,"TAG"] #)) == 3){
#      tab[i,"TAG"] <- paste("0", tab[i,"TAG"], sep="")}
#}

tree <- merge(tab, growth, by = "TAG", all.x = T)
# repérer les arbres en double (fouches = twins)
pb_TAG <- tree[duplicated(tree$TAG),"TAG"]
nbtwin <- dim(tree[tree$TAG %in% pb_TAG,])[1] ## nombre de twins
## changer leurs noms
tree$TAG <- as.character(tree$TAG)
tree[tree$TAG %in% pb_TAG,"TAG"] <- paste("twin", 1:nbtwin ,sep="_")
tree$TAG <- as.factor(tree$TAG)

# Supression des données de croissance pour les twins
# car on ne sait pas auquel des deux twin attribuer
# la croissance
tree[substr(tree$TAG, 1, 4)=="twin",]
tree[substr(tree$TAG, 1, 4)=="twin",6:ncol(tree)] <- NA
tree[substr(tree$TAG, 1, 4)=="twin",]

# Supression des données de croissance pour les ind.
# de la mauvaise espèce (les ind. sans espèces sont exclus
# automatiquement)
tree[tree$TAG==10094,]
tree[tree$TAG==10094,6:ncol(tree)] <- NA
tree[tree$TAG==10094,]


####################################################
##                    diameters                   ##
####################################################

### faire varier les diamètres des arbres carottés au cours de la période d'étude
coldiam <- paste("DHP",seq(beg,end,1), sep="")
diam <- matrix(nrow=nrow(tree), ncol=length(coldiam))
diam <- cbind(diam, as.character(tree$TAG))
colnames(diam)<- c(coldiam,"TAG")

tree <- merge(tree, diam, by = "TAG")
tree$DHP2011 <- tree$DHP11

for (i in 2012:end){
  tree[,paste("DHP", i, sep="")] <- apply(tree[,c(paste("DHP", i-1, sep=""), i)], 1, sum, na.rm=T)
}

tree[,as.character(beg:2011)] <- tree[,as.character(beg:2011)]*-1
for (i in 2010:beg){
  tree[,paste("DHP", i, sep="")] <- apply(tree[,c(paste("DHP", i+1, sep=""), i+1)], 1, sum, na.rm=T)
}
tree[,as.character(beg:2011)] <- tree[,as.character(beg:2011)]*-1


blabla <- tree[!is.na(tree[,"2011"]),]  ### pour vérifier que les croissances sont bien répercutées sur les diamètres


####################################################
##                  save               ##
####################################################

write.csv(tree, "~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/input_chron/AbSUTtree.csv")
