# Delete all objects in the work space
rm(list=ls())


####################################################
##          set species and temporal window       ##
####################################################

# Load Package
library(dplR)
library(ggplot2)
# Import growth
rw <- read.tucson("~/ownCloud/Scan & measures/RWL/To1823epsme.rwl")
# import diameter & coordinates & IDs
dc <- read.csv("~/ownCloud/Work_directory/Data/Sampling/Duparquet/1823/1823_titre_point.csv", sep=";")
# columns format
dc$Y <- as.numeric(as.character(dc$Y))
dc$DHP11 <- as.numeric(as.character(dc$DHP11))*10

# set temporal window
beg <- 1991
end <- 2013
growth <- rw[rownames(rw)>=beg & rownames(rw)<=end ,]

####################################################
##                  growth data                   ##
####################################################
# changes individuals names in rw
colnames(growth) <- substr(colnames(growth), 4,7)

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
unique(dc$Sp)
dc[dc$Sp=="TOC ","Sp"] <- "TOC"
unique(dc$Sp)

# only living trees (no DBH measures for dead trees)
#tab <- dc[dc[,"MORT"] == "0",c("TAG", "Sp","DHP11", "X", "Y")
tab <- dc[dc[,"Statut11"] == "VC" | dc[,"Statut11"] == "VD" | dc[,"Statut11"] == "VP",c("TAG", "Sp","DHP11", "X", "Y")]
tab$TAG <- as.factor(tab$TAG)

# réduire la zone d'étude
tab <- tab[!is.na(tab$X),] # on supprime d'abord les arbres sans coordonnées
tab <- tab[!is.na(tab$Y),]

ggplot(data=tab, aes(X,Y))+
geom_point(aes(size=DHP11), alpha=0.5)



####################################################
##                   gather data                  ##
####################################################
## modifying tree names
growth[substr(growth$TAG,1,1) == 0,"TAG"] <- substr(growth[substr(growth$TAG,1,1) == 0,"TAG"],2,4)
growth[substr(growth$TAG,1,1) == 0,"TAG"] <- substr(growth[substr(growth$TAG,1,1) == 0,"TAG"],2,3)
growth[substr(growth$TAG,1,2) == "I0","TAG"] <- paste("I",substr(growth[substr(growth$TAG,1,2) == "I0","TAG"],3,4),sep="")
growth[substr(growth$TAG,1,2) == "I0","TAG"] <- paste("I",substr(growth[substr(growth$TAG,1,2) == "I0","TAG"],3,3),sep="")

growth[growth$TAG=="T42","TAG"] <- 3931
growth[growth$TAG=="T29","TAG"] <- 3433

tab$TAG <- as.character(tab$TAG)
tab[nchar(tab$TAG)>4,"TAG"] <- substr(tab[nchar(tab$TAG)>4,"TAG"],1,4)
tab$TAG <- as.factor(tab$TAG)

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
treei[treei$Sp=="ABA",]

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
  #
  #tree[tree$TAG==3229,]
  #tree[tree$TAG==3229,6:ncol(tree)] <- NA
  #tree[tree$TAG==3229,]
  #
# le ABA en trop dans les PGL: 3229


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


blabla <- tree[!is.na(tree[,"2000"]),]  ### pour vérifier que les croissances sont bien répercutées sur les diamètres



####################################################
##                  save               ##
####################################################

write.csv(tree, "~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/input_chron/To1823tree.csv")
