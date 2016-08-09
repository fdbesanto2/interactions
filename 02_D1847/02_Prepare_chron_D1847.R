# Delete all objects in the work space
rm(list=ls())


####################################################
##          set species and temporal window       ##
####################################################

# Load Package
library(dplR)
library(ggplot2)
# Import growth
rw <- read.tucson("~/ownCloud/Scan & measures/RWL/Pg1847epsme.rwl")
# import diameter & coordinates & IDs
dc <- read.csv("~/ownCloud/Work_directory/Data/Sampling/Duparquet/1847/1847_diez_point.csv", sep=";")
# columns format
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

# only living trees (no DBH measures for dead trees)
#tab <- dc[dc[,"MORT"] == "0",c("TAG", "Sp","DHP11", "X", "Y")
tab <- dc[dc[,"Statut11"] == "VC" | dc[,"Statut11"] == "VD" | dc[,"Statut11"] == "VP" | dc[,"Statut11"] == "VCASSÉ" ,c("TAG", "Sp","DHP11", "X", "Y")]
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


# pour les PGL
growth[growth$TAG=="226","TAG"] <- "V226"
growth[growth$TAG=="496","TAG"] <- "V496"
growth[growth$TAG=="371","TAG"] <- "V371"
growth[growth$TAG=="373","TAG"] <- "V373"
growth[growth$TAG=="95","TAG"] <- "V95"
growth[growth$TAG=="390","TAG"] <- "V390"
growth[growth$TAG=="100","TAG"] <- "V100"
growth[growth$TAG=="396","TAG"] <- "V396"
growth[growth$TAG=="147","TAG"] <- "V147"
growth[growth$TAG=="409","TAG"] <- "V409"
growth[growth$TAG=="437","TAG"] <- "V437"
growth[growth$TAG=="467","TAG"] <- "V467"
growth[growth$TAG=="217","TAG"] <- "V217"
growth[growth$TAG=="495","TAG"] <- "V495"
growth[growth$TAG=="711","TAG"] <- "V711"
growth[growth$TAG=="393","TAG"] <- "V393"
growth[growth$TAG=="864","TAG"] <- "V864"
growth[growth$TAG=="387","TAG"] <- "V387"
growth[growth$TAG=="773","TAG"] <- "V773"
growth[growth$TAG=="860","TAG"] <- "V860"
growth[growth$TAG=="792","TAG"] <- "V792"
growth[growth$TAG=="400","TAG"] <- "V400"
growth[growth$TAG=="765","TAG"] <- "V765"
growth[growth$TAG=="489","TAG"] <- "V489"
growth[growth$TAG=="173","TAG"] <- "V173"
growth[growth$TAG=="879","TAG"] <- "K879"
growth[growth$TAG=="817","TAG"] <- "K817"
growth[growth$TAG=="826","TAG"] <- "K826"
growth[growth$TAG=="199","TAG"] <- "V199"
growth[growth$TAG=="279","TAG"] <- "V279"
growth[growth$TAG=="735","TAG"] <- "V735"
growth[growth$TAG=="804","TAG"] <- "V804"
growth[growth$TAG=="853","TAG"] <- "K853"
growth[growth$TAG=="857","TAG"] <- "V857"
growth[growth$TAG=="863","TAG"] <- "V863"
growth[growth$TAG=="912","TAG"] <- "V912"
growth[growth$TAG== "945","TAG"] <- "K945"
growth[growth$TAG== "720","TAG"] <- "V720"
growth[growth$TAG== "793","TAG"] <- "V793"


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


# # pour les PGL
tree[tree$TAG=="183" | tree$TAG=="907" | tree$TAG=="K879" | tree$TAG=="V860",6:ncol(tree)] <- NA

# pour les PGL
# tree[tree$TAG=="K826",]
# tree[tree$TAG=="K826",6:ncol(tree)] <- NA
# tree[tree$TAG=="K826",]


# # pour les TOC
# tree[tree$TAG=="809",]
# tree[tree$TAG=="809",6:ncol(tree)] <- NA
# tree[tree$TAG=="809",]
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


blabla <- tree[!is.na(tree[,"1992"]),]  ### pour vérifier que les croissances sont bien répercutées sur les diamètres



####################################################
##                  save               ##
####################################################

write.csv(tree, "~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/input_chron/Pg1847tree.csv")
