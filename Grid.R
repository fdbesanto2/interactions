# Delete all objects in the work space
rm(list=ls(all=TRUE))
options(digits=22)

####################################################
##                  Data & Packages               ##
####################################################
library(ggplot2)
setwd("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/RUN_MODEL")

####################################################
##             alpha & beta estimation            ##
####################################################

## limite de l'espace alpha/beta
p <- 5000 # nb de point voulu sur la grille
alpha1 <- seq(0, 10, length.out = round(sqrt(p)))
beta1 <- seq(0, 10, length.out = round(sqrt(p)))
d1 <- expand.grid(alpha1 = alpha1, beta1 = beta1)
d1[d1==0]<- 0.0000001  # pas de 0 car sinon ca plante...

write.table(d1, file="grid.txt", row.names=F, col.names=F, quote=F, sep=",")
grid <- read.table("grid.txt", sep=",")

ggplot(data=grid, aes(V1,V2))+
geom_point()+
scale_x_continuous(breaks=seq(1,10,1))+
scale_y_continuous(breaks=seq(1,10,1))

# précision de la grille (nombre de digits après la virgule)
grid[grid[,1]==2 & grid[,2]>=1 & grid[,2]<=2, ]
