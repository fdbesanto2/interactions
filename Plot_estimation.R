####################################################
##             Plot estimation results            ##
####################################################
# Delete all objects in the work space
rm(list=ls(all=TRUE))
options(digits=22)

setwd("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/RUN MODEL")
load("./result_estimate.Rdata")

#load("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/04 - BIC/output/PtBICsh.rdata")
#load("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/RUN MODEL/ArBICsh10.rdata")

####################################################
##                       raster                   ##
####################################################

library(ggplot2)
out <- as.data.frame(out)
colnames(out) <- c("alpha", "beta", "ll")
point <- paste("alpha =",round(out[out[,3]==max(out[,3]),"alpha"],3), " beta =",  round(out[out[,3]==max(out[,3]),"beta"],3))

ggplot(data=out, aes(alpha, beta))+
geom_raster(aes(fill=ll), interpolate=T)+
scale_fill_gradientn(colours=c("blue","green","red"))+
geom_point(data= out[out[,3]==max(out[,3]),], color="black")+
ggtitle(point)

#ggsave("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/04 - BIC/output/PtBICsh.jpg", width=5, height= 5, dpi=600)


####################################################
##                       3d                       ##
####################################################

stock <- out

stock[stock[,3]==max(stock[,3]),]  ## meilleur valeur
qsdf <- stock[order(-stock[,3]),]
head(qsdf)

## plot en 3D
library(akima)
library(rgl)
s <- interp(stock[,1],stock[,2],stock[,3])
dim(s$z)
## gradient de couleur
zlim <- range(s$z)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- topo.colors(zlen,alpha=0) # height color lookup table
col <- colorlut[ s$z-zlim[1]+1 ] # assign colors to heights for each point

a <- 1 # pour applatir les graphiques trop hauts
surface3d(s$x,s$y,s$z/a, color=col, alpha=1, back="fill")
axes3d(c('x-+','y-+','z--'))
title3d(xlab="alpha", ylab="beta", zlab="loglik")
box3d()
grid3d(c("z++", "z--", "x", "y"), at = NULL, col = "gray", lwd = 1, lty = 1, n = 10)
play3d(spin3d(axis = c(0, 0, 1),rpm=10), duration=Inf)
