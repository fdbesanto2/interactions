####################################################
####################################################

# Delete all objects in the work space
rm(list=ls(all=TRUE))

####################################################
##                  Data & Packages               ##
####################################################

setwd("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/output_model")
library(effects)
library(grid)
library(ggplot2)

####################################################
# open all ".rdata" files
#file_names=as.list(dir(pattern="*.rdata"))
#file_names
#for(i in 1:length(file_names)) load(file_names[[i]])
#assign(namepar, par) pour assigner un nom

####################################################
##               extraction function              ##
####################################################

#############################  "effects" plot  ###############################
funck <- function(signiflevel, plot, sp){
  nbpar <- dim(bsimfixef)[2]
  signif <- as.data.frame(matrix(nrow=nbpar, ncol=2))
  colnames(signif) <- c("Inf0", "Sup0")
  for (i in 1:nbpar){
    rownames(signif)[i] <- colnames(bsimfixef)[i]
    signif[i,"Inf0"] <- sum(bsimfixef[,i]<0)/1000 ### specific hypothesis
    signif[i,"Sup0"] <- sum(bsimfixef[,i]>0)/1000 ### specific hypothesis
  }
  signif$nchar <- nchar(rownames(signif))   # count the number of characters
  signif$lastchar <- substr(rownames(signif), signif$nchar-2, signif$nchar)
  signifNCI <- signif[signif$lastchar == "NCI",]
  signifNCI <- signifNCI[signifNCI$Inf0 <= signiflevel | signifNCI$Sup0 <= signiflevel,]
  if (dim(signifNCI)[1]>0){
    listparam <- substr(rownames(signifNCI), 1, signifNCI$nchar-4)

    effplot <- data.frame()
    for (i in listparam){
      eff <- as.data.frame(Effect(c("NCI",i),mod, se=T))
      eff$NCI <- as.factor(eff$NCI)
      #levels(eff$NCI) <- grid.pretty(range(saveNCI)) ## descale NCI
      colnames(eff)[colnames(eff)=="fit"] <- "lBAI"
      eff[,i] <- (eff[,i]*sd(saveclim[,i])) + mean(saveclim[,i]) ## descale clim
      colnames(eff)[2] <- "clim"
      eff$plotsppar <- paste(plot, sp, i , sep="")
      effplot <- rbind(effplot, eff)

    }
  } else {effplot <- "ns"}

return (effplot)
}   #fin de la fonction


######################## BIC #######################
###################################################### Ab
p <- "BIC"
s <- "Ab"
load("saveNCIBICAb.rdata")
load("saveclimBICAb.rdata")
load("modBICAb.rdata")
load("data_newBICAb.rdata")
load("bsimfixefBICAb.rdata")
load("par2BICAb.rdata")
par <- round(apply(bsimfixef,2,quantile, prob=c(0.025, 0.5, 0.975)), digits=3)
par2
eff <- funck(signiflevel = 0.05, plot = p, sp = s)
assign(paste("eff",p,s,sep=""), eff)
alleff <- eff
###################################################### Ar
p <- "BIC"
s <- "Ar"
load("saveNCIBICAr.rdata")
load("saveclimBICAr.rdata")
load("modBICAr.rdata")
load("data_newBICAr.rdata")
load("bsimfixefBICAr.rdata")
load("par2BICAr.rdata")
par <- round(apply(bsimfixef,2,quantile, prob=c(0.025, 0.5, 0.975)), digits=3)
par2
eff <- funck(signiflevel = 0.05, plot = p, sp = s)
assign(paste("eff",p,s,sep=""), eff)
alleff <- rbind(alleff, eff)
######################################################



####################################################
##                     tables                     ##
####################################################



####################################################
##                     figures                    ##
####################################################

ggplot(alleff, aes(x=clim,y=lBAI,group=NCI))+
geom_line(aes(linetype=NCI))+
facet_grid(. ~ plotsppar, scales = "free", space = "free") +  # free / fixed / free_x / ..._y
theme_bw()









#
#
# ##############################################################
# ##############################################################
# ##############################################################
#
#
# ### prendre seulement les significatifs à 10% ou 5%...
# ### ggplot facets?
#
#
#
# #  voir documentation package effects pour les plots
#
# allEffects(mod)
# plot(allEffects(mod), multiline=T)
# plot(allEffects(mod), multiline=F)
#
# eff <- Effect(c("P6","NCI"),mod, se=T)
# summary(eff)
# #plot(eff, multiline=F)
# plot(eff, multiline=T)#, ci.style="bands")
#
#
#
#
#
#
# ## avec ggplot2
# eff <- Effect(c("NCI","P6"),mod, se=T)
# eff <- as.data.frame(eff)
# eff$NCI<- as.factor(eff$NCI)
# levels(eff$NCI) <- grid.pretty(range(saveNCI)) ## descale NCI
# colnames(eff)[colnames(eff)=="fit"] <- "lBAI"
# eff$P6 <- (eff$P6*sd(saveclim$P6)) + mean(saveclim$P6) ## descale clim
#
#
# ggplot(effplot, aes(x=clim,y=lBAI,group=NCI))+
# geom_line(aes(linetype=NCI))+
# theme_bw()
#
#
#
# # ## colour=NA suppresses edges of the ribbon
# # geom_ribbon(color=NA,alpha=0.1,aes(ymin=lower,ymax=upper))+
# # ## add rug plot based on original data
# #geom_rug(data=eff$data,aes(y=NULL),sides="b")
