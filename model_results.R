####################################################
####################################################

# Delete all objects in the work space
rm(list=ls(all=TRUE))

####################################################
##                  Data & Packages               ##
####################################################

setwd("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/output_model")
library(effects)


####################################################
##                parameters values               ##
####################################################
#file_names=as.list(dir(pattern="*.rdata"))
#file_names
#for(i in 1:length(file_names)) load(file_names[[i]])
#assign(namepar, par) pour assigner un nom

######################## BIC #######################
                                               ## Ab
load("bsimfixefBICAb.rdata")
par <- round(apply(bsimfixef,2,quantile, prob=c(0.025, 0.5, 0.975)), digits=3)
par
nbpar <- dim(bsimfixef)[2]
nbpar

signif <- matrix(nrow=nbpar, ncol=2)
signif <- as.data.frame(signif)
colnames(signif) <- c("Inf0", "Sup0")
for (i in 1:nbpar){
  rownames(signif)[i] <- colnames(bsimfixef)[i]
  signif[i,"Inf0"] <- sum(bsimfixef[,i]<0)/1000 ### specific hypothesis
  signif[i,"Sup0"] <- sum(bsimfixef[,i]>0)/1000 ### specific hypothesis

}

## une colonne signification: if <=0.1..... *, **, ***, ......




load("mod1BICAb.rdata")
summary(mod1)
adjr2 <- summary(mod1)$adj.r.squared
adjr2




####################################################
##            model for "effects" plots           ##
####################################################

######################## BIC #######################
                                               ## Ab
load("saveNCIBICAb.rdata")
load("saveclimBICAb.rdata")
load("modBICAb.rdata")
load("data_newBICAb.rdata")



### prendre seulement les significatifs à 10% ou 5%...
### ggplot facets?


## avec ggplot2
eff <- Effect(c("NCI","Tp9"),mod, se=T)
eff <- as.data.frame(eff)
eff$NCI<- as.factor(eff$NCI)
levels(eff$NCI) <- grid.pretty(range(saveNCI)) ## descale NCI
colnames(eff)[colnames(eff)=="fit"] <- "lBAI"
eff$Tp9 <- (eff$Tp9*sd(saveclim$Tp9)) + mean(saveclim$Tp9) ## descale clim






#  voir documentation package effects pour les plots

allEffects(mod)
plot(allEffects(mod), multiline=T)
plot(allEffects(mod), multiline=F)

eff <- Effect(c("P6","NCI"),mod, se=T)
summary(eff)
#plot(eff, multiline=F)
plot(eff, multiline=T)#, ci.style="bands")






## avec ggplot2
eff <- Effect(c("NCI","P6"),mod, se=T)
eff <- as.data.frame(eff)
eff$NCI<- as.factor(eff$NCI)
levels(eff$NCI) <- grid.pretty(range(saveNCI)) ## descale NCI
colnames(eff)[colnames(eff)=="fit"] <- "lBAI"
eff$P6 <- (eff$P6*sd(saveclim$P6)) + mean(saveclim$P6) ## descale clim


ggplot(eff, aes(x=Tp9,y=lBAI,group=NCI))+
geom_line(aes(linetype=NCI))+
theme_bw()



# ## colour=NA suppresses edges of the ribbon
# geom_ribbon(color=NA,alpha=0.1,aes(ymin=lower,ymax=upper))+
# ## add rug plot based on original data
#geom_rug(data=eff$data,aes(y=NULL),sides="b")
