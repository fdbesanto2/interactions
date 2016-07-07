####################################################
####################################################

# Delete all objects in the work space
rm(list=ls(all=TRUE))
options(digits=22)

####################################################
##                  Data & Packages               ##
####################################################

setwd("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/output_model")


load("modESSAI")


summary(mod)


###### predict

#  voir documentation package effects pour les plots

allEffects(mod)
plot(allEffects(mod), multiline=T)
plot(allEffects(mod), multiline=F)

eff <- Effect(c("Tp9","NCI"),mod, se=T)
summary(eff)
#plot(eff, multiline=F)
plot(eff, multiline=T)#, ci.style="bands")






## avec ggplot2
eff <- Effect(c("NCI","Tp9"),mod, se=T)
eff <- as.data.frame(eff)
eff$NCI<- as.factor(eff$NCI)
levels(eff$NCI) <- grid.pretty(range(saveNCI)) ## descale NCI
colnames(eff)[colnames(eff)=="fit"] <- "lBAI"
eff$Tp9 <- (eff$Tp9*sd(saveclim$Tp9)) + mean(saveclim$Tp9) ## descale clim


ggplot(eff, aes(x=Tp9,y=lBAI,group=NCI))+
geom_line(aes(linetype=NCI))+
theme_bw()



# ## colour=NA suppresses edges of the ribbon
# geom_ribbon(color=NA,alpha=0.1,aes(ymin=lower,ymax=upper))+
# ## add rug plot based on original data
#geom_rug(data=eff$data,aes(y=NULL),sides="b")
