####################################################
####################################################
#                                                  #
#                 softwood / hardwood              #
#                                                  #
####################################################
####################################################

# Delete all objects in the work space
rm(list=ls(all=TRUE))
options(digits=22)

####################################################
##                  Data & Packages               ##
####################################################
# Load Package
library(arm)
library(ggplot2)
setwd("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/RUN_MODEL")
load("./clim.Rdata")
load("./chrono.Rdata")
load("./dist_inv.Rdata")
load("./dbh.Rdata")
hist(stockdist_inv[stockdist_inv>0])

#### sparse matrix (pour ne pas enregistrer les 0 = gagner temps + mémoire notamment lors des calculs matriciels)
stockdist_inv <- Matrix(stockdist_inv, sparse=T)

####################################################
##                      Model                     ##
####################################################
# set alpha & beta
alpha <-  8.3
beta  <- 1.7

# dbh or dbh^2 (dbh2) in the model?
diam <- "dbh"

# which sp/site
spsite <- "SUTBa" # ex: BICAr / D1823To


####################################################
# model formula
# BIC
if (spsite=="BICAs"){
  form <- lBAI ~ DC7+P6+NCIhard+NCIsoft+DBH+NCI:DC7+NCI:P6+(DC7+P6+NCIhard+NCIsoft+DBH+NCI:DC7+NCI:P6|TAG)
} else if (spsite=="BICAr"){
  form <- lBAI ~ DC7+P6+S3+NCIhard+NCIsoft+DBH+NCI:DC7+NCI:P6+NCI:S3+(DC7+P6+S3+NCIhard+NCIsoft+DBH+NCI:DC7+NCI:P6+NCI:S3|TAG)
} else if (spsite=="BICPt"){
  form <- lBAI ~ DC7+P6+S3+NCIhard+NCIsoft+DBH+NCI:DC7+NCI:P6+NCI:S3+(DC7+P6+S3+NCIhard+NCIsoft+DBH+NCI:DC7+NCI:P6+NCI:S3|TAG)
} else if (spsite=="BICAb"){
    form <- lBAI ~ GSLp+Tp9+Sp9+NCIhard+NCIsoft+DBH+NCI:GSLp+NCI:Tp9+NCI:Sp9+(GSLp+Tp9+Sp9+NCIhard+NCIsoft+DBH+NCI:GSLp+NCI:Tp9+NCI:Sp9|TAG)

# SUT
} else if (spsite=="SUTAb"){
    form <- lBAI ~ Sp11+S6+NCIhard+NCIsoft+DBH+NCI:Sp11+NCI:S6+(Sp11+S6+NCIhard+NCIsoft+DBH+NCI:Sp11+NCI:S6|TAG)
} else if (spsite=="SUTAs"){
    form <- lBAI ~ DC6+NCIhard+NCIsoft+DBH+NCI:DC6+(DC6+NCIhard+NCIsoft+DBH+NCI:DC6|TAG)
} else if (spsite=="SUTBa"){
    form <- lBAI ~ Pp10+S4+NCIhard+NCIsoft+DBH+NCI:Pp10+NCI:S4+(Pp10+S4+NCIhard+NCIsoft+DBH+NCI:Pp10+NCI:S4|TAG)

# ABI
} else if (spsite=="ABIAb"){
    form <- lBAI ~ GSLp+NCIhard+NCIsoft+DBH+NCI:GSLp+(GSLp+NCIhard+NCIsoft+DBH+NCI:GSLp|TAG)
} else if (spsite=="ABIAr"){
    form <- lBAI ~ Tp7+T3+NCIhard+NCIsoft+DBH+NCI:Tp7+NCI:T3+(Tp7+T3+NCIhard+NCIsoft+DBH+NCI:Tp7+NCI:T3|TAG)
} else if (spsite=="ABIAs"){
    form <- lBAI ~ T5+S6+DC7+P7+DC8+NCIhard+NCIsoft+DBH+NCI:T5+NCI:S6+NCI:DC7+NCI:P7+NCI:DC8+(T5+S6+DC7+P7+DC8+NCIhard+NCIsoft+DBH+NCI:T5+NCI:S6+NCI:DC7+NCI:P7+NCI:DC8|TAG)
} else if (spsite=="ABIPg"){
    form <- lBAI ~ DC5+T6+T8+NCIhard+NCIsoft+DBH+NCI:DC5+NCI:T6+NCI:T8+(DC5+T6+T8+NCIhard+NCIsoft+DBH+NCI:DC5+NCI:T6+NCI:T8|TAG)
} else if (spsite=="ABITo"){
    form <- lBAI ~ Pp10+T6+P8+NCIhard+NCIsoft+DBH+NCI:Pp10+NCI:T6+NCI:P8+(Pp10+T6+P8+NCIhard+NCIsoft+DBH+NCI:Pp10+NCI:T6+NCI:P8|TAG)

# D1823
} else if (spsite=="D1823Ab"){
    form <- lBAI ~ DC7+NCIhard+NCIsoft+DBH+NCI:DC7+(DC7+NCIhard+NCIsoft+DBH+NCI:DC7|TAG)
} else if (spsite=="D1823Pg"){
    form <- lBAI ~ DCp7+DCp9+P6+DC7+NCIhard+NCIsoft+DBH+NCI:DCp7+NCI:DCp9+NCI:P6+NCI:DC7+(DCp7+DCp9+P6+DC7+NCIhard+NCIsoft+DBH+NCI:DCp7+NCI:DCp9+NCI:P6+NCI:DC7|TAG)
} else if (spsite=="D1823Pt"){
    form <- lBAI ~ Tp8+S3+NCIhard+NCIsoft+DBH+NCI:Tp8+NCI:S3+(Tp8+S3+NCIhard+NCIsoft+DBH+NCI:Tp8+NCI:S3|TAG)
} else if (spsite=="D1823To"){
    form <- lBAI ~ T2+T6+T8+NCIhard+NCIsoft+DBH+NCI:T2+NCI:T6+NCI:T8+(T2+T6+T8+NCIhard+NCIsoft+DBH+NCI:T2+NCI:T6+NCI:T8|TAG)

# D1847
} else if (spsite=="D1847Ab"){
    form <- lBAI ~ Pp6+T4+NCIhard+NCIsoft+DBH+NCI:Pp6+NCI:T4+(Pp6+T4+NCIhard+NCIsoft+DBH+NCI:Pp6+NCI:T4|TAG)
} else if (spsite=="D1847Pg"){
    form <- lBAI ~ Pp7+T4+NCIhard+NCIsoft+DBH+NCI:Pp7+NCI:T4+(Pp7+T4+NCIhard+NCIsoft+DBH+NCI:Pp7+NCI:T4|TAG)
} else if (spsite=="D1847Pt"){
    form <- lBAI ~ DCp8+S3+NCIhard+NCIsoft+DBH+NCI:DCp8+NCI:S3+(DCp8+S3+NCIhard+NCIsoft+DBH+NCI:DCp8+NCI:S3|TAG)
} else if (spsite=="D1847To"){
    form <- lBAI ~ Tp8+Sp9+S1+T6+P8+NCIhard+NCIsoft+DBH+NCI:Tp8+NCI:Sp9+NCI:S1+NCI:T6+NCI:P8+(Tp8+Sp9+S1+T6+P8+NCIhard+NCIsoft+DBH+NCI:Tp8+NCI:Sp9+NCI:S1+NCI:T6+NCI:P8|TAG)
}


####################################################
##                 NCI Calculation                ##
####################################################
ptm <- proc.time() # start clock
dbh <- stockdbh^alpha ## on élève le dbh à la puissance Alpha
dist_inv <- stockdist_inv^beta ## on élève la distance a la puissance Beta
hist(dist_inv[dist_inv>0])
data[is.na(data$Sp),"Sp"] <- "IND"
sp <- unique(data$Sp)    ## nom des espèces sur la placette

for (s in sp) {
  dbh_tmp <- dbh
  dbh_tmp[data[rownames(dist_inv),"Sp"]!=s] <- 0
  NCI_dbh_dist <- as.vector(dist_inv %*% dbh_tmp)
  data <- cbind(data,NCI_dbh_dist)
}
proc.time() - ptm # stop clock

data$NNN <- rowSums(data[,substr(colnames(data),1,3)=="NCI"])## calcul du NCI produit par toutes les individus autour de l'arbre focal
colnames(data)[substr(colnames(data),1,3)=="NCI"] <- sp ## nom des colonnes de compétition pour chaque espèce
colnames(data)[length(colnames(data))] <- "NCI"


if (substr(spsite,1,3)=="BIC"){
  ## hardwood competition
  data$NCIhard <- rowSums(data[,c("BEPA","POTR","ACRU","ACSA","ACPE","SOAU","POBA","POGR","ACSP","SOAM", "QURU","SODE")])
  ## softwood competition
  data$NCIsoft <- rowSums(data[,c("ABBA","PIGL","THOC","PIRU")])
} else if (substr(spsite,1,3)=="SUT"){
  ## hardwood competition
  data$NCIhard <- rowSums(data[,c("BEPA","BEAL","FAGR","AMSP","ACRU","ACSA","ACPE","ACSP","SOAM","SODE","PRPE")])
  ## softwood competition
  data$NCIsoft <- rowSums(data[,c("ABBA","PIGL","PIRU","TSCA")])
} else if (substr(spsite,1,3)=="ABI"){
  ## hardwood competition
  data$NCIhard <- rowSums(data[,c("BEPA","ACRU","ACSP","ACSA","BEAL","PRPE")])
  ## softwood competition
  data$NCIsoft <- rowSums(data[,c("ABBA","THOC","PIST","PIGL")])
} else if (substr(spsite,1,5)=="D1823"){
  ## hardwood competition
  data$NCIhard <- rowSums(data[,c("BPA","PTR")])
  ## softwood competition
  data$NCIsoft <- rowSums(data[,c("ABA","TOC","PGL","PMA")])
} else if (substr(spsite,1,5)=="D1847"){
  ## hardwood competition
  data$NCIhard <- rowSums(data[,c("BPA","PTR")])
  ## softwood competition
  data$NCIsoft <- rowSums(data[,c("ABA","TOC","PGL")])
}

data_mod = data

####################################################
##                 Plot Competition               ##
####################################################

# data pour l'espèce a l'étude (semection des arbres qui ont des données de croissance pour 1991 ou 2000 ou 2010)
data_sp <- data[!is.na(data$X1991) | !is.na(data$X2000) | !is.na(data$X2010),]
# valeur median de competition pour échelle de couleurs
mid <- median(data_mod$NCI)

# plot
ggplot(data=data_mod, aes(X,Y))+
geom_point(aes(size=DHP11/10, color=NCI), show.legend=T)+
scale_color_gradient2(midpoint=mid,low="blue",mid="green",high="red",space="Lab")+
geom_point(aes(size=DHP11/10), alpha=0.25, show.legend=T, data=data_sp)+
theme(legend.position="right")

####################################################
####################################################
####################################################

beg <- min(clim$X)
end <- max(clim$X)

# empiler les données par année
data_new=cbind(data_mod[,c("TAG",paste("DHP", beg, sep=""),paste("X", beg, sep=""),sp,"NCI","NCIhard","NCIsoft")],1991)
colnames(data_new)[2]="DHP"
colnames(data_new)[3]="AGR"
colnames(data_new)[ncol(data_new)]="Year"

for(y in ((beg+1):end)){
  tmp=cbind(data_mod[,c("TAG",paste("DHP", y, sep=""),paste("X", y, sep=""),sp,"NCI","NCIhard","NCIsoft")],y)
  colnames(tmp)[2]="DHP"
  colnames(tmp)[3]="AGR"
  colnames(tmp)[ncol(data_new)]="Year"
  data_new=rbind(data_new,tmp)
}
data_new=data_new[!is.na(data_new$AGR),]     ## enlever les années pour lesquelles l'arbres n'a pas de largeur de cerne


####################################################
##          Remove Insect Outreak periods         ##
####################################################

if (spsite=="BICPt"){
  data_new <- data_new[data_new$Year< 2007 | data_new$Year> 2010 ,]
} else if (spsite=="BICAb"){
  data_new <- data_new[data_new$Year< 1977 | data_new$Year> 1985 ,]
} else if (spsite=="SUTAb"){
  data_new <- data_new[data_new$Year< 1972 | data_new$Year> 1983 ,]
} else if (spsite=="ABIPg"){
  data_new <- data_new[data_new$Year< 1981 | data_new$Year> 1991 ,]
} else if (spsite=="ABIAb"){
  data_new <- data_new[data_new$Year< 1981 | data_new$Year> 1991 ,]
} else if (spsite=="D1823Ab"){
  data_new <- data_new[data_new$Year< 1980 | data_new$Year> 1990 ,]
} else if (spsite=="D1823Pg"){
  data_new <- data_new[data_new$Year< 1974 | data_new$Year> 1987 ,]
} else if (spsite=="D1823Pt"){
  data_new <- data_new[data_new$Year< 2001 | data_new$Year> 2004 ,]
} else if (spsite=="D1847Ab"){
  data_new <- data_new[data_new$Year< 1974 | data_new$Year> 1990 ,]
} else if (spsite=="D1847Pg"){
  data_new <- data_new[data_new$Year< 1974 | data_new$Year> 1987 ,]
} else if (spsite=="D1847Pt"){
  data_new <- data_new[data_new$Year< 2001 | data_new$Year> 2005 ,]
}

####################################################
##               Basal Area Increment             ##
####################################################

## calcul de la surface terrière
data_new$BA <- pi*(data_new$DHP^2)/4
## croissance en diamètre (diamètre après croissance)
data_new$DHP2 <- data_new$DHP + data_new$AGR*2
## surface terrière après croissance
data_new$BA2 <- pi*(data_new$DHP2^2)/4
data_new$BAI <- data_new$BA2 - data_new$BA

####################################################
##                      Plot                      ##
####################################################

## individuals' growth (BAI)
ggplot(data=data_new, aes(Year,BAI, group=TAG))+
geom_line(aes(colour=TAG), show.legend=F)

## individuals' growth (radial growth)
ggplot(data=data_new, aes(Year,AGR, group=TAG))+
geom_line(aes(colour=TAG), show.legend=F)

## size time series (BA)
ggplot(data=data_new, aes(Year,BA, group=TAG))+
geom_line(aes(colour=TAG), show.legend=F)

## size time series (DHP)
ggplot(data=data_new, aes(Year,DHP, group=TAG))+
geom_line(aes(colour=TAG), show.legend=F)

## effect of DHP on BAI
ggplot(data=data_new, aes(DHP,BAI))+
geom_line(aes(colour=TAG), show.legend=F)

## effect of DHP on radial growth
ggplot(data=data_new, aes(DHP,AGR))+
geom_line(aes(colour=TAG), show.legend=F)

## effect of NCI on BAI
ggplot(data=data_new, aes(NCI,BAI))+
geom_line(aes(colour=TAG), show.legend=F)

####################################################
##       Transforming variables (log, ^2...)      ##
####################################################

data_new$DBH2 <- data_new$DHP^2  ## dbh au carré (pour les modèles)
data_new$lAGR <- log(data_new$AGR + 1)   ## log+1 de growth
data_new$lBAI <- log(data_new$BAI + 1)   ## log+1 de BAI
nb <- ncol(data_new)
## ajouter le climat
data_new=merge(data_new,clim,by.x="Year",by.y="X")

## effect of DHP on lBAI
ggplot(data=data_new, aes(DHP,lBAI))+
geom_line(aes(colour=TAG), show.legend=F)

## effect of DBH2 on lBAI
ggplot(data=data_new, aes(DBH2,lBAI))+
geom_line(aes(colour=TAG), show.legend=F)

## effect of NCI on lBAI
ggplot(data=data_new, aes(NCI,lBAI))+
geom_line(aes(colour=TAG), show.legend=F)

## effect of NCIhard on lBAI
ggplot(data=data_new, aes(NCIhard,lBAI))+
geom_line(aes(colour=TAG), show.legend=F)

## effect of NCIsoft on lBAI
ggplot(data=data_new, aes(NCIsoft,lBAI))+
geom_line(aes(colour=TAG), show.legend=F)

####################################################
##              Scaling variables                 ##
####################################################
## save initial values to descale the variables for effects plots
saveNCI <- data_new$NCI
saveclim <- data_new[,(nb+1):ncol(data_new)]

## scale
data_new$NCI <- scale(data_new$NCI)
data_new$NCIhard <- scale(data_new$NCIhard)
data_new$NCIsoft <- scale(data_new$NCIsoft)
data_new$DHP <- scale(data_new$DHP)
data_new$DBH2 <- scale(data_new$DBH2)

# clim variables
data_new[,(nb+1):ncol(data_new)] <- scale(data_new[,(nb+1):ncol(data_new)])
#
# data_new$DC5 <- scale(data_new$DC5)
# data_new$DC6 <- scale(data_new$DC6)
# data_new$DC7 <- scale(data_new$DC7)
# data_new$DC8 <- scale(data_new$DC8)
# data_new$DCp7 <- scale(data_new$DCp7)
# data_new$DCp8 <- scale(data_new$DCp8)
# data_new$DCp9 <- scale(data_new$DCp9)
# data_new$GSLp <- scale(data_new$GSLp)
# data_new$P6 <- scale(data_new$P6)
# data_new$P7 <- scale(data_new$P7)
# data_new$P8 <- scale(data_new$P8)
# data_new$Pp10 <- scale(data_new$Pp10)
# data_new$Pp6 <- scale(data_new$Pp6)
# data_new$Pp7 <- scale(data_new$Pp7)
# data_new$S1 <- scale(data_new$S1)
# data_new$S3 <- scale(data_new$S3)
# data_new$S4 <- scale(data_new$S4)
# data_new$S6 <- scale(data_new$S6)
# data_new$Sp11 <- scale(data_new$Sp11)
# data_new$Sp9 <- scale(data_new$Sp9)
# data_new$T2 <- scale(data_new$T2)
# data_new$T3 <- scale(data_new$T3)
# data_new$T4 <- scale(data_new$T4)
# data_new$T5 <- scale(data_new$T5)
# data_new$T6 <- scale(data_new$T6)
# data_new$T8 <- scale(data_new$T8)
# data_new$Tp7 <- scale(data_new$Tp7)
# data_new$Tp8 <- scale(data_new$Tp8)
# data_new$Tp9 <- scale(data_new$Tp9)

# sp
#data_new$ABBA <- scale(data_new$ABBA)
#data_new$POTR <- scale(data_new$POTR)
#data_new$ACSA <- scale(data_new$ACSA)
#data_new$ACRU <- scale(data_new$ACRU)
#data_new$BEPA <- scale(data_new$BEPA)

####################################################
##                      Model                     ##
####################################################
# set dbh or sbh^2
if (diam=="dbh"){
  data_new$DBH <- data_new$DHP
} else if (diam=="dbh2"){
  data_new$DBH <- data_new$DBH2
}

ptm <- proc.time() # start clock
mod <-  lmer(form, data=data_new, REML = FALSE)
proc.time() - ptm # stop clock
summary(mod)

# variance des effets mixtes
variance <- as.data.frame(VarCorr(mod))
variance <- variance[is.na(variance$var2),c("var1","vcov")]
variance$percvar <- variance$vcov*100/sum(variance$vcov)
variance <- cbind(variance[,1],round(variance[,2:3], digits=5))


# bayesian estimation of the parameters
nsim <- 1000
bsim <- sim(mod, n.sim=nsim)
#str(bsim)
par <- round(apply(bsim@fixef,2,quantile, prob=c(0.025, 0.5, 0.975)), digits=3)
par

####################################################
##         Plot posterior distribution            ##
####################################################

# representer les chaines (avec ou sans burnin)
nbpar <- dim(bsim@fixef)[2] ## nb of parameters

tab <- data.frame()
tab1 <- data.frame()
for (i in 1:nbpar){
  tab <- rbind(tab, as.data.frame(bsim@fixef[,i]))
  tab1 <- rbind(tab1, as.data.frame(rep(colnames(bsim@fixef)[i],nsim)))
}
posterior <- cbind(tab, tab1)
colnames(posterior) <- c("distrib", "parameter")

options(digits=2) ### seul moyen pour arrondir les étiquettes des axes
ggplot(data=posterior, aes(distrib, group=parameter))+
geom_histogram()+
geom_vline(xintercept = 0, color="red", size=0.25)+
facet_wrap(~parameter, scales="free")+
theme(axis.text.x=element_text(angle=30, hjust=1, size=7))
options(digits=7) ### retour à la valeur par défaut


####################################################
####################################################
####################################################

### specific hypothesis
sum(bsim@fixef[,7]>0)/nsim

# Plot predictions
## effect of NCI on BAI
data_pred <- data_new[,c("TAG", "Year", "lBAI")]
data_pred$meas <- seq(nrow(data_new))
predict <- predict(mod)
data_pred$predict <- predict


ggplot(data=data_pred, aes(meas, lBAI))+
geom_line()+
geom_line(aes(meas, predict), color="red", alpha=0.5)


ggplot(data=data_pred, aes(predict, lBAI))+
geom_point()+
geom_smooth(method=lm, color="red", size=0.5, se=T, fill="grey30")


mod1 <- lm(predict ~ lBAI, data=data_pred)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
par(mfrow=c(1,1))


####################################################
##            save parameters values              ##
####################################################
setwd("~/ownCloud/Work_directory/Analysis/chapitre_2/interactions/output_model")

## bayesian parameters estimates
namebsim <- paste("bsimfixef", spsite, ".rdata", sep="")
bsimfixef <- as.data.frame(bsim@fixef)
save(bsimfixef, file = namebsim)

## adjusted r2 (mod1)
par2 <- as.data.frame(matrix(ncol=3, nrow=1))
colnames(par2) <- c("adjr2", "alpha", "beta")
par2[,1] <- summary(mod1)$adj.r.squared
par2[,2] <- alpha
par2[,3] <- beta
namempar2 <- paste("par2", spsite, ".rdata", sep="")
save(par2, file = namempar2)

####################################################
##         save models for "effects" plots        ##
####################################################
# model
namemod <- paste("mod", spsite, ".rdata", sep="")
save(mod, file = namemod)

# data_new
namedata <- paste("data_new", spsite, ".rdata", sep="")
save(data_new, file = namedata)

# non scaled Clim & NCI
namesaveclim <- paste("saveclim", spsite, ".rdata", sep="")
save(saveclim, file = namesaveclim)

namesaveNCI <- paste("saveNCI", spsite, ".rdata", sep="")
save(saveNCI, file = namesaveNCI)
