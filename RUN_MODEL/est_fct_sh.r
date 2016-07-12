####################################################
##                     function                   ##
####################################################

options(digits=22)


fn <- function(alpha=1,beta=1){

  #### changer de formule, type --> fn = function(par, data)

  ####################################################
  ##                      Model                     ##
  ####################################################

  # dbh or dbh^2 (dbh2) in the model?
  diam <- "dbh"  # dbh or dbh2

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
  }

  ####################################################
  ##                 NCI Calculation                ##
  ####################################################
  dbh <- stockdbh^alpha ## on élève le dbh à la puissance Alpha
  dist_inv <- stockdist_inv^beta ## on élève la distance a la puissance Beta

  sp <- unique(meas$Sp)    ## nom des espèces sur la placette

  for (s in sp) {
    dbh_tmp <- dbh
    dbh_tmp[meas[rownames(dist_inv),"Sp"]!=s] <- 0
    NCI_dbh_dist <- as.vector(dist_inv %*% dbh_tmp)
    meas <- cbind(meas,NCI_dbh_dist)
  }

  meas$NNN <- rowSums(meas[,substr(colnames(meas),1,3)=="NCI"]) ## calcul du NCI produit par toutes les individus autour de l'arbre focal
  colnames(meas)[substr(colnames(meas),1,3)=="NCI"] <- sp ## nom des colonnes de compétition pour chaque espèce
  colnames(meas)[length(colnames(meas))] <- "NCI"


  if (substr(spsite,1,3)=="BIC"){
    ## hardwood competition
    meas$NCIhard <- rowSums(meas[,c("BEPA","POTR","ACRU","ACSA","ACPE","SOAU","POBA","POGR","ACSP","SOAM", "QURU","SODE")])

    ## softwood competition
    meas$NCIsoft <- rowSums(meas[,c("ABBA","PIGL","THOC","PIRU")])
  } else if (substr(spsite,1,3)=="SUT"){
    ## hardwood competition
    meas$NCIhard <- rowSums(meas[,c("BEPA","BEAL","FAGR","AMSP","ACRU","ACSA","ACPE","ACSP","SOAM","SODE","PRPE")])

    ## softwood competition
    meas$NCIsoft <- rowSums(meas[,c("ABBA","PIGL","PIRU","TSCA")])
  }
  data_mod = meas

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
  }else if (spsite=="SUTAb"){
    data_new <- data_new[data_new$Year< 1972 | data_new$Year> 1983 ,]
  }

  ####################################################
  ##               Basal Area Increment             ##
  ####################################################

  ## calcul de la surface terrière
  data_new$BA <- pi*(data_new$DHP^2)/4
  ## croissance en diamètre (diamètre après croissance)
  data_new$DHP2 <- data_new$DHP + data_new$AGR*2
  ## surface terrière après diamètre
  data_new$BA2 <- pi*(data_new$DHP2^2)/4
  data_new$BAI <- data_new$BA2 - data_new$BA

  ####################################################
  ##       Transforming variables (log, ^2...)      ##
  ####################################################

  data_new$DBH2 <- data_new$DHP^2  ## dbh au carré (pour les modèles)
  data_new$lAGR <- log(data_new$AGR + 1)   ## log+1 de growth
  data_new$lBAI <- log(data_new$BAI + 1)   ## log+1 de BAI
  nb <- ncol(data_new)
  ## ajouter le climat
  data_new=merge(data_new,clim,by.x="Year",by.y="X")

  ####################################################
  ##              Scaling variables                 ##
  ####################################################

  ## scale
  data_new$NCI <- scale(data_new$NCI)
  data_new$NCIhard <- scale(data_new$NCIhard)
  data_new$NCIsoft <- scale(data_new$NCIsoft)
  data_new$DHP <- scale(data_new$DHP)
  data_new$DBH2 <- scale(data_new$DBH2)

  # clim variables
  data_new[,(nb+1):ncol(data_new)] <- scale(data_new[,(nb+1):ncol(data_new)])

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


  suppressWarnings(mod <-  lmer(form, data=data_new, REML = FALSE))
  return (c(alpha, beta, logLik(mod))) # *-1 = inverse de la log-likelihood car GenSA recherche les minimums
}
