# Delete all objects in the work space
rm(list=ls(all=TRUE))
options(digits=22)

####################################################
##                  Packages               ##
####################################################
# Load Package
library('doParallel')
library('Matrix')
source('./est_fct_sh.r')


load("./clim.Rdata")
load("./chrono.Rdata")
load("./dist_inv.Rdata")
load("./dbh.Rdata")

#### sparse matrix (pour ne pas enregistrer les 0 = gagner temps + mémoire notamment lors des calculs matriciels)
# stockdist_inv <- Matrix(stockdist_inv, sparse=T)

meas <- data

# Get node addresses given by Torque and assigned to your task
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile,header=FALSE) # Node adresses

# Get the number of cores by node
nCores <- detectCores()

# Open the clusters
cl <- makePSOCKcluster(rep(as.character(hostlist$V1),nCores))

# We replicated the addresses nCores times.
registerDoParallel(cl)

# load parameters
input <- read.table('grid.txt',sep=',')
names(input) <- c('alpha','beta')

out <- foreach(i=1:nrow(input),.packages='lme4',.export=c('meas'),.combine='rbind') %dopar% {
  fn(alpha = input$alpha[i], beta = input$beta[i])
}

save(out,file='./result_estimate.rdata')
