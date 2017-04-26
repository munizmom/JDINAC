library(mvtnorm) # data generate
library(Matrix) # data generate
library(KernSmooth)
library(akima)
library(glmnet)
library("randomForest")


PreS4 <- function(size0=300,size1=300,nP=100,Pdim=100,lim=2,s=2,e=1,times=100,thres=0,nsplit=10,nfolds=10,seed=100){

  set.seed(seed)
  seed2 <- sample(100:2016000,times)
  class0 <- array(runif(size0*nP,-lim,lim),dim=c(size0,nP))
  class0[,(Pdim/2+1):Pdim] <- exp(class0[,(Pdim/2+1):Pdim]) + rmvnorm(n=size0,sigma=diag(e^2,Pdim/2)) - s
  class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
  cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
  class1[,1:(Pdim/2)] <- cla
  class1[,(Pdim/2+1):Pdim] <- exp(cla) + rmvnorm(n=size1,sigma=diag(e^2,Pdim/2)) - s
  ylabel <- c(rep(0,size0),rep(1,size1))
  # edgeid <- list()
  cPre <- list()
  PLR <- list()
  for (i in 1:times){
    set.seed(seed2[i])   
    netData <- rbind(class0,class1)
    deltaR <- abs(cor(class0) - cor(class1))
    deltaR[lower.tri(deltaR)] <- 0
    Tid <- array(1:Pdim,dim=c(Pdim/2,2))
    deltaR[Tid] <- 1
    EDGE <- which(deltaR> thres, arr.ind=T)
    # nrow(EDGE)
    # edgeid[[i]] <- nrow(deltaR)*(EDGE[,2] - 1) + EDGE[,1]
    # sum(EDGEid %in% edgeid)

    class0 <- array(runif(size0*nP,-lim,lim),dim=c(size0,nP))
    class0[,(Pdim/2+1):Pdim] <- exp(class0[,(Pdim/2+1):Pdim]) + rmvnorm(n=size0,sigma=diag(e^2,Pdim/2)) - s
    class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
    cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
    class1[,1:(Pdim/2)] <- cla
    class1[,(Pdim/2+1):Pdim] <- exp(cla) + rmvnorm(n=size1,sigma=diag(e^2,Pdim/2)) - s
    DataPre <- rbind(class0,class1)
    
    cPre[[i]] <- classPre(EDGE=EDGE,classLabel=ylabel,DataFit=netData,DataPre=DataPre,nsplit=nsplit,nfolds=nfolds)
    
    # penPLR
    crossX <- (netData[,EDGE[,1]]-colMeans(netData[,EDGE[,1]])) * (netData[,EDGE[,2]]-colMeans(netData[,EDGE[,2]]))
    crossX2 <- (DataPre[,EDGE[,1]]-colMeans(DataPre[,EDGE[,1]])) * (DataPre[,EDGE[,2]]-colMeans(DataPre[,EDGE[,2]]))
    cv.fit3 <- cv.glmnet(x=crossX, y=ylabel, family = "binomial", nfolds=nfolds) 
    yp3 <- predict(cv.fit3,newx=crossX2, s="lambda.min",type="response") 
    coefs <- which(coef(cv.fit3,s="lambda.min")[-1] !=0)
    PLR[[i]] <- list(yPre=c(yp3),v=EDGE[coefs,])
}
 list(cPre=cPre,PLR=PLR)
}


PreS4RF <- function(size0=300,size1=300,nP=100,Pdim=100,lim=2,s=2,e=1,times=100,seed=100){

  set.seed(seed)
  seed2 <- sample(100:2016000,times)
  class0 <- array(runif(size0*nP,-lim,lim),dim=c(size0,nP))
  class0[,(Pdim/2+1):Pdim] <- exp(class0[,(Pdim/2+1):Pdim]) + rmvnorm(n = size0,sigma=diag(e^2,Pdim/2)) - s
  class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
  cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
  class1[,1:(Pdim/2)] <- cla
  class1[,(Pdim/2+1):Pdim] <- exp(cla) + rmvnorm(n = size1,sigma=diag(e^2,Pdim/2)) - s
  ylabel <- c(rep(0,size0),rep(1,size1))
  cPre <- list()
  PLR <- list()
  for (i in 1:times){
    set.seed(seed2[i])   
    netData <- rbind(class0,class1)

    class0 <- array(runif(size0*nP,-lim,lim),dim=c(size0,nP))
    class0[,(Pdim/2+1):Pdim] <- exp(class0[,(Pdim/2+1):Pdim]) + rmvnorm(n = size0,sigma=diag(e^2,Pdim/2)) - s
    class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
    cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
    class1[,1:(Pdim/2)] <- cla
    class1[,(Pdim/2+1):Pdim] <- exp(cla) + rmvnorm(n = size1,sigma=diag(e^2,Pdim/2)) - s
    DataPre <- rbind(class0,class1)

    RF <- randomForest(x=netData, y=as.factor(ylabel))
    votes <- predict(RF,DataPre,type="vote")
    cPre[[i]] <- list(yPre=votes[,2])    
   
    # penPLR2
    cv.fit3 <- cv.glmnet(x=netData, y=ylabel, family = "binomial", nfolds=5) 
    yp3 <- predict(cv.fit3,newx=DataPre, s="lambda.min",type="response") 
    PLR[[i]] <- list(yPre=c(yp3))
}
 list(cPre=cPre,PLR=PLR)
}


RFPreSimu4 <- function(size0=300,size1=300,nP=100,Pdim=100,lim=2,s=2,e=1,times=100,thres=0,seed=100){

  set.seed(seed)
  seed2 <- sample(100:2016000,times)
  class0 <- array(runif(size0*nP,-lim,lim),dim=c(size0,nP))
  class0[,(Pdim/2+1):Pdim] <- exp(class0[,(Pdim/2+1):Pdim]) + rmvnorm(n = size0,sigma=diag(e^2,Pdim/2)) - s
  class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
  cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
  class1[,1:(Pdim/2)] <- cla
  class1[,(Pdim/2+1):Pdim] <- exp(cla) + rmvnorm(n = size1,sigma=diag(e^2,Pdim/2)) - s
  ylabel <- c(rep(0,size0),rep(1,size1))
  # edgeid <- list()
  cPre <- list()
  for (i in 1:times){
    set.seed(seed2[i])   
    netData <- rbind(class0,class1)
    deltaR <- abs(cor(class0) - cor(class1))
    deltaR[lower.tri(deltaR)] <- 0
    Tid <- array(1:Pdim,dim=c(Pdim/2,2))
    deltaR[Tid] <- 1
    EDGE <- which(deltaR> thres, arr.ind=T)
    # nrow(EDGE)
    # edgeid[[i]] <- nrow(deltaR)*(EDGE[,2] - 1) + EDGE[,1]
    # sum(EDGEid %in% edgeid)

    class0 <- array(runif(size0*nP,-lim,lim),dim=c(size0,nP))
    class0[,(Pdim/2+1):Pdim] <- exp(class0[,(Pdim/2+1):Pdim]) + rmvnorm(n = size0,sigma=diag(e^2,Pdim/2)) - s
    class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
    cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
    class1[,1:(Pdim/2)] <- cla
    class1[,(Pdim/2+1):Pdim] <- exp(cla) + rmvnorm(n = size1,sigma=diag(e^2,Pdim/2)) - s
    DataPre <- rbind(class0,class1)
    
    cPre[[i]] <- classPreRF(EDGE=EDGE,classLabel=ylabel,DataFit=netData,DataPre=DataPre)
    
}
 list(cPre=cPre)
}

NBPreSimu4 <- function(size0=300,size1=300,nP=100,Pdim=100,lim=2,s=2,e=1,times=100,thres=0,seed=100){

  set.seed(seed)
  seed2 <- sample(100:2016000,times)
  class0 <- array(runif(size0*nP,-lim,lim),dim=c(size0,nP))
  class0[,(Pdim/2+1):Pdim] <- exp(class0[,(Pdim/2+1):Pdim]) + rmvnorm(n = size0,sigma=diag(e^2,Pdim/2)) - s
  class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
  cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
  class1[,1:(Pdim/2)] <- cla
  class1[,(Pdim/2+1):Pdim] <- exp(cla) + rmvnorm(n = size1,sigma=diag(e^2,Pdim/2)) - s
  ylabel <- c(rep(0,size0),rep(1,size1))
  # edgeid <- list()
  cPre <- list()
  for (i in 1:times){
    set.seed(seed2[i])   
    netData <- rbind(class0,class1)
    deltaR <- abs(cor(class0) - cor(class1))
    deltaR[lower.tri(deltaR)] <- 0
    Tid <- array(1:Pdim,dim=c(Pdim/2,2))
    deltaR[Tid] <- 1
    EDGE <- which(deltaR> thres, arr.ind=T)
    # nrow(EDGE)
    # edgeid[[i]] <- nrow(deltaR)*(EDGE[,2] - 1) + EDGE[,1]
    # sum(EDGEid %in% edgeid)

    class0 <- array(runif(size0*nP,-lim,lim),dim=c(size0,nP))
    class0[,(Pdim/2+1):Pdim] <- exp(class0[,(Pdim/2+1):Pdim]) + rmvnorm(n = size0,sigma=diag(e^2,Pdim/2)) - s
    class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
    cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
    class1[,1:(Pdim/2)] <- cla
    class1[,(Pdim/2+1):Pdim] <- exp(cla) + rmvnorm(n = size1,sigma=diag(e^2,Pdim/2)) - s
    DataPre <- rbind(class0,class1)
    
    cPre[[i]] <- classPreNB(EDGE=EDGE,classLabel=ylabel,DataFit=netData,DataPre=DataPre)
    
}
 list(cPre=cPre)
}
