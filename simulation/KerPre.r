library(mvtnorm) # data generate
library(Matrix) # data generate
library(KernSmooth)
library(akima)
library(glmnet)

denPre2D <- function(edge,classLabel,DataFit,newData,method=c("integers","bandwidth")[2]){
  
  Index0 <- which(classLabel==0)
  Index1 <- which(classLabel==1)
  x <- DataFit[c(Index0,Index1),edge[1]]
  y <- DataFit[c(Index0,Index1),edge[2]]
  x0 <- newData[,edge[1]]
  y0 <- newData[,edge[2]]
  N <- max(c(x,y))
  bgI <- 1:length(Index0)
  fgI <- length(Index0)+(1:length(Index1))
  
  fgWidth <- c(bw.nrd0(x[fgI]),bw.nrd0(y[fgI]))
  bgWidth <- c(bw.nrd0(x[bgI]),bw.nrd0(y[bgI]))
  gridSize <- switch(method,
                     integers  = c(N, N),
                     bandwidth = ceiling(N / c(min(fgWidth[1], bgWidth[1]), 
                                               min(fgWidth[2], bgWidth[2]))))
  gridSize <- pmax(gridSize,10) # make sure there are at least 100 points in total

  fgSmooth <- bkde2D(x=cbind(x[fgI], y[fgI]), bandwidth=fgWidth,  gridsize=gridSize)
  fgP <- fgSmooth$fhat
  bgSmooth <- bkde2D(x=cbind(x[bgI], y[bgI]), bandwidth=bgWidth, gridsize=gridSize)
  bgP <- bgSmooth$fhat

  # make sure there are no zeros in the smooth function (since we will take a log of that)
  # fgP[fgP==0] <- min(fgP[fgP>0])/100
  fgP <- pmax(fgP, 1e-10)
  # bgP[bgP==0] <- min(bgP[bgP>0])/100
  bgP <- pmax(bgP, 1e-10)

  fgfit <- bicubic(x=fgSmooth$x1, y=fgSmooth$x2, z=fgP, x0=x0,y0=y0)
  bgfit <- bicubic(x=bgSmooth$x1, y=bgSmooth$x2, z=bgP, x0=x0,y0=y0)
  fgfit$z <- pmax(fgfit$z, 1e-10)
  bgfit$z <- pmax(bgfit$z, 1e-10)
  denPre <- log(fgfit$z/bgfit$z)
  denPre
}

classPre <- function(EDGE,classLabel,DataFit,DataPre,nsplit=10,nfolds=10){
  
  if(missing(DataPre)) 
    DataPre <- DataFit
  
  if(nsplit==0){       
    denX <- apply(EDGE,1,denPre2D,classLabel=classLabel,DataFit=DataFit,
                   newData=rbind(DataFit,DataPre))
    y <- classLabel
    cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
    yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
    yPre <- c(yp)
    Vid <- which(coef(cv.fit,s="lambda.min")[-1] !=0)
    Eset <- EDGE[Vid,]
    colnames(Eset) <- c("row","col")    
  } else {
    preY <- NULL
    vset <- NULL
    for(i in 1:nsplit){
      size0 <- sum(classLabel==0)
      size1 <- sum(classLabel==1)
      sn0 <- round(size0/2)
      sn1 <- round(size1/2)
      splitid <- c(sample(1:size0,sn0),sample((1:size1)+size0,sn1))
    
      cLabel <- classLabel[splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[splitid,],
                    newData=rbind(DataFit[-splitid,],DataPre))
      y <- classLabel[-splitid]
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-1] !=0)
      vset <- c(vset,coefs)
    
      cLabel <- classLabel[-splitid]        
      denX <- apply(EDGE,1,denPre2D,classLabel=cLabel,DataFit=DataFit[-splitid,],
                    newData=rbind(DataFit[splitid,],DataPre))
      y <- classLabel[splitid]
      cv.fit <- cv.glmnet(x=denX[1:length(y),], y=y, family = "binomial", nfolds=nfolds) 
      yp <- predict(cv.fit,newx=denX[-(1:length(y)),], s="lambda.min",type="response")
      preY <- cbind(preY,yp) 
      coefs <- which(coef(cv.fit,s="lambda.min")[-1] !=0)
      vset <- c(vset,coefs)     
    } 
    yPre <- rowMeans(preY) 
    numb <- table(vset)
    Vid <- as.numeric(rownames(numb))  
    Eset <- cbind(EDGE[Vid,],numb)
    Eset <- Eset[order(Eset[,3],decreasing=T),]
    colnames(Eset) <- c("row","col","numb")
  } 
  list(yPre=yPre,Eset=Eset) 
}

PreSimu <- function(size0=300,size1=300,Pdim=100,Pblock=10,nblock=2,rho=.5,rD=.5,rC=.5,
                    thres=0,times=100,nsplit=10,nfolds=10,seed=100){

  Sigma00 <- rho^abs(outer(1:(Pdim-Pblock*nblock),1:(Pdim-Pblock*nblock),"-"))
  Sigma01 <- rC* (-1)^abs(outer(1:(Pblock),1:(Pblock),"-"))
  diag(Sigma01) <- 1
  SS0 <- rep(list(Sigma01),nblock)
  SS0 <- bdiag(SS0)
  Sigma0 <- as.matrix(bdiag(Sigma00, SS0))

  Sigma10 <- rho^abs(outer(1:(Pdim-Pblock*nblock),1:(Pdim-Pblock*nblock),"-"))
  Sigma11 <- array(rD,dim=c(Pblock,Pblock))
  diag(Sigma11) <- 1
  SS1 <- rep(list(Sigma11),nblock)
  SS1 <- bdiag(SS1)
  Sigma1 <- as.matrix(bdiag(Sigma10, SS1))
  ylabel <- c(rep(0,size0),rep(1,size1))
  difR <- abs(Sigma1-Sigma0)
  difR[lower.tri(difR)] <- 0
  EDGEset <- which(difR> 0,arr.ind=T)
  EDGEid <- which(difR> 0)
  set.seed(seed)
  seed2 <- sample(100:2016000,times)
  class0 <- rmvnorm( n = size0, sigma = Sigma0,method = "svd")
  class1 <- rmvnorm( n = size1, sigma = Sigma1,method = "svd" ) 

  edgeid <- list()
  cPre <- list()
  PLR <- list()
  for (i in 1:times){
    set.seed(seed2[i])   
    netData <- rbind(class0,class1)
    deltaR <- abs(cor(class0) - cor(class1))
    deltaR[lower.tri(deltaR)] <- 0
    EDGE <- which(deltaR> thres, arr.ind=T)
    # nrow(EDGE)
    edgeid[[i]] <- nrow(deltaR)*(EDGE[,2] - 1) + EDGE[,1]
    # sum(EDGEid %in% edgeid)

    class0 <- rmvnorm( n = size0, sigma = Sigma0,method = "svd")
    class1 <- rmvnorm( n = size1, sigma = Sigma1,method = "svd" ) 
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
