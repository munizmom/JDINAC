library(mvtnorm) # data generate
library(Matrix) # data generate
library(KernSmooth)
library(akima)


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



classPreNB <- function(EDGE,classLabel,DataFit,DataPre,nsplit=1,nfolds=10){
  
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
      yp <- rowSums(denX[-(1:length(y)),])
      preY <- cbind(preY,yp) 
    } 
    yPre <- rowMeans(preY) 
  } 
  list(yPre=yPre) 
}



NBPreSimu3 <- function(size0=300,size1=300,nP=100,Pdim=100,lim=2,s=4/3,e=.5,times=100,thres=0,nsplit=1,nfolds=10,seed=100){

  set.seed(seed)
  seed2 <- sample(100:2016000,times)
  class0 <- array(runif(size0*nP,-lim,lim),dim=c(size0,nP))
  class0[,(Pdim/2+1):Pdim] <- class0[,(Pdim/2+1):Pdim]^2 + rmvnorm(n = size0,sigma=diag(e^2,Pdim/2)) - s
  class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
  cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
  class1[,1:(Pdim/2)] <- cla
  class1[,(Pdim/2+1):Pdim] <- cla^2 + rmvnorm(n = size1,sigma=diag(e^2,Pdim/2)) - s
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
    class0[,(Pdim/2+1):Pdim] <- class0[,(Pdim/2+1):Pdim]^2 + rmvnorm(n = size0,sigma=diag(e^2,Pdim/2)) - s
    class1 <- array(runif(size1*nP,-lim,lim),dim=c(size1,nP))
    cla <- array(runif(size1*Pdim/2,-lim,lim),dim=c(size1,Pdim/2))
    class1[,1:(Pdim/2)] <- cla
    class1[,(Pdim/2+1):Pdim] <- cla^2 + rmvnorm(n = size1,sigma=diag(e^2,Pdim/2)) - s
    DataPre <- rbind(class0,class1)
    
    cPre[[i]] <- classPreNB(EDGE=EDGE,classLabel=ylabel,DataFit=netData,DataPre=DataPre)
    
}
 list(cPre=cPre)
}
