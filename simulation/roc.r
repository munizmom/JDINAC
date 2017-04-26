library(glmnet)

roc_cPre <- function(obj,size0,size1,cutn){
  y <- c(rep(0,size0),rep(1,size1))
  auc <- c(0)
  times <- length(obj$cPre)
  sen <- array(0,dim=c(times,cutn))
  spe <- array(0,dim=c(times,cutn))
  err <- array(0,dim=c(times,cutn))
  maxY <- ceiling(max(unlist(lapply(obj$cPre,function(x){max(abs(x$yPre))}))))
  for(i in 1:times){
    yPre <- obj$cPre[[i]]$yPre
    if(max(yPre)>1) yPre <- yPre/(2*maxY) + 0.5
    auc[i] <- auc(y,yPre)
    cutoff <- seq(0, 1, length.out = cutn)
    for(j in 1:cutn){
      sen[i,j] <- sum((y==1) & (yPre>=cutoff[j]))/sum(y==1)
      spe[i,j] <- sum((y==0) & (yPre<cutoff[j]))/sum(y==0)
      err[i,j] <- mean(((yPre > cutoff[j]) - y)^2)
    }
  }
  AUC <- c(mean(auc),sd(auc))
  names(AUC) <- c("auc","sd")
  SEN <- colMeans(sen)  #Sensitivity
  SPE <- colMeans(spe)  #Specificity
  ERR <- colMeans(err)  #classification error
  ERRsd <- apply(err,2,sd)
  list(AUC=AUC,SEN=SEN,SPE=SPE,ERR=ERR,ERRsd=ERRsd)
}

roc_PLR <- function(obj,size0,size1,cutn){
  y <- c(rep(0,size0),rep(1,size1))
  auc <- c(0)
  times <- length(obj$PLR)
  sen <- array(0,dim=c(times,cutn))
  spe <- array(0,dim=c(times,cutn))
  err <- array(0,dim=c(times,cutn))
  maxY <- ceiling(max(unlist(lapply(obj$cPre,function(x){max(abs(x$yPre))}))))
  for(i in 1:times){
    yPre <- obj$PLR[[i]]$yPre
    if(max(yPre)>1) yPre <- yPre/(2*maxY) + 0.5
    auc[i] <- auc(y,obj$PLR[[i]]$yPre)
    cutoff <- seq(0, 1, length.out = cutn)
    for(j in 1:cutn){
      sen[i,j] <- sum((y==1) & (yPre>=cutoff[j]))/sum(y==1)
      spe[i,j] <- sum((y==0) & (yPre<cutoff[j]))/sum(y==0)
      err[i,j] <- mean(((yPre > cutoff[j]) - y)^2)
    }
  }
  AUC <- c(mean(auc),sd(auc))
  names(AUC) <- c("auc","sd")
  SEN <- colMeans(sen)  #Sensitivity
  SPE <- colMeans(spe)  #Specificity
  ERR <- colMeans(err)  #classification error
  ERRsd <- apply(err,2,sd)
  list(AUC=AUC,SEN=SEN,SPE=SPE,ERR=ERR,ERRsd=ERRsd)
}

optCut <- function(roc1,roc2,roc3,roc4,roc5,cutn){

  cutoff <- seq(0, 1, length.out = cutn)

  # optimal cutoff by maximum Sensitivity + Specificity
  maxSScut <- cutoff[c(which.max(roc1$SEN + roc1$SPE),which.max(roc2$SEN + roc2$SPE),
  which.max(roc3$SEN + roc3$SPE),which.max(roc4$SEN + roc4$SPE),which.max(roc5$SEN + roc5$SPE))]
  maxSScut

  # optimal cutoff by minimum classification error
  minEcut <- cutoff[c(which.min(roc1$ERR),which.min(roc2$ERR),which.min(roc3$ERR),
  which.min(roc4$ERR),which.min(roc5$ERR))]
  minEcut

  optcut <- cbind(maxSScut,minEcut)
  rownames(optcut) <- c("JDINAC","cPLR","RF","NB","oPLR")
  optcut
}
