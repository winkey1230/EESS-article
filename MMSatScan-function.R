#' @description This function is to implement the estimation-error-based scan \r
#' statistic (EESS).
#' @author Wei Wang
#' @param id A vector: the identification numbers for regions.
#' @param geo A dataframe or matrix including two columns referring to plane coordinates.
#' @param value A matrix: each row defines a estimated region-specific ERR.
#' @param slist A list of matrix: each element indicates the covariance of the corresponding ERR.
#' @param maxsize A percent: the maximum window size. 
#' @param mcmc An integer indicating the resampling number of Monte carlo test.
#' @param replica Logic: indicating whether the replicate clusters are reported.
#' @param ncore A integer referring to the number of cores for parallel computation.
MMSatScan <- function(id,geo,value,slist,maxsize = NULL,mcmc = 99,replica = F,ncore = 4){
  if(is.null(maxsize)) maxsize <- 0.5
  n <- length(id)
  maxn <- floor(n*maxsize)
  ins <- lapply(slist, function(x) chol2inv(chol(x))) # the inverse matrix of s
  dets <- unlist(lapply(slist, function(x) determinant(x)[[1]])) # ln(s)
  metabeta <- function(a,blist){
    # a: beta  blist: inverse cov
    covbeta0 <- chol2inv(chol(Reduce("+",blist))) # inverse of the sum of ins
    sbeta0 <- lapply(1:length(blist),function(x) blist[[x]]%*%a[x,]) # s * beta
    beta0 <- (covbeta0 %*% Reduce("+",sbeta0))[,1]
    beta0
  }
  llfun <- function(a,blist){
    aa <- sapply(1:length(blist), function(x) a[x,]%*%blist[[x]]%*%a[x,])
    -sum(aa)
  }
  
  distance <- as.matrix(dist(geo,diag = T,method = "euclidean"))
  idname <- row.names(distance)
  
  # find the maxll and cluster with replicate units
  clusterlist <- list()
  maxlli <- NULL
  ll0 <- llfun(value,ins)
  for (j in 1:n) {
    orderj <- sort(distance[j,])
    idj <- match(names(orderj),idname)
    maxllj <- NULL
    for (k in 1:maxn) {
      id1 <- idj[1:k] # id in cluster
      id2 <- idj[-(1:k)] # id not in cluster
      if(k==1) ll1 <- 0 else{
        beta1 <- metabeta(value[id1,],blist = ins[id1])
        betaeta1 <- t(t(value[id1,]) - beta1)
        ll1 <- llfun(betaeta1,ins[id1])
      }
      beta2 <- metabeta(value[id2,],blist = ins[id2])
      betaeta2 <- t(t(value[id2,]) - beta2)
      ll2 <- llfun(betaeta2,ins[id2])
      maxllj <- c(maxllj,ll1+ll2)
    }
    clusterlist[[j]] <- id[idj[1:which.max(maxllj)]]  
    maxlli <- c(maxlli,max(maxllj)-ll0)
  }
  orderid <- order(maxlli,decreasing = T) 
  maxlli <- maxlli[orderid]
  clusterlist <- clusterlist[orderid]
  
    
  # mcmc for ll
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  #clusterExport(cl = cl,c("mvrnorm"))
  pb <- txtProgressBar(max = mcmc, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  xx <- foreach(i = 1:mcmc,.combine = "c",.multicombine = T,
                .options.snow = opts) %dopar%{
                  set.seed(i)
                  repermutaion <- sample(n,n,replace = F)
                  betasim <- value[repermutaion,]
                  pins <- ins[repermutaion]
                  maxllsim <- -Inf
                  for (j in 1:n) {
                    orderj <- sort(distance[j,])
                    idj <- match(names(orderj),idname)  
                    for (k in 1:maxn) {
                      id1 <- idj[1:k] # id in cluster
                      id2 <- idj[-(1:k)] # id not in cluster
                      if(k==1) ll1 <- 0 else{
                        beta1 <- metabeta(betasim[id1,],blist = pins[id1])
                        betaeta1 <- t(t(betasim[id1,]) - beta1)
                        ll1 <- llfun(betaeta1,pins[id1])
                      }
                      beta2 <- metabeta(betasim[id2,],blist = pins[id2])
                      betaeta2 <- t(t(betasim[id2,]) - beta2)
                      ll2 <- llfun(betaeta2,pins[id2])
                      maxllsim <- max(maxllsim,ll1+ll2-ll0)
                    }
                  }
                  maxllsim
                }
  close(pb)
  stopCluster(cl)
  
  pvalue <- sapply(maxlli, function(x) mean(x<=c(xx,x)))
  clusterlist <- lapply(1:n, function(x) {
    attr(clusterlist[[x]],"p_value") <- pvalue[x]
    attr(clusterlist[[x]],"ll_ratio") <- maxlli[x]
    clusterlist[[x]]
  }) 
  if(!replica) {
    cluster <- clusterlist[1]
    for (i in 2:n) {
      aa <- intersect(unlist(cluster),clusterlist[[i]])
      if(length(aa)==0) cluster <- c(cluster,clusterlist[i]) 
    }
  }
  clusters <- NULL
  for (i in 1:length(cluster)) {
    bb <- cluster[[i]]
    aa <- cbind(i,bb,attr(bb,"p_value"),attr(bb,"ll_ratio"))
    clusters <- rbind(clusters,aa)
  }
  colnames(clusters) <- c("i","id","pvalue","llratio")
  list(clusters = clusters,maxllsim = xx)
}

library(doSNOW)
library(MASS)




