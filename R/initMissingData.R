###
###Missing Data
###
set.seed(seed)
if(MissData){
  Yaggreg <- array(0,dim(Y[,,1]))
  for(tt in 1:TT) Yaggreg <- Yaggreg + Y[,,tt]
  Yaggreg[which(Yaggreg>1,arr.ind=TRUE)] <- 1
  SPaggreg <- shortest.paths(graph=graph.adjacency(Yaggreg),mode="all")
  SPaggreg[which(SPaggreg==Inf,arr.ind=TRUE)] <- 5
  
  outDeg <- inDeg <- denom <- numeric(n)
  for(tt in 1:TT){
    denom[c(1:n)[-Missing[[tt]]]] <- denom[c(1:n)[-Missing[[tt]]]] + 1
    if(length(Missing[[tt]])==0) denom <- denom + 1
    inDeg <- inDeg + colSums(Y[,,tt])/(n-1-length(Missing[[tt]])+as.numeric(1:n%in%Missing[[tt]]))*(n-1)
    outDeg <- outDeg + rowSums(Y[,,tt])
  }
  inDeg= inDeg/TT
  outDeg <- round(outDeg/denom)
  
  for(tt in 1:TT){
    for(i in Missing[[tt]]){
      Probs <- inDeg[-i]/SPaggreg[i,-i]
      ind <- sample(size=outDeg[i],x=c(1:n)[-i],prob=Probs,replace=FALSE)
      Y[ind,i,tt] <- Y[i,ind,tt]  <- 1
    }
  }
}