#' Log-likelihood approximation subsampling
#' 
#' @param Y an n x n x T array of relational matrices, 
#' where the third dimension corresponds to different time periods.
#' @param n0 n0
#' @param seed random seed
#' @param N number of MCMC iterations
#' @usage initWeights( Y, N ) 
#' @return returns list of starting values for log likelihood approximation:
#' \item{dInMax}
#' \item{dOutMax}
#' \item{n0}
#' \item{elOut}
#' \item{elIn}
#' \item{subseq}
#' \item{degree}
#' \item{edgeList}
#' @export getStartingValues
#' 

initLogLikeApprox <- function(Y, n0, seed){

	#
	n <- dim(Y)[1]
	T <- dim(Y)[3]		

	degree <- array(0,c(n,2,T))
	for(t in 1:T){
		degree[,1,t] <- colSums(Y[,,t])
		degree[,2,t] <- rowSums(Y[,,t])
	}
	dInMax <- max(degree[,1,])
	dOutMax <- max(degree[,2,])
	elIn <- array(0,c(n,dInMax,T))
	elOut <- array(0,c(n,dOutMax,T))

	#
	for(t in 1:T){
		for(i in 1:n){
			ind <- which(Y[,i,t]==1)
			if(length(ind)>0){ elIn[i,1:length(ind),t] <- ind }
			ind <- which(Y[i,,t]==1)
			if(length(ind)>0){ elOut[i,1:length(ind),t] <- ind }
		}
	}
	n0<-max(n0,dInMax+10,dOutMax+10)

	#
	edgeList <- lapply(1:n, function(i){
		uniqueEdges <- unique( which(Y[i,,]==1|Y[,i,]==1,arr.ind=TRUE)[,1] )
		return(uniqueEdges) })

	#
	subseq <- matrix(0,n,n0)
	for(i in 1:n){
		nOnes <- round(length(edgeList[[i]])/n*n0) #stratified sampling
		if(length(edgeList[[i]])>0){
			nOnes <- max(nOnes,1)
			set.seed(seed) ; subseq[i,1:nOnes] <- sample(edgeList[[i]],size=nOnes,replace=TRUE) # should replace be false?      
		}
		set.seed(seed) ; subseq[i,(nOnes+1):n0] <- sample(c(1:n)[-c(i,edgeList[[i]])],size=n0-nOnes,replace=TRUE)
	}

	return( list(
		dInMax=dInMax, dOutMax=dOutMax, n0=n0, elOut=elOut, 
		elIn=elIn, subseq=subseq, degree=degree, edgeList=edgeList
		) )
}


