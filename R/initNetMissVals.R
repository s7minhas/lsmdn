#' Initial values for missing entries
#' 
#' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
#' @param missing List with information on missingness
#' @param seed random seed
#' @usage initNetMissVals( Y ) 
#' @return Inputted Y network with missingness filled in
#' @export getStartingValues

initNetMissVals <- function(Y, missing, seed=6886){
	
	#
	n <- dim(Y)[1]
	T <- dim(Y)[3]	
	yAggreg <- apply(Y, c(1,2), sum)
	
	# 
	yAggreg[which(yAggreg>1), arr.ind=TRUE] <- 1

	# 
	yGraph <- igraph::graph.adjacency(yAggreg)
	spAggreg <- igraph::shortest.paths(graph=yGraph, mode="all")
	spAggreg[which(spAggreg==Inf, arr.ind=TRUE)] <- 5

	#
	outDeg <- inDeg <- denom <- numeric(n)
	for(t in 1:T){
		denom[c(1:n)[-missing[[t]]]] <- denom[c(1:n)[-missing[[t]]]] + 1
		if(length(missing[[t]])==0) denom <- denom + 1
		inDeg <- inDeg + colSums(Y[,,t])/(n-1-length(missing[[t]])+as.numeric(1:n%in%missing[[t]]))*(n-1)
		outDeg <- outDeg + rowSums(Y[,,t])
	}
	inDeg= inDeg/T
	outDeg = outDeg/T ### the later seems it should only be executed if the actor is missing for the entire year
	# outDeg <- round(outDeg/denom)

	# 
	for(t in 1:T){
		for(i in missing[[t]]){
			probs <- inDeg[-i]/spAggreg[i,-i]
			ind <- sample(size=outDeg[i],x=c(1:n)[-i],prob=probs,replace=FALSE)
			Y[ind,i,t] <- Y[i,ind,t]  <- 1
		}
	}

	# 
	return(Y)	
}