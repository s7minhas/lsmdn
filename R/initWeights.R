#' Initial weights
#' 
#' @param Y an n x n x T array of relational matrices, 
#' where the third dimension corresponds to different time periods.
#' @param N number of MCMC iterations
#' @usage initWeights( Y, N ) 
#' @return returns n x N matrix of weights
#' @export getStartingValues
#' 

initWeights <- function(Y, N){

	#
	n <- dim(Y)[1]
	T <- dim(Y)[3]		

	#
	w <- matrix(0,n,N)	
	for(t in 1:T){
		w[,1] <- w[,1] + apply(Y[,,t],1,sum) + apply(Y[,,t],2,sum)
	}
	w[,1] <- w[,1]/sum(Y)/2

	#
	if(sum(w==0)>0){
		w[,1] <- w[,1]+1e-5
		w[,1] <- w[,1]/sum(w[,1])
	}

	#
	w[,1] <- w[,1]/sum(w[,1])

	#
	return(w)
}