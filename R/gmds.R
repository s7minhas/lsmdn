#' Initial Latent Positions (GMDS, Sarkar and Moore, 2005)
#' 
#' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
#' @param w weights
#' @param p number of latent dimensions
#' @param family type of model to run. Options include 'normal', 'nonNegNormal', 'poisson', 'binomial'. 
#' @usage gmds( Y ) 
#' @return Latent positions of actors based on Sarkar and Moore
#' @export
#' 

###
gmds <- function(Y, w, p, family){
	#
	n <- dim(Y)[1]
	T <- dim(Y)[3]	

	dissim <- array(0,dim=dim(Y))
	if(family=='binomial'){
		for(t in 1:T){ dissim[,,t] <- igraph::shortest.paths(graph=igraph::graph.adjacency(Y[,,t]),mode="all") }
		dissim[which(dissim==Inf,arr.ind=TRUE)] <- max(c(dissim[which(dissim!=Inf,arr.ind=TRUE)]))
	} else {
		Ybin = Y > 0 
		dsim = function(i,j,t, Y, w ){
			con = Y[i,j,t] + Y[j,i,t]
			if(con == 2){ret = min(w[i], w[j])/2}
			if(con == 1){ret = (w[i] + w[j])/2}
			if(con == 0){ret = 3*(w[i] + w[j])/2}
			return(ret)	
		}

		i = 1:dim(Y)[1] ; j = 1:dim(Y)[2]
		for(t in 1:T){
			dissim[,,t] = sapply(i, function(i){ sapply(j, dsim, i = i, t = t, Y = Ybin, w = w[,1]) } )
			diag(dissim[,,t]) = 0 }
	}
	
	X <- list()
	X[[1]] <- array(0,c(p,T,n))
	X[[1]][,1,] <- t(cmdscale(d=dissim[,,1],k=p))
	temp.lambda <- 10
	H <- matrix(-1/n,n,n)+diag(n)

	for(t in 2:T){
		temp <- 1/(1+temp.lambda)*H%*%(-1/2*dissim[,,t]^2)%*%H + 
			temp.lambda/(1+temp.lambda)*t(X[[1]][,t-1,])%*%X[[1]][,t-1,]
		temp <- eigen(temp)
		X[[1]][,t,] <- t(temp$vectors[,1:p]%*%diag(temp$values[1:p])^(1/2))
		X[[1]][,t,] <- t(vegan::procrustes(X=t(X[[1]][,t-1,]),Y=t(X[[1]][,t,]),scale=FALSE)$Yrot)
	}

	return(X)
}
