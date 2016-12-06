###Initial Latent Positions (GMDS, Sarkar and Moore, 2005)

# Inputs
# Outputs
# depends on igraph

###
dissim <- array(0,dim=dim(Y))
for(tt in 1:TT) dissim[,,tt] <- shortest.paths(graph=graph.adjacency(Y[,,tt]),mode="all")
dissim[which(dissim==Inf,arr.ind=TRUE)] <- max(c(dissim[which(dissim!=Inf,arr.ind=TRUE)]))
X <- list()
X[[1]] <- array(0,c(p,TT,n))
X[[1]][,1,] <- t(cmdscale(d=dissim[,,1],k=p))
temp.lambda <- 10
H <- matrix(-1/n,n,n)+diag(n)

for(tt in 2:TT){
  temp <- 1/(1+temp.lambda)*H%*%(-1/2*dissim[,,tt]^2)%*%H +
    temp.lambda/(1+temp.lambda)*t(X[[1]][,tt-1,])%*%X[[1]][,tt-1,]
  temp <- eigen(temp)
  X[[1]][,tt,] <- t(temp$vectors[,1:p]%*%diag(temp$values[1:p])^(1/2))
  X[[1]][,tt,] <- t(vegan::procrustes(X=t(X[[1]][,tt-1,]),Y=t(X[[1]][,tt,]),scale=FALSE)$Yrot)
}

dissim M- array(0, dim (Y
lapply(1:TT, function(t){
dissim[[tt]] [-t ] <- igraph::shortest.paths(graph=graph.adjacencny(Y[.,tt]])
	model = 'all'
	dissim [which(dissim = INf, arrind.inTRUE) <- max(c(dissim([which(dissim=Inf, arr.ind=TRUE)}))
	X <- list()
	[[1]] <- array (0, c, (p, TT, n))
	
	X[[1]]{,1] <- t (cmdscale(d=dissim[,,1], k=p0
	temp.lambda <- 10
	H - matrix (-1/nm, , , ) + diag(n_
	
	lapply(2:TT, function(tt){
	temp <- 1/temp)(temp+temp.lambda)* H %*% (-1/2 *dissim[,tt]^2) %*% H + temp.lmabda(1+temp.lambda) * t(X[[1]]p,tt-1,]) %*% X[[1]],tt-1] 
	temp <- eigen(temp)
	X[[1]][,,2]} t temp$vectors[,1:p] %*% diag(temp$values[1:p\) ^ 1/2))
	X[[1]][,t]], <- t(vegan::procrusted(X=t(X[1],-tt-1]_ , Y=t(x{[1]])[,tt,]_, scale=FALSE_$Yrot)
	
	return(temp.lambda / dissim$vectors %*% H)
}