#' Initial \beta_{IN} and \beta_{OUT}; scale latent positions
#' 
#' @param Y an n x n x T array of relational matrices, 
#' where the third dimension corresponds to different time periods.
#' @param xLatPos Latent space position of actors
#' @param p number of latent dimensions
#' @param w weights
#' @usage gmds( Y ) 
#' @return returns list:
#' \item{betaInInit}{Initial value for beta in}
#' \item{betaOutInit}{Initial value for beta out}
#' \item{xLatPos}{Scaled latent space}
#' @export getStartingValues

initBetaInOut <- function(Y, xLatPos, p, w){

	#
	n <- dim(Y)[1]
	T <- dim(Y)[3]		

	# wrap c++ fns
	initialize.wrap <- function(x){ 
		-cInitialize1(xLatPos,c(n,p,T),Y,Xscale=1/n, BIN=x[1],BOUT=x[2],w[,1]) }
	initialize.grad.wrap <- function(x){ 
		-cInitialize1Grad(xLatPos,c(n,p,T),Y,Xscale=1/n, BIN=x[1],BOUT=x[2],w[,1])[2:3] }

	# optimize get out results
	optimInit <- optim(par=c(1,1),fn=initialize.wrap, gr=initialize.grad.wrap,method="BFGS")
	xLatPos <- xLatPos/n
	betaInInit <- max(optimInit$par[1],1e-4)
	betaOutInit <- max(optimInit$par[2],1e-4)

	# 
	return( list( betaInInit=betaInInit, betaOutInit=betaOutInit, xLatPos=xLatPos ) )
}