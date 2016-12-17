#' Initial betaIn and betaOUT; scale latent positions
#' 
#' @param Y an n x n x T array of relational matrices, 
#' where the third dimension corresponds to different time periods.
#' @param xLatPos Latent space position of actors
#' @param p number of latent dimensions
#' @param w weights or radii
#' @usage gmds( Y ) 
#' @return returns list:
#' \item{betaInInit}{Initial value for $\beta_{IN}$}
#' \item{betaOutInit}{Initial value for $\beta_{OUT}$}
#' \item{xLatPos}{Scaled latent space}
#' @export
#' 

initBetaInOut <- function(Y, xLatPos, p, w){

	#
	n <- dim(Y)[1]
	T <- dim(Y)[3]		

	# wrap c++ fns
	initializeWrap <- function(x){ 
		-initializeBinom(xLatPos,c(n,p,T),Y,Xscale=1/n, BIN=x[1],BOUT=x[2],w[,1]) }
	initializeGradWrap <- function(x){ 
		-initializeBinomGrad(xLatPos,c(n,p,T),Y,Xscale=1/n, BIN=x[1],BOUT=x[2],w[,1])[2:3] }

	# optimize get out results
	optimInit <- optim(par=c(1,1), fn=initializeWrap, gr=initializeGradWrap, method="BFGS")
	xLatPos <- xLatPos/n
	betaInInit <- max(optimInit$par[1],1e-4)
	betaOutInit <- max(optimInit$par[2],1e-4)

	# 
	return( list( betaInInit=betaInInit, betaOutInit=betaOutInit, xLatPos=xLatPos ) )
}