#' Log-likelihood approximation subsampling
#' 
#' @param Y an n x n x T array of relational matrices, 
#' where the third dimension corresponds to different time periods.
#' @param p number of latent dimensions
#' @param family type of model to run. Options include 'normal', 'nonNegNormal', 'poisson', 'binomial'. 
#' @param llApprox logical indicating whether or not to utilize log-likelihood 
#' approximationg. Only available for binomial model types.
#' @param missData logical indicating whether to impute missing data.
#' @param N number of MCMC iterations.
#' @param seed random seed
#' @usage getStartingValues( Y, p=2, family='binomial', llApprox=FALSE, missData=FALSE, N=1000, seed=6886) 
#' @return returns list of starting values:
#' \item{w}{weights}
#' \item{X}{initial actor latent space positions calculated via GMDS}
#' \item{betaIn}
#' \item{betaOut}
#' \item{nuIn}
#' \item{nuOut}
#' \item{xiIn}
#' \item{t2}
#' \item{shapeT2}
#' \item{scaleT2}
#' \item{s2}
#' \item{shapeS2}
#' \item{scaleS2}
#' if llApprox=TRUE, also returns
#' \item{dInMax}
#' \item{dOutMax}
#' \item{n0}
#' \item{elOut}
#' \item{elIn}
#' \item{subseq}
#' \item{degree}
#' \item{edgeList}
#' @export lsmdn
#' 

getStartingValues <- function(Y, p=2, family='binomial', llAprox=FALSE, missData=FALSE, N=1000, seed=6886){

	# Starting values for parameters
	n <- dim(Y)[1]
	T <- dim(Y)[3]		
	betaOut <- betaIn <- numeric(N)
	s2 <- t2 <- numeric(N)
	AccRate <- numeric(3+n*T)
	names(AccRate) <- c("betaIn","betaOut","weights", paste("X",rep(1:n,T),rep(1:T,each=n),sep=","))

	# init values for missingness in Y
	tmp = which(is.na(Y),arr.ind=TRUE)
	missing <- lapply( 1:T, function(t){ return( unique(tmp[which(tmp[,3]==t),1]) ) } )
	Y[tmp] = 0; rm(tmp)			
	if(missData){ Y <- initNetMissVals( Y, missing ) }

	# Weights
	w <- initWeights(Y, N)

	# Initial latent space positions via Sarkar & Moore 2005
	X <- gmds( Y )

	# initial values for betaIn and betaOut
	tmp <- initBetaInOut(Y, xLatPos=X[[1]], p=p, w=w)
	X[[1]] <- tmp$'xLatPos'
	betaIn[1] <- tmp$'betaInInit'
	betaOut[1] <- tmp$'betaOutInit'
	rm(tmp)

	# Priors and further initializations
	nuIn <- betaIn[1]
	nuOut <- betaOut[1]
	xiIn <- xiOut <- 100

	#Tau^2
	t2[1] <- sum(X[[1]][,1,]*X[[1]][,1,])/(n*p)
	shapeT2 <- 2.05
	scaleT2 <- (shapeT2-1)*t2[1]

	#Sigma^2
	s2[1] <- 0.001
	shapeS2 <- 9
	scaleS2 <- 1.5	

	# init values for log-likelihood approximation subsampling
	if(llAprox){
		tmp <- initLogLikeApprox(Y, n0, seed)
		return( 
			list(
				w=w, X=X, betaIn=betaIn, betaOut=betaOut, nuIn=nuIn, nuOut=nuOut,
				xiIn=xiIn, xiOut=xiOut, t2=t2, shapeT2=shapeT2, scaleT2=scaleT2,
				s2=s2, shapeS2=shapeS2, scaleS2=scaleS2, 
				dInMax=tmp$dInMax, dOutMax=tmp$dOutMax, n0=tmp$n0, 
				elOut=tmp$elOut, elIn=tmp$elIn, subseq=tmp$subseq, degree=tmp$degree,
				edgeList=tmp$edgeList
				)
			)
	}

	#
	return(
		list(
			w=w, X=X, betaIn=betaIn, betaOut=betaOut, nuIn=nuIn, nuOut=nuOut,
			xiIn=xiIn, xiOut=xiOut, t2=t2, shapeT2=shapeT2, scaleT2=scaleT2,
			s2=s2, shapeS2=shapeS2, scaleS2=scaleS2, 
			dInMax=tmp$dInMax, dOutMax=tmp$dOutMax, n0=tmp$n0			
			)
		)
}

