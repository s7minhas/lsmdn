#' Initial parameter values for lsmdn
#' 
#' @param Y an n x n x T array of relational matrices, 
#' where the third dimension corresponds to different time periods.
#' @param W an n x n x T x p design array, where the third dimension corresponds to
#' time and the fourth to different covariates.
#' @param p number of latent dimensions
#' @param family type of model to run. Options include 'normal', 'nonNegNormal', 'poisson', 'binomial'. 
#' @param llApprox logical indicating whether or not to utilize log-likelihood 
#' approximationg. Only available for binomial model types.
#' @param missData logical indicating whether to impute missing data.
#' @param s2Init starting value for s2
#' @param t2Init starting value for t2
#' @param xLatPos starting actor positions in latent space
#' @param betaInInit starting value for betaIn
#' @param betaOutInit starting value for betaOut
#' @param nuIn prior mean for betaIn
#' @param nuOut prior variance for betaOut
#' @param xiIn prior variance for betaIn
#' @param xiOut prior variance for betaOut
#' @param shapeT2 shape parameter for t2
#' @param scaleT2 scale parameter for t2
#' @param shapeS2 shape parameter for s2
#' @param scaleS2 shape parameter for s2
#' @param N number of MCMC iterations.
#' @param seed random seed
#' @usage getStartingValues( Y, p=2, family='binomial', llApprox=FALSE, missData=FALSE, N=1000, seed=6886) 
#' @return returns list of starting values:
#' \item{w}{weights}
#' \item{X}{initial actor latent space positions calculated via GMDS}
#' \item{betaIn}{Starting value for importance of popularity}
#' \item{betaOut}{Relative value for importance of activity }
#' \item{nuIn}{Prior mean for beta-in}
#' \item{nuOut}{Prior mean for beta-out}
#' \item{xiIn}{Prior variance for beta-in}
#' \item{xiOut}{Prior varianerce for beta-out}
#' \item{t2}{Variance for latent positions in initial time period}
#' \item{shapeT2}{Shape parameter for prior distribution of tau2}
#' \item{scaleT2}{Scale parameter for prior distribution of tau2}
#' \item{s2}{Variance of movement in the latent space in subsequent periods}
#' \item{shapeS2}{Shape parameter for prior distribution of sigma2}
#' \item{scaleS2}{Scale parameter for prior distribution of sigma2}
#' if llApprox=TRUE, also returns
#' \item{dInMax}{add desc}
#' \item{dOutMax}{add desc}
#' \item{n0}{add desc}
#' \item{elOut}{add desc}
#' \item{elIn}{add desc}
#' \item{degree}{add desc}
#' \item{edgeList}{add desc}
#' @export
#' 

getStartingValues <- function(
	Y, W, p, family, llApprox, missData, N, seed,
	s2Init=NULL, t2Init=NULL, xLatPos=NULL, betaInInit=NULL, betaOutInit=NULL,
	nuIn=NULL, nuOut=NULL, xiIn=NULL, xiOut=NULL, shapeT2=NULL, scaleT2=NULL, 
	shapeS2=NULL, scaleS2=NULL, g2=NULL, shapeG2=NULL, scaleG2=NULL, sdLambda = NULL
	){

	# Starting values for parameters
	n <- dim(Y)[1]
	T <- dim(Y)[3]		
	betaOut <- betaIn <- numeric(N)
	s2 <- t2 <- numeric(N)
	accRate <- numeric(3+n*T)
	names(accRate) <- c("betaIn","betaOut","weights", paste("X",rep(1:n,T),rep(1:T,each=n),sep=","))

	# init values for missingness in Y
	if(missData){
		tmp = which(is.na(Y),arr.ind=TRUE)
		missing <- lapply( 1:T, function(t){ return( unique(tmp[which(tmp[,3]==t),1]) ) } )
		Y[tmp] = 0; rm(tmp)					
		Y <- initNetMissVals( Y, missing )
		Y <- Y
	}

	# set up lambda obj for covars
	lambda <- matrix(0, ncol = 1, nrow = dim(W)[4])
	sdLambda <- 50

	# Weights
	w <- initWeights(Y, N)

	# Initial latent space positions via Sarkar & Moore 2005
	X <- gmds( Y, w, p, family )

	# initial values for betaIn and betaOut
	if( family=='binomial' ){
		tmp <- initBetaInOut(Y, xLatPos=X[[1]], p=p, w=w)
		X[[1]] <- tmp$'xLatPos'
		betaIn[1] <- tmp$'betaInInit'
		betaOut[1] <- tmp$'betaOutInit'
		rm(tmp)
	} else { # could condition this on dist
		X[[1]] <- X[[1]]/n
		betaIn[1] <- 10
		betaOut[1] <- 10
	}
	alpha = (betaIn[1] + betaOut[1])*(1 - 2*rbinom(1,1,.5))
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

	#Gamma^2
	if(family=='nonNegNormal'){
		g2 <- numeric(N)
		g2[1] <- 25
		shapeG2 <- 2.05
		scaleG2 <- (1.05)*25 }

	# init values for log-likelihood approximation subsampling
	if( llApprox & family=='binomial' ){
		tmp <- initLogLikeApprox(Y, n0, seed)
		return( 
			list(
				Y=Y, w=w, X=X, betaIn=betaIn, betaOut=betaOut, nuIn=nuIn, nuOut=nuOut,
				xiIn=xiIn, xiOut=xiOut, t2=t2, shapeT2=shapeT2, scaleT2=scaleT2,
				s2=s2, shapeS2=shapeS2, scaleS2=scaleS2, 
				dInMax=tmp$dInMax, dOutMax=tmp$dOutMax, n0=tmp$n0, 
				elOut=tmp$elOut, elIn=tmp$elIn, degree=tmp$degree,
				edgeList=tmp$edgeList, accRate=tmp$accRate, alpha = alpha
				)
			)
	} else { # could condition this on dist
		return(
			list(
				Y=Y, w=w, X=X, betaIn=betaIn, betaOut=betaOut, nuIn=nuIn, nuOut=nuOut,
				xiIn=xiIn, xiOut=xiOut,
				t2=t2, shapeT2=shapeT2, scaleT2=scaleT2,
				s2=s2, shapeS2=shapeS2, scaleS2=scaleS2, 
				g2=g2, shapeG2=shapeG2, scaleG2=scaleG2, 
				accRate=accRate, lambda = lambda, sdLambda = sdLambda, alpha = alpha
				)
			)
	}

}

