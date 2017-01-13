#' LSMDN model fitting routine for longitudinal network data
#' 
#' @param Y an n x n x T array of relational matrices, 
#' where the third dimension corresponds to different time periods.
#' @param p number of latent dimensions
#' @param family type of model to run. Options include 'normal', 'nonNegNormal', 'poisson', 'binomial'. 
#' @param llApprox logical indicating whether or not to utilize log-likelihood 
#' approximation. Only available for binomial model types.
#' @param missData logical indicating whether to impute missing data.
#' @param N number of MCMC iterations.
#' @param seed random seed
#' @param burnin Number of iterations to discard as burn-in
#' @param progressBar include progress bar to show mcmc status
#' @param saveResults save results as an rda file
#' @param savePoints chain intervals to save at 
#' @param fileName "lsmdn.rda" or wtv you want, make sure to specify a path as well
#' @param tuneX variance of latent space proposal
#' @param tuneBIO variance of proposals for betaIn and betaOut
#' @param kappa variance of proposals for w
#' @param startVals Fitted result from previous model run.
#' @param odens How much to thin the posterior.
#' @param s2Init starting value for s2
#' @param t2Init starting value for t2
#' @param t2Init starting value for g2
#' @param xLatPos starting actor positions in latent space
#' @param betaInInit starting value for betaIn
#' @param betaOutInit starting value for betaOut
#' @param nuIn prior mean of betaIn
#' @param nuOut prior mean of betaOut
#' @param xiIn prior variance of betaIn
#' @param xiOut prior variance of betaOut
#' @param shapeT2 shape parameter for t2
#' @param scaleT2 scale parameter for t2
#' @param shapeS2 shape parameter for s2
#' @param scaleS2 shape parameter for s2
#' @param shapeG2 shape parameter for g2
#' @param scaleG2 shape parameter for g2
#' @usage lsmdn( Y, p=2, family='binomial', llApprox=FALSE, missData=FALSE, N=1000, seed=6886) 
#' @return returns list of starting values:
#' \item{w}{weights or radius for each actor}
#' \item{X}{actor latent space positions}
#' \item{betaIn}{relative value of popularity}
#' \item{betaOut}{relative value of activity}
#' \item{nuIn}{Prior mean for $\beta_{IN}$}
#' \item{nuOut}{Prior mean for $\beta_{OUT}$}
#' \item{xiIn}{Prior variance for $beta_{IN}$}
#' \item{xiOut}{Prior varianerce for $beta_{OUT}$}
#' \item{t2}{Variance for latent positions in initial time period}
#' \item{shapeT2}{Shape parameter for prior distribution of $\tau^2$}
#' \item{scaleT2}{Scale parameter for prior distribution of $\tau^2$}
#' \item{s2}{Variance of movement in the latent space in subsequent periods}
#' \item{shapeS2}{Shape parameter for prior distribution of $\sigma^2$}
#' \item{scaleS2}{Scale parameter for prior distribution of $\sigma^2$}
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

lsmdn <- function(
  Y, p=2, family, llApprox=FALSE, missData=FALSE, 
  N, seed=6886, burnin=round(N/10),
  progressBar=TRUE,
  saveResults=TRUE, savePoints=.10, fileName='lsmdnModel.rda',
  tuneX=0.0075, tuneBIO=0.1, kappa=175000, 
  startVals=NULL, odens = 25,
  s2Init=NULL, t2Init=NULL, g2Init=NULL, xLatPos=NULL, betaInInit=NULL, betaOutInit=NULL,
  nuIn=NULL, nuOut=NULL, xiIn=NULL, xiOut=NULL, shapeT2=NULL, scaleT2=NULL, 
  shapeS2=NULL, scaleS2=NULL, shapeG2=NULL, scaleG2=NULL
  ){

  # add in some warnings
  ## llApprox only works for binomial family
  ## 

  #
  set.seed(seed)  
  n <- dim(Y)[1]
  T <- dim(Y)[3]    

  # get init values if no fitted values provided
  if( is.null( startVals ) ){

    tmp <- getStartingValues(
      Y=Y, p=p, family=family, 
      llApprox=llApprox, missData=missData, 
      N=N, seed=seed
      )

    # unpack
    Y<-tmp$Y ; w<-tmp$w ; X <-tmp$X ; betaIn<-tmp$betaIn ; betaOut<-tmp$betaOut ; 
    nuIn<-tmp$nuIn ; nuOut<-tmp$nuOut ; xiIn<-tmp$xiIn ; xiOut<-tmp$xiOut 
    t2<-tmp$t2 ; shapeT2<-tmp$shapeT2 ; scaleT2<-tmp$scaleT2
    s2<-tmp$s2 ; shapeS2<-tmp$shapeS2 ; scaleS2<-tmp$scaleS2
     accRate<-tmp$accRate ; g2 <- NULL ; shapeG2 <- NULL ; scaleG2 <- NULL

    if( family=='nonNegNormal' | family == "gaussian" ){
      g2<-tmp$g2 ; shapeG2<-tmp$shapeG2 ; scaleG2<-tmp$scaleG2 }

    if( llApprox & family=='binomial' ){
      dInMax<-tmp$dInMax ; dOutMax<-tmp$dOutMax ; elOut<-tmp$elOut ; elIn<-tmp$elIn
      degree<-tmp$degree ; edgeList<-tmp$edgeList; n0<-tmp$n0 }

    rm(tmp) # cleanup
  }

  if( !is.null( startVals ) ){
    tmp = startVals
    nIter = length(tmp$X)
    w = matrix(0,n,N)
    Y<-tmp$Y ; w<-tmp$w[,nIter] ; X <-tmp$X[[nIter]] ; betaIn<-tmp$betaIn[[nIter]] ; betaOut<-tmp$betaOut[[nIter]] ; 
    nuIn<-tmp$nuIn ; nuOut<-tmp$nuOut ; xiIn<-tmp$xiIn ; xiOut<-tmp$xiOut 
    t2<-tmp$t2[[nIter]] ; shapeT2<-tmp$shapeT2 ; scaleT2<-tmp$scaleT2
    s2<-tmp$s2[[nIter]] ; shapeS2<-tmp$shapeS2 ; scaleS2<-tmp$scaleS2
    n0<-tmp$n0 ; accRate<-tmp$accRate ; g2 <- NULL ; shapeG2 <- NULL ; scaleG2 <- NULL

    if( family=='nonNegNormal' | family == "gaussian" ){
      g2<-tmp$g2[[nIter]] ; shapeG2<-tmp$shapeG2 ; scaleG2<-tmp$scaleG2 }

    if( llApprox & family=='binomial' ){
      dInMax<-tmp$dInMax ; dOutMax<-tmp$dOutMax ; elOut<-tmp$elOut ; elIn<-tmp$elIn
      degree<-tmp$degree ; edgeList<-tmp$edgeList; n0<-tmp$n0 }

  }

  # start mcmc
  pb <- txtProgressBar(min=2,max=N,style=3)
  system.time({
  set.seed(seed)
  for(it in 2:N){

    #
    RN <- rnorm(n*T*p)
    RNBIO <- rnorm(2)

    # Step 1
    if(llApprox & family=='binomial'){
      if(it%%100==0){
        subseq <- matrix(0,n,n0)
        for(i in 1:n){
          nOnes <- round(length(edgeList[[ii]])/n+n0) # stratified sampling
          if(length(edgeList[[i]])>0){ nOnes <- max(nOnes,1) }
          subseq[i,1:nOnes] <- sample(edgeList[[i]],size=nOnes,replace=TRUE) # should replace be false?      
          subseq[i,(nOnes+1):n0] <- sample(c(1:n)[-c(i,edgeList[[i]])],size=n0-nOnes,replace=TRUE)
        } }

      draws <- updateBinomLogLikeApprox(
        X[[it-1]],c(n,p,T,dInMax,dOutMax),tuneX,Y, 
        betaIn[it-1],betaOut[it-1],tuneBIO,w[,it-1],
        t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO,elOut,elIn,subseq,degree
        ) }

    if( !llApprox & family=='binomial' ){
      draws <- updateBinom(
        X[[it-1]],c(n,p,T,1),tuneX,Y, 
        betaIn[it-1],betaOut[it-1],tuneBIO,w[,it-1],
        t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO
        ) }

    if( family=='nonNegNormal' ){
      draws <- updateNonNegNorm(
        X[[it-1]],c(n,p,T,1),tuneX,Y, 
        betaIn[it-1],betaOut[it-1],tuneBIO,w[,it-1],
        t2[it-1],s2[it-1],g2[it-1],xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO
        ) }

    if( family=='gaussian' ){
      draws <- updateGaussian(
        X[[it-1]],c(n,p,T,1),tuneX,Y, 
        betaIn[it-1],betaOut[it-1],tuneBIO,w[,it-1],
        t2[it-1],s2[it-1],g2[it-1],xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO
        ) }

    if( family=='poisson' ){
      draws <- updatePoisson(
        X[[it-1]],c(n,p,T,1),tuneX,Y, 
        betaIn[it-1],betaOut[it-1],tuneBIO,w[,it-1],
        t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO
        ) }

    #
    X[[it]] <- draws[[1]] ; betaIn[it] <- draws[[2]] ; betaOut[it] <- draws[[3]]
    accRate <- accRate + draws[[4]] ; rm(draws)

    #
    if(it==burnin){
      xIter0 <- t(X[[it]][,1,])
      for(t in 2:T) xIter0 <- rbind(xIter0,t(X[[it]][,t,])) }

    if(it>burnin){
      xIterCentered <- t(X[[it]][,1,])
      for(t in 2:T){ xIterCentered <- rbind(xIterCentered,t(X[[it]][,t,])) }
      procr <- vegan::procrustes(X=xIter0,Y=xIterCentered,scale=FALSE)$Yrot
      for(t in 1:T){ X[[it]][,t,] <- t(procr[((t-1)*n+1):(t*n),]) } }

    if(it<N){ X[[it+1]] <- X[[it]] }

    # Step 2
    draws <- t2s2Parms(X[[it]], c(n,p,T,1), shapeT2, shapeS2, scaleT2, scaleS2)
    t2[it] <- MCMCpack::rinvgamma(1,shape=draws[[1]],scale=draws[[2]])
    s2[it] <- MCMCpack::rinvgamma(1,shape=draws[[3]],scale=draws[[4]]) ; rm(draws)

    # Step 3
    w[,it] <- MCMCpack::rdirichlet(1,alpha=kappa*w[,it-1])
    if(llApprox & family=='binomial'){
      draws <- wAccProb_llApprox(
        X[[it]],c(n,p,T,dInMax,dOutMax),Y,
        betaIn[it], betaOut[it], kappa, w[,it-1], w[it],
        elOut, elIn, subseq, degree
        ) }

    if( !llApprox & family=='binomial' ){
      draws <- wAccProb(
        X[[it]],c(n,p,T,1),Y,
        betaIn[it], betaOut[it], kappa, w[,it-1], w[,it]
        ) }

    if( family=='nonNegNormal' ){
      draws <- wAccProbNonNegNormal(
        X[[it]],c(n,p,T,1),Y,
        betaIn[it], betaOut[it], kappa, w[,it-1], w[,it], g2[it-1]
        ) }

    if( family=='gaussian' ){
      draws <- wAccProbGaussian(
        X[[it]],c(n,p,T,1),Y,
        betaIn[it], betaOut[it], kappa, w[,it-1], w[,it], g2[it-1]
        ) }

    if( family=='poisson' ){
      draws <- wAccProbPoisson(
        X[[it]],c(n,p,T,1),Y,
        betaIn[it], betaOut[it], kappa, w[,it-1], w[,it]
        ) }

    w[,it] <- draws[[1]] ; accRate[3] <- accRate[3] + draws[[2]] ; rm(draws)

    if( family=='nonNegNormal' ){
      g2[it] <- MCMCpack::rinvgamma(1,shape=shapeG2,scale=scaleG2)
      draws <- gammaAccProb(
        X[[it]],c(n,p,T,1),Y, 
        betaIn[it],betaOut[it],shapeG2,scaleG2, 
        w[,it],g2[it-1], g2[it]
        )
      g2[it] <- draws[[1]] ; accRate[4] <- accRate[4] + draws[[2]] ; rm(draws)
    }

    if( family=='gaussian' ){
      g2[it] <- MCMCpack::rinvgamma(1,shape=shapeG2,scale=scaleG2)
      draws <- gammaAccProbGaussian(
        X[[it]],c(n,p,T,1),Y, 
        betaIn[it],betaOut[it],shapeG2,scaleG2, 
        w[,it],g2[it-1], g2[it]
        )
      g2[it] <- draws[[1]] ; accRate[4] <- accRate[4] + draws[[2]] ; rm(draws)
    }

    # Step 4
    if(missData){
      for(t in 1:T){
        if(family == 'binomial'){
        Y <- imputeMissingBinomial(
          X[[it]], c(n,p,T), MM=missing[[t]]-1, Y, Ttt=t,
          BIN=betaIn[it], BOUT=betaOut[it], ww=w[,it]
          )}
        if(family == 'poisson'){
        Y <- imputeMissingPoisson(
          X[[it]], c(n,p,T), MM=missing[[t]]-1, Y, Ttt=t,
          BIN=betaIn[it], BOUT=betaOut[it], ww=w[,it]
          )}
        if(family == 'gaussian'){
        Y <- imputeMissingGaussian(
          X[[it]], c(n,p,T), MM=missing[[t]]-1, Y, Ttt=t,
          BIN=betaIn[it], BOUT=betaOut[it], ww=w[,it], g2 = g2[it]
          )}
      }
    }

    keeps = seq(burnin + 1, N, odens)
    # save results
    if(it > burnin){
      if(saveResults){
      if( it %in% round( quantile( (burnin+1):N, probs=seq(0,1,savePoints)) ) ){
        result <- list( Y=ifelse(missData, Y[keeps], Y), X=X[keeps], p=p, betaIn=betaIn[keeps], betaOut=betaOut[keeps], t2=t2[keeps], s2=s2[keeps], g2=g2[keeps],
          shapeT2=shapeT2, shapeS2=shapeS2, scaleT2=scaleT2, scaleS2=scaleS2,
          shapeG2=shapeG2, scaleG2=scaleG2, nuIn=nuIn, nuOut=nuOut,
          xiIn=xiIn, xiOut=xiOut, w=w[,keeps], accRate=accRate )
        save( result , file=fileName ) ; if(it!=N){ rm(result) }
      }
      }      
    } 

  if(progressBar){ setTxtProgressBar(pb,it) }
  } # end mcmc
  close(pb)
  })

  # output
        result <- list( Y=ifelse(missData, Y[keeps], Y), X=X[keeps], p=p, betaIn=betaIn[keeps], betaOut=betaOut[keeps], t2=t2[keeps], s2=s2[keeps], g2=g2[keeps],
          shapeT2=shapeT2, shapeS2=shapeS2, scaleT2=scaleT2, scaleS2=scaleS2,
          shapeG2=shapeG2, scaleG2=scaleG2, nuIn=nuIn, nuOut=nuOut,
          xiIn=xiIn, xiOut=xiOut, w=w[,keeps], accRate=accRate )
  return( result )

} # end function
