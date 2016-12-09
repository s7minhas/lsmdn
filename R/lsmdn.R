#' LSMDN model fitting routine for longitudinal network data
#' 
#' @param Y an n x n x T array of relational matrices, 
#' where the third dimension corresponds to different time periods.
#' @param p number of latent dimensions
#' @param family type of model to run. Options include 'normal', 'nonNegNormal', 'poisson', 'binomial'. 
#' @param llApprox logical indicating whether or not to utilize log-likelihood 
#' approximationg. Only available for binomial model types.
#' @param missData logical indicating whether to impute missing data.
#' @param tuneX tuneX
#' @param tuneBIO tuneBIO
#' @param kappa kappa
#' @param s2Init starting value for s2
#' @param t2Init starting value for t2
#' @param xLatPos starting actor positions in latent space
#' @param betaInInit starting value for betaIn
#' @param betaOutInit starting value for betaOut
#' @param nuIn starting value for nuIn
#' @param nuOut starting value for nuOut
#' @param xiIn starting value for xiIn
#' @param xiOut starting value for xiOut
#' @param shapeT2 shape parameter for t2
#' @param scaleT2 scale parameter for t2
#' @param shapeS2 shape parameter for s2
#' @param scaleS2 shape parameter for s2
#' @param burnin burnin
#' @param N number of MCMC iterations.
#' @param seed random seed
#' @param startVals Fitted result from previous model run.
#' @param savePoints chain intervals to save at 
#' @param fileName "lsmdn.rda" or wtv you want, make sure to specify a path as well
#' @usage lsmdn( Y, p=2, family='binomial', llApprox=FALSE, missData=FALSE, N=1000, seed=6886) 
#' @return returns list of starting values:
#' \item{w}{weights}
#' \item{X}{initial actor latent space positions calculated via GMDS}
#' \item{betaIn}{add desc}
#' \item{betaOut}{add desc}
#' \item{nuIn}{add desc}
#' \item{nuOut}{add desc}
#' \item{xiIn}{add desc}
#' \item{xiOut}{add desc}
#' \item{t2}{add desc}
#' \item{shapeT2}{add desc}
#' \item{scaleT2}{add desc}
#' \item{s2}{add desc}
#' \item{shapeS2}{add desc}
#' \item{scaleS2}{add desc}
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
  Y, p=2, family='binomial', llApprox=FALSE, missData=FALSE, 
  tuneX=0.0075, tuneBIO=0.1, kappa=175000, burnin=round(N/10),
  s2Init=NULL, t2Init=NULL, g2Init=NULL, xLatPos=NULL, betaInInit=NULL, betaOutInit=NULL,
  nuIn=NULL, nuOut=NULL, xiIn=NULL, xiOut=NULL, shapeT2=NULL, scaleT2=NULL, 
  shapeS2=NULL, scaleS2=NULL, shapeG2=NULL, scaleG2=NULL, 
  N=1000, seed=6886, startVals=NULL, 
  savePoints=.10, fileName='lsmdnModel.rda'){

  # add in some warnings
  ## llApprox only works for binomial family
  ## 

  #
  set.seed(seed)  
  n <- dim(Y)[1]
  T <- dim(Y)[3]    

  # get init values if no fitted values provided
  if( is.null( startVals ) ){

    tmp <- getStartingValues(Y, p, family, llApprox, missData, N, seed)

    # unpack
    Y<-tmp$Y ; w<-tmp$w ; X <-tmp$X ; betaIn<-tmp$betaIn ; betaOut<-tmp$betaOut ; 
    nuIn<-tmp$nuIn ; nuOut<-tmp$nuOut ; xiIn<-tmp$xiIn ; xiOut<-tmp$xiOut 
    t2<-tmp$t2 ; shapeT2<-tmp$shapeT2 ; scaleT2<-tmp$scaleT2
    s2<-tmp$s2 ; shapeS2<-tmp$shapeS2 ; scaleS2<-tmp$scaleS2
    n0<-tmp$n0 ; accRate<-tmp$accRate

    if( family=='nonNegNormal' ){
      g2<-tmp$g2 ; shapeG2<-tmp$shapeG2 ; scaleG2<-tmp$scaleG2 }

    if( llApprox & family=='binomial' ){
      dInMax<-tmp$dInMax ; dOutMax<-tmp$dOutMax ; elOut<-tmp$elOut ; elIn<-tmp$elIn
      degree<-tmp$degree ; edgeList<-tmp$edgeList }

    rm(tmp) # cleanup
  }

  if( !is.null( startVals ) ){
    print('need to add code for unpacking')
  }

  # start mcmc
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
          set.seed(seed) ; subseq[i,1:nOnes] <- sample(edgeList[[i]],size=nOnes,replace=TRUE) # should replace be false?      
          set.seed(seed) ; subseq[i,(nOnes+1):n0] <- sample(c(1:n)[-c(i,edgeList[[i]])],size=n0-nOnes,replace=TRUE)
        } }

      draws <- updateBinomLogLikeApprox(
        X[[it-1]],c(n,p,T,dInMax,dOutMax),tuneX,Y, 
        betaIn[it-1],betaOut[it-1],tuneBIO,w[,it-1],
        t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO,elOut,elIn,subseq,degree ) }

    if( !llApprox & family=='binomial' ){
      draws <- updateBinom(X[[it-1]],c(n,p,T,1),tuneX,Y, 
        betaIn[it-1],betaOut[it-1],tuneBIO,w[,it-1],
        t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO) }

    if( family=='nonNegNormal' ){
      draws <- updateNonNegNorm(X[[it-1]],c(n,p,T,1),tuneX,Y, 
        betaIn[it-1],betaOut[it-1],tuneBIO,w[,it-1],
        t2[it-1],s2[it-1],g2[it-1],xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO) }

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
    draws2 <- t2s2Parms(X[[it]], c(n,p,T,1), shapeT2, shapeS2, scaleT2, scaleS2)
    shapeT2<-draws2[[1]] ; scaleT2<-draws2[[2]]
    shapeS2<-draws2[[3]] ; scaleS2<-draws2[[4]] ; rm(draws2)
    set.seed(seed) ; t2[it] <- MCMCpack::rinvgamma(1,shape=shapeT2,scale=scaleT2)
    set.seed(seed) ; s2[it] <- MCMCpack::rinvgamma(1,shape=shapeS2,scale=scaleS2)  

    # Step 3
    set.seed(seed) ; w[,it] <- MCMCpack::rdirichlet(1,alpha=kappa*w[,it-1])
    if(llApprox & family=='binomial'){
      draws3 <- wAccProb_llApprox(X[[it]],c(n,p,T,dInMax,dOutMax),Y,
        betaIn[it], betaOut[it], kappa, w[,it-1], w[it],
        elOut, elIn, subseq, degree ) }

    if( !llApprox & family=='binomial' ){
      draws3 <- wAccProb(X[[it]],c(n,p,T),Y,
        betaIn[it], betaOut[it], kappa, w[,it-1], w[,it]) }

    if( family=='nonNegNormal' ){
      draws3 <- wAccProb_nnn(X[[it]],c(n,p,T),Y,
        betaIn[it], betaOut[it], kappa, w[,it-1], w[,it], g2[it-1]) }

    w[,it] <- draws3[[1]] ; accRate[3] <- accRate[3] + draws3[[2]] ; rm(draws3)

    if( family=='nonNegNormal' ){
      g2[it] <- MCMCpack::rinvgamma(1,shape=shapeG2,scale=scaleG2)
      draws3 <- gammaAccProb(X[[it]],c(n,p,T,1),Y, 
        betaIn[it],betaOut[it],shapeG2,scaleG2, 
        w[,it-1],g2[it-1], g2[it])
      g2[it] <- draws3[[1]] ; accRate[4] <- accRate[4] + draws3[[2]] ; rm(draws3)
    }

    # Step 4
    if(missData){
      for(t in 1:T){
        Y <- imputeMissingNet(X[[it]], c(n,p,T), MM=missing[[t]]-1, Y, Ttt=t,
          BIN=betaIn[it], BOUT=betaOut[it], ww=w[,it])
      }
    }

    # save results
    if(it > burnin){
      if( it %in% round(quantile(1:(N-burnin), probs=seq(0,1,savePoints))[-1] ) ){
        result <- list( Y=Y, X=X, p=p, betaIn=betaIn, betaOut=betaOut, t2=t2, s2=s2, g2=g2,
          shapeT2=shapeT2, shapeS2=shapeS2, scaleT2=scaleT2, scaleS2=scaleS2,
          shapeG2=shapeG2, scaleG2=scaleG2, nuIn=nuIn, nuOut=nuOut,
          xiIn=xiIn, xiOut=xiOut, w=w, accRate=accRate )
        save( result , file=fileName ) ; rm(result)
      }      
    } 

  } # end mcmc

  # output
  result <- list(Y=Y, X=X, p=p, betaIn=betaIn, betaOut=betaOut, t2=t2, s2=s2, g2=g2, 
    shapeT2=shapeT2, shapeS2=shapeS2, scaleT2=scaleT2, scaleS2=scaleS2,
    shapeG2=shapeG2, scaleG2=scaleG2, nuIn=nuIn, nuOut=nuOut,
    xiIn=xiIn, xiOut=xiOut, w=w, accRate=accRate )
  return( result )

} # end function
