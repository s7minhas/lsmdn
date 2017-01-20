lsmdnSv <- function(
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
  keep = length(seq(burnin + 1, N, odens)) + 1
  w = matrix(0,keep,N)
  X = list()
  t2 = numeric(keep)
  s2 = numeric(keep)
      if( family=='nonNegNormal' | family == "gaussian" ){
  g2 = numeric(keep)}
  betaIn = numeric(keep)
  betaOut = numeric(keep)
  # get init values if no fitted values provided
  if( is.null( startVals ) ){

    tmp <- getStartingValues(
      Y=Y, p=p, family=family, 
      llApprox=llApprox, missData=missData, 
      N=N, seed=seed
      )

    # unpack
    Y<-tmp$Y ; w0<-tmp$w[,1] ; X0 <-tmp$X ; betaIn0<-tmp$betaIn ; betaOut0<-tmp$betaOut ; 
    nuIn<-tmp$nuIn ; nuOut<-tmp$nuOut ; xiIn<-tmp$xiIn ; xiOut<-tmp$xiOut 
    t20<-tmp$t2 ; shapeT2<-tmp$shapeT2 ; scaleT2<-tmp$scaleT2
    s20<-tmp$s2 ; shapeS2<-tmp$shapeS2 ; scaleS2<-tmp$scaleS2
     accRate<-tmp$accRate ; g20 <- NULL ; shapeG2 <- NULL ; scaleG2 <- NULL

    if( family=='nonNegNormal' | family == "gaussian" ){
      g20<-tmp$g2 ; shapeG2<-tmp$shapeG2 ; scaleG2<-tmp$scaleG2 }

    if( llApprox & family=='binomial' ){
      dInMax<-tmp$dInMax ; dOutMax<-tmp$dOutMax ; elOut<-tmp$elOut ; elIn<-tmp$elIn
      degree<-tmp$degree ; edgeList<-tmp$edgeList; n0<-tmp$n0 }

    rm(tmp) # cleanup
  }

  if( !is.null( startVals ) ){
    tmp = startVals
    nIter = length(tmp$X)
    Y<-tmp$Y ; w0<-tmp$w[,nIter] ; X0 <-tmp$X[[nIter]] ; betaIn0<-tmp$betaIn[nIter] ; betaOut0<-tmp$betaOut[nIter] ; 
    nuIn<-tmp$nuIn ; nuOut<-tmp$nuOut ; xiIn<-tmp$xiIn ; xiOut<-tmp$xiOut 
    t20<-tmp$t2[nIter] ; shapeT2<-tmp$shapeT2 ; scaleT2<-tmp$scaleT2
    s20<-tmp$s2[nIter] ; shapeS2<-tmp$shapeS2 ; scaleS2<-tmp$scaleS2
    n0<-tmp$n0 ; accRate<-tmp$accRate ; g20 <- NULL ; shapeG2 <- NULL ; scaleG2 <- NULL

    if( family=='nonNegNormal' | family == "gaussian" ){
      g20<-tmp$g2[nIter] ; shapeG2<-tmp$shapeG2 ; scaleG2<-tmp$scaleG2 }

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
        X0,c(n,p,T,dInMax,dOutMax),tuneX,Y, 
        betaIn0,betaOut0,tuneBIO,w0,
        t20,s20,xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO,elOut,elIn,subseq,degree
        ) }

    if( !llApprox & family=='binomial' ){
      draws <- updateBinom(
        X0,c(n,p,T,1),tuneX,Y, 
        betaIn0,betaOut0,tuneBIO,w0,
        t20,s20,xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO
        ) }

    if( family=='nonNegNormal' ){
      draws <- updateNonNegNorm(
        X0,c(n,p,T,1),tuneX,Y, 
        betaIn0,betaOut0,tuneBIO,w0,
        t20,s20,g20,xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO
        ) }

    if( family=='gaussian' ){
      draws <- updateGaussian(
        X0,c(n,p,T,1),tuneX,Y, 
        betaIn0,betaOut0,tuneBIO,w0,
        t20,s20,g20,xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO
        ) }

    if( family=='poisson' ){
      draws <- updatePoisson(
        X0,c(n,p,T,1),tuneX,Y, 
        betaIn0,betaOut0,tuneBIO,w0,
        t20,s20,xiIn,xiOut,nuIn,
        nuOut,Cauchy=0,RN,RNBIO
        ) }

    #
    X0 <- draws[[1]] ; betaIn0 <- draws[[2]] ; betaOut0 <- draws[[3]]
    accRate <- accRate + draws[[4]] 
    if((it - burnin) %% odens == 0){
       ind = (it - burnin)/odens + 1
       X[[ind]] <- X0 ; betaIn[ind] <- betaIn0 ; betaOut[ind] <- betaOut0
       }
    rm(draws)
    #
    if(it==burnin){
      xIter0 <- t(X[[1]][,1,])
      for(t in 2:T) xIter0 <- rbind(xIter0,t(X[[1]][,t,])) }

    if(it>burnin & ((it - burnin) %% odens == 0)){
      ind = (it - burnin)/odens + 1
      xIterCentered <- t(X[[ind]][,1,])
      for(t in 2:T){ xIterCentered <- rbind(xIterCentered,t(X[[ind]][,t,])) }
      procr <- vegan::procrustes(X=xIter0,Y=xIterCentered,scale=FALSE)$Yrot
      for(t in 1:T){ X[[ind]][,t,] <- t(procr[((t-1)*n+1):(t*n),]) } }

    if(it< (length(seq(burnin + 1, N, odens)) * odens + burnin) & it>burnin & ((it - burnin) %% odens == 0)){ 
      ind = (it - burnin)/odens + 1
      X[[ind+1]] <- X[[ind]] }

    # Step 2
    draws <- t2s2Parms(X[[it]], c(n,p,T,1), shapeT2, shapeS2, scaleT2, scaleS2)
    t20 <- MCMCpack::rinvgamma(1,shape=draws[[1]],scale=draws[[2]])
    s20 <- MCMCpack::rinvgamma(1,shape=draws[[3]],scale=draws[[4]]) ; rm(draws)
    if((it - burnin) %% odens == 0){
       ind = (it - burnin)/odens + 1
       t2[ind] <- t20 ; s2[ind] <- s20
       }


    # Step 3
    wProp <- MCMCpack::rdirichlet(1,alpha=kappa*w0)
    if(llApprox & family=='binomial'){
      draws <- wAccProb_llApprox(
        X0,c(n,p,T,dInMax,dOutMax),Y,
        betaIn0, betaOut0, kappa, w0, wProp,
        elOut, elIn, subseq, degree
        ) }

    if( !llApprox & family=='binomial' ){
      draws <- wAccProb(
        X0,c(n,p,T,1),Y,
        betaIn0, betaOut0, kappa, w0, wProp
        ) }

    if( family=='nonNegNormal' ){
      draws <- wAccProbNonNegNormal(
        X0,c(n,p,T,1),Y,
        betaIn0, betaOut0, kappa, w0, wProp, g20
        ) }

    if( family=='gaussian' ){
      draws <- wAccProbGaussian(
        X0,c(n,p,T,1),Y,
        betaIn0, betaOut0, kappa, w0, wProp, g20
        ) }

    if( family=='poisson' ){
      draws <- wAccProbPoisson(
        X0,c(n,p,T,1),Y,
        betaIn0, betaOut0, kappa, w0, wProp
        ) }

    w0 <- draws[[1]] ; accRate[3] <- accRate[3] + draws[[2]] ; rm(draws)
    if((it - burnin) %% odens == 0){
       ind = (it - burnin)/odens + 1
       w[,ind] <- w0
       }

    if( family=='nonNegNormal' ){
      g2Prop <- MCMCpack::rinvgamma(1,shape=shapeG2,scale=scaleG2)
      draws <- gammaAccProb(
        X0,c(n,p,T,1),Y, 
        betaIn0,betaOut0,shapeG2,scaleG2, 
        w0,g20, g2Prop
        )
      g20 <- draws[[1]] ; accRate[4] <- accRate[4] + draws[[2]] ; rm(draws)
        if((it - burnin) %% odens == 0){
       ind = (it - burnin)/odens + 1
       g2[ind] <- g20
       }

    }

    if( family=='gaussian' ){
      g2Prop <- MCMCpack::rinvgamma(1,shape=shapeG2,scale=scaleG2)
      draws <- gammaAccProbGaussian(
        X0,c(n,p,T,1),Y, 
        betaIn0,betaOut0,shapeG2,scaleG2, 
        w0,g20, g2Prop
        )
      g20 <- draws[[1]] ; accRate[4] <- accRate[4] + draws[[2]] ; rm(draws)
        if((it - burnin) %% odens == 0){
       ind = (it - burnin)/odens + 1
       g2[ind] <- g20
       }
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

    # save results
    if(it > burnin){
      if(saveResults){
      if( it %in% round( quantile( (burnin+1):N, probs=seq(0,1,savePoints)) ) ){
        result <- list( Y=Y, X=X, p=p, betaIn=betaIn, betaOut=betaOut, t2=t2, s2=s2, g2=g2,
          shapeT2=shapeT2, shapeS2=shapeS2, scaleT2=scaleT2, scaleS2=scaleS2,
          shapeG2=shapeG2, scaleG2=scaleG2, nuIn=nuIn, nuOut=nuOut,
          xiIn=xiIn, xiOut=xiOut, w=w, accRate=accRate )
        save( result , file=fileName ) ; if(it!=N){ rm(result) }
      }
      }      
    } 

  if(progressBar){ setTxtProgressBar(pb,it) }
  } # end mcmc
  close(pb)
  })

  # output
        result <- list( Y= Y, X=X, p=p, betaIn=betaIn, betaOut=betaOut, t2=t2, s2=s2, g2=g2,
          shapeT2=shapeT2, shapeS2=shapeS2, scaleT2=scaleT2, scaleS2=scaleS2,
          shapeG2=shapeG2, scaleG2=scaleG2, nuIn=nuIn, nuOut=nuOut,
          xiIn=xiIn, xiOut=xiOut, w=w, accRate=accRate )
  return( result )

} # end function
