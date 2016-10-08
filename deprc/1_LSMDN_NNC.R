#Number of MCMC iterations
N=1000000
#Dimension of the Euclidean latent space
p=2
#Use log likelihood approximation (BOOLEAN)?
#If TRUE, how large a subsample n0?
llApprox = FALSE
if(llApprox) n0 = 100
#Are there missing edges?
MissData = FALSE
#If TRUE, construct Missing: Missing[[tt]]= c(1,2,3) => at time tt we have no row data on actors 1,2&3 
if(MissData){
  Missing <- list() #Enter as list manually, or if NAs are in data run the following:
  temp = which(is.na(Y),arr.ind=TRUE)
  for(tt in 1:dim(Y)[3]){
    Missing[[tt]] = unique(temp[which(temp[,3]==tt),1])
  }
  Y[temp] = 0;rm(temp)
}

#MCMC tuning parameters
tuneX <-   0.0075
tuneBIO <- 0.1
Kappa <-   175000
burnin = round(N/10)

# Load pkgs and functions
lsmdnPkgs = c('igraph', 'MCMCpack', 'inline',
  'RcppArmadillo','vegan')
loadPkg(lsmdnPkgs)
source(paste0(rFuncs, "functions_NNC.R"))
source(paste0(rFuncs, "initialize_NNC.R"))

###
###Run MCMC
###
pb <- txtProgressBar(min=2,max=N,style=3)
system.time({
  set.seed(1)
  for(it in 2:N){
      
    RN <- rnorm(n*TT*p)
    RNBIO <- rnorm(2)
    if(llApprox){
      if(it%%100==0){
        SUBSEQ = matrix(0,n,n0)
        for(i in 1:n){
          nOnes <- round(length(edgeList[[i]])/n*n0) #stratified sampling
          if(length(edgeList[[i]])>0){ nOnes <- max(nOnes,1) }
          SUBSEQ[i,1:nOnes] <- sample(edgeList[[i]],size=nOnes,replace=FALSE)
          SUBSEQ[i,(nOnes+1):n0] <- sample(c(1:n)[-c(i,edgeList[[i]])],size=n0-nOnes,replace=FALSE)
        }
      }
    }
    if(llApprox){
      Draws <- c.update2(X[[it-1]],c(n,p,TT,dinmax,doutmax),tuneX,Y,
                         Bin[it-1],Bout[it-1],tuneBIO,w[,it-1],
                         t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
                         nuOut,CAUCHY=0,RN,RNBIO,ELOUT,ELIN,SUBSEQ,DEGREE)
    }else{
      Draws <- c.update1(X[[it-1]],c(n,p,TT,1),tuneX,Y,
                         Bin[it-1],Bout[it-1],tuneBIO,w[,it-1],
                         t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
                         nuOut,CAUCHY=0,RN,RNBIO, g2[it-1])
    }
    X[[it]] <- Draws[[1]]
    Bin[it] <- Draws[[2]]
    Bout[it] <- Draws[[3]]
    AccRate <- AccRate + Draws[[4]]
    
    if(it==burnin){
      Xit0 <- t(X[[it]][,1,])
      for(tt in 2:TT) Xit0 <- rbind(Xit0,t(X[[it]][,tt,]))
    }
    if(it>burnin){
      XitCentered <- t(X[[it]][,1,])
      for(tt in 2:TT) XitCentered <- rbind(XitCentered,t(X[[it]][,tt,]))
      procr <- vegan::procrustes(X=Xit0,Y=XitCentered,scale=FALSE)$Yrot
      for(tt in 1:TT){
        X[[it]][,tt,] <- t(procr[((tt-1)*n+1):(tt*n),])
      }
    }
    if(it < N) X[[it+1]] <- X[[it]]
    
    #------------------Step 2--------------------------------
    Draws1 <- c.t2s2Parms(X[[it]],c(n,p,TT,1),shapeT2,
                          shapeS2,scaleT2,scaleS2)
    t2[it] <- rinvgamma(1,shape=Draws1[[1]],scale=Draws1[[2]])
    s2[it] <- rinvgamma(1,shape=Draws1[[3]],scale=Draws1[[4]])
    
    #------------------Step 3--------------------------------
    
    w[,it] <- rdirichlet(1,alpha=Kappa*w[,it-1])
    if(llApprox){
      Draws2 <- c.WAccProb2(X[[it]],c(n,p,TT,dinmax,doutmax),Y,
                            Bin[it],Bout[it],Kappa,w[,it-1],w[,it],
                            ELOUT,ELIN,SUBSEQ,DEGREE)
    }else{
      Draws2 <- c.WAccProb1(X[[it]],c(n,p,TT,1),Y,
                            Bin[it],Bout[it],Kappa,w[,it-1],w[,it], g2[it - 1])
    }
    w[,it] <- Draws2[[1]]
    AccRate[3] <- AccRate[3] + Draws2[[2]]
    
    #------------------Step 3.5--------------------------------
    g2[it] = rinvgamma(1,shape=shapeG2,scale=scaleG2)
    Draws3 = c.GammaAccProb1(X[[it]],c(n,p,TT,1),Y,
            Bin[it],Bout[it],shapeG2,scaleG2,
            w[,it-1],g2[it - 1], g2[it])
    g2[it] = Draws3[[1]]
    AccRate[4] = AccRate[4] + Draws3[[2]]
    #------------------Step 4--------------------------------
    
    if(MissData){
      for(tt in 1:TT){
        Y <- c.missing(X[[it]],c(n,p,TT),MMM=Missing[[tt]]-1,Y,Ttt=tt,
                       BETAIN=Bin[it],BETAOUT=Bout[it],WW=w[,it])
      }
    }
    
    if(it%%1000==0) save.image(outFileName)
    
    setTxtProgressBar(pb,it)
  }
  close(pb)
})
AccRate[1:4]/(it-1)
summary(AccRate[-c(1:3)]/(it-1))