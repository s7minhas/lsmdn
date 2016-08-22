####
if(Sys.info()['user']=='janus829' | Sys.info()['user']=='s7m'){
  source('~/Research/NetworkEvolution/Code/setup.R') }
####  

####  
# Load lsdm results
pick=list.files(outPath)[grepl('.rda',list.files(outPath))]
load( paste0(outPath, 'LSMDN_matlCoop2_quarterly.rda' ) ) 
####  

####  
# ran 268000 iterations before crashing
N = 268000
# eval full chain
par(mfrow=c(1,2))
pdf(file=paste0(graphicsPath, 'matlCoop2_quarterly_fullChain_conv.pdf'))
plot(1:N, Bin[1:N])
plot(1:N, Bout[1:N])
dev.off()
# with default burnin
pdf(file=paste0(graphicsPath, 'matlCoop2_quarterly_defaultBurnin_conv.pdf'))
plot((burnin+1):N, Bin[(burnin+1):N])
plot((burnin+1):N, Bout[(burnin+1):N])
dev.off()
# mod burnin
burnin2 = 200000
pdf(file=paste0(graphicsPath, 'matlCoop2_quarterly_modBurnin_conv.pdf'))
plot((burnin2+1):N, Bin[(burnin2+1):N])
plot((burnin2+1):N, Bout[(burnin2+1):N])
dev.off()
par(mfrow=c(1,1))
# try with mod burnin
burnin = burnin2
####  

####  
# Determine radii
postRad = rowMeans(w[,(burnin+1):N])
priorRad = w[,1]
radii = data.frame(ccode=dimnames(Y)[[2]], postRad, priorRad, ncol=3, stringsAsFactors=FALSE)
radii$postRadRank = rank(radii$postRad)
radii$priorRadRank = rank(radii$priorRad)
radii$cname = panel$cname[match(radii$ccode, panel$ccode)]
radii$diff = radii$postRad - radii$priorRad
radii$rankDiff = radii$postRadRank - radii$priorRadRank
radii[order(abs(radii$diff), decreasing=TRUE),][1:10,]
####  

####  
# Load pkgs and functions
lsmdnPkgs = c('igraph', 'MCMCpack', 'inline',
  'RcppArmadillo','vegan', 'circular')
loadPkg(lsmdnPkgs)
source(paste0(rFuncs, "functions_NNC.R"))
####  

####  
# Organize lat space array
N1 <- N-burnin

EBin <- mean(Bin[(burnin+1):N]) # global parameter reflecting popularity (p. 1648, S&C '15)
EBout <- mean(Bout[(burnin+1):N]) # global parameter reflecting social activity (p. 1648)
Ew <- apply(w[,(burnin+1):N],1,mean) # radii for countries, represents social reach (p. 1648)
EX <- array(0,dim=c(p,TT,n))
for(it in c(1:N1)){ EX <- EX + X[[burnin+it]] } # Expected positiong of countries
EX <- EX/N1  
####  

####  
# Nodal Influence

# Parameters
lambda <- 0.05 # tuning parameter
p0 <- 0.5 # prior probability zero
kappaij <- matrix(0,n,n)
temp <- numeric(TT-1)
X1 <- X2 <- array(0,c(n,TT,N-burnin+1))

for(it in burnin:N){
  X1[,,it-burnin+1] <- t(X[[it]][1,,])
  X2[,,it-burnin+1] <- t(X[[it]][2,,])
}

# Calculate kappaij
## kappa is the ratio of there is an interaction or there is not an interaction
count=0
pb <- txtProgressBar(min=1,max=n*(n-1)/2,style=3)
for(i in 1:(n-1)){
  for( j in (i+1):n){ 
    kappaij[i,j] <- c.postzeroprob(X1[i,,],X2[i,,],X1[j,,],X2[j,,],s2[burnin:N],lambda,p0)
    kappaij[j,i] <- c.postzeroprob(X1[j,,],X2[j,,],X1[i,,],X2[i,,],s2[burnin:N],lambda,p0)
    count = count+1
    setTxtProgressBar(pb,count)
  }
}
close(pb)
save(kappaij, file=paste0(outPath, 'kappa_', 'LSMDN_matlCoop2_quarterly.rda'))
# load(paste0(outPath, 'kappa_', 'LSMDN_matlCoop2_quarterly.rda'))

# Eqn 11
Prob0 <- matrix(0,n,n)
for(i in 1:n){
  for(j in c(1:n)[-i]){
    Prob0[i,j] <- 1/(1+kappaij[i,j])
  }
}

# Here we are just calculating whether i and j are moving in the same direction
NegInfl <- matrix(0,n,n)
for( i in 1:n){
  for(j in c(1:n)[-i]){
    temp <- numeric(TT-1)
    for(tt in 2:TT){
      temp[tt-1] <- atan2(EX[2,tt,j]-EX[2,tt-1,i],EX[1,tt,j]-EX[1,tt-1,i])-
        atan2(EX[2,tt,i]-EX[2,tt-1,i],EX[1,tt,i]-EX[1,tt-1,i])
    }
    if(mean(abs(temp))>pi/2){ Prob0[i,j] <- 1;NegInfl[i,j] <- 1;}#
  }
}

# Basically, this block of code identifies the influencers and influencees
index <- which(apply(Prob0,1,function(x) sum(x < 0.5)-1)>0) # identifies the countries that are influencers in the network
Results <- apply(Prob0[index,],1,function(x) return(which(x < 0.5))) # find the targets of influence by row
for(i in 1:length(index)){
  Results[[i]] <- Results[[i]][which(Results[[i]] != index[i])]
}

# look across each influencer and their influencees
# look through time and check to see if the distance between the points is ever less than the
# radius of the influencer, and if it is then you add it in.
Refined <- list()
for(i in 1:length(index)){
  temp <- NULL
  for(j in 1:length(Results[[i]])){
    for(tt in 1:TT){
      dijt <- sqrt(t(EX[,tt,index[i]]-EX[,tt,Results[[i]][j]])%*%(EX[,tt,index[i]]-EX[,tt,Results[[i]][j]]))
      if(dijt < Ew[index[i]] | dijt < Ew[Results[[i]][j]]){
        temp <- c(temp,Results[[i]][j])
      }
      toAdd = unique(temp)
      if(is.null(toAdd)){ toAdd=NA } 
      Refined[[i]]<-toAdd
    }
  }
}
length(index) ; length(Refined)
indexOrig = index
RefinedOrig = Refined
index = index[-which( lapply(Refined, function(x){ unique(is.na(x)) }) == TRUE )]
Refined = Refined[-which( lapply(Refined, function(x){ unique(is.na(x)) }) == TRUE )]
length(index) ; length(Refined)
print(index)
print(Refined)

cntries = dimnames(Y)[[2]]
tpds = dimnames(Y)[[3]]
muAngArr = array(0, dim=c(n,n,2,TT-1), 
  dimnames=list(cntries,cntries,c('angle','mu'),tpds[2:length(tpds)]) )
for(ii in 1:length(index)){
  for(ind in 1:length(Refined[[ii]])){
    for(tt in 2:TT){
      angle <- -atan2(EX[2,tt,Refined[[ii]][ind]]-EX[2,tt-1,index[ii]],
                                EX[1,tt,Refined[[ii]][ind]]-EX[1,tt-1,index[ii]])
      Rotate <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),2,2)
      mu <- (Rotate%*%(EX[,tt,index[ii]]-EX[,tt-1,index[ii]]))[1]
      muAngArr[cntries[ii],cntries[Refined[[ii]][ind]],'angle',tpds[tt]] = angle
      muAngArr[cntries[ii],cntries[Refined[[ii]][ind]],'mu',tpds[tt]] = mu
    }
  }
}
save(muAngArr, file=paste0(outPath, 'muAngArr.rda'))

# Now calculate mu
conc <- list()
temp <- numeric(TT-1)
for( i in 1:length(index)){
  conc[[i]] <- matrix(0,TT-1,3)
  muX <- numeric(length(Refined[[i]]))
  Angle <- matrix(0,TT-1,length(muX))
  for(ind in 1:length(Refined[[i]])){
    for(tt in 2:TT){
      Angle[tt-1,ind] <- -atan2(EX[2,tt,Refined[[i]][ind]]-EX[2,tt-1,index[i]],
                                EX[1,tt,Refined[[i]][ind]]-EX[1,tt-1,index[i]])
      Rotate <- matrix(c(cos(Angle[tt-1,ind]),sin(Angle[tt-1,ind]),-sin(Angle[tt-1,ind]),cos(Angle[tt-1,ind])),2,2)
      temp[tt-1] <- (Rotate%*%(EX[,tt,index[i]]-EX[,tt-1,index[i]]))[1]
    }
    muX[ind] <- max(mean(temp),0)
  }
  ind <- which.max(muX)
  conc[[i]][,3] <- rep(Refined[[i]][ind],TT-1)
  muX <- max(muX)
  for(tt in 2:TT){
    temp[tt-1] <- sqrt(t(EX[,tt,index[i]]-EX[,tt-1,index[i]])%*%(EX[,tt,index[i]]-EX[,tt-1,index[i]]))
    conc[[i]][tt-1,1] <- muX*temp[tt-1]
  }
  conc[[i]][,1] <- conc[[i]][,1]/mean(s2[burnin:N])
  conc[[i]][,2] <- -Angle[,ind]
}

# Plot stuff with von mises
yl <- xl <- matrix(0,length(index),2)
for(i in 1:length(index)){
  tempX <- mean(range(EX[1,,c(Refined[[i]],index[i])]))
  tempY <- mean(range(EX[2,,c(Refined[[i]],index[i])]))
  rangeX <- range(EX[1,,c(Refined[[i]],index[i])])
  rangeY <- range(EX[2,,c(Refined[[i]],index[i])])
  Size <- max(diff(rangeX),diff(rangeY))
  xl[i,] <- tempX + c(-1,1)*Size/2
  yl[i,] <- tempY + c(-1,1)*Size/2
  xl[i,] <- xl[i,] + c(-1,1)*0.05*Size
  yl[i,] <- yl[i,] + c(-1,1)*0.05*Size
}

Theta <- seq(from=0,to=2*pi,length=2001)[-2001]
MAX <- NULL
for(i in 1:length(index)){
  for(tt in 1:(TT-1)){
    MAX <- c(MAX,dvonmises(Theta,mu=circular(conc[[i]][tt,2]),kappa=conc[[i]][tt,1]))
  }
}

MAX <- max(MAX)
RSc <- 0.05*apply(xl,1,diff)

for(i in 1:length(index)){
  cols <- "black"; 
  cols <- c(cols,rep("blue",length(Refined[[i]])));
  PCH <- 16;
  PCH <- c(PCH, rep(17,length(Refined[[i]])))
  
  temp <- NULL
  for(tt in 1:TT) temp <- rbind(temp,t(EX[,tt,c(index[i],Refined[[i]])]))
  jpeg(file=paste0(graphicsPath, "nodeInfluence/NodalInfluence_Actor",index[i],".jpg"))
  par(mar=c(1,1,1,1)+0.1)
  plot(temp,col="white",ylab="",xlab="",xlim=xl[i,],ylim=yl[i,],xaxt="n",yaxt="n")
  for(tt in 1:TT){ 
    points(t(EX[,tt,c(index[i],Refined[[i]])]),pch=PCH,col=cols,cex=2)
    if(tt>1) arrows(EX[1,tt-1,c(index[i],Refined[[i]])],EX[2,tt-1,c(index[i],Refined[[i]])],
                    EX[1,tt,c(index[i],Refined[[i]])],EX[2,tt,c(index[i],Refined[[i]])],
                    length=0.2,col=cols)
    if(tt < TT){
      dTheta <- dvonmises(Theta,mu=circular(conc[[i]][tt,2]),kappa=conc[[i]][tt,1])
      ind <- 1000*(MAX-dTheta)/MAX; 
      if(min(ind)<1) ind[which.min(ind)]=1; 
      if(max(ind)>1000) ind[which.max(ind)]=1000
      #MAX=max(dTheta)
      #Cols <- gray(1:1000/1000)[1000*(max(dTheta)-dTheta)/max(dTheta)]
      Cols <- gray(1:1000/1000)[ind]
      tempSeq <- seq(from=1.75,to=0.000001,length=1000)
      #CEX <- tempSeq[1000*(max(dTheta)-dTheta)/max(dTheta)]
      CEX <- tempSeq[ind]
      par(new=TRUE)
      plot(x=EX[1,tt,index[i]]+RSc[i]*cos(Theta),y=EX[2,tt,index[i]]+RSc[i]*sin(Theta),
           col=Cols,type="b",pch=16,ylab="",xlab="",cex=CEX,
           ylim=yl[i,],xlim=xl[i,],xaxt="n",yaxt="n",bty="n")
    }
  }
  iName = panel$cname[match(dimnames(Y)[[1]][index[i]], panel$ccode)]
  jIndex = unique(conc[[i]][,3])
  jName = panel$cname[match(dimnames(Y)[[1]][jIndex], panel$ccode)]
  text(x=mean(xl[i,]), y=yl[i,][2], iName, col='black')
  text(x=mean(xl[i,]), y=yl[i,][1], jName, col='blue')
dev.off()
}

#---

Theta <- seq(from=-pi,to=pi,length=500)
conc <- list()
temp <- numeric(TT-1)
for( i in 1:length(index)){
  conc[[i]] <- numeric(length(Refined[[i]]))
  for(j in 1:length(conc[[i]])){
    ind <- Refined[[i]][j]
    for(tt in 2:TT){
      Angle <- -atan2(EX[2,tt,ind]-EX[2,tt-1,index[i]],EX[1,tt,ind]-EX[1,tt-1,index[i]])
      Rotate <- matrix(c(cos(Angle),sin(Angle),-sin(Angle),cos(Angle)),2,2)
      temp[tt-1] <- (Rotate%*%(EX[,tt,index[i]]-EX[,tt-1,index[i]]))[1]
    }
    muX <- max(mean(temp),0)
    for(tt in 2:TT){
      temp[tt-1] <- sqrt((EX[,tt,index[i]]-EX[,tt-1,index[i]])%*%(EX[,tt,index[i]]-EX[,tt-1,index[i]]))
    }
    conc[[i]][j] <- muX*mean(temp)/mean(s2[burnin:N])
  }
}

for(i in 1:length(index)){
  MAX <- dvonmises(0,mu=circular(0),kappa=max(conc[[i]]))
  
  jpeg(filename=paste0(graphicsPath, "nodeInfluenceCircle/NodInflAggregateWrapped_Actor",index[i],".jpg"))
  par(mar=c(2, 2, 2, 2))
  for(j in 1:length(Refined[[i]])){
    dTheta <- dvonmises(Theta,mu=circular(0),kappa=conc[[i]][order(conc[[i]])][j])
    ind <- 1000*(MAX-dTheta)/MAX 
    if(min(ind)<1) ind[which.min(ind)]=1; 
    if(max(ind)>1000) ind[which.max(ind)]=1000
    tempSeq <- seq(from=3,to=0.0001,length=1000)
    CEX <- tempSeq[ind]
    #     Cols <- gray(1:2000/2000)[2000*(max(dTheta)-dTheta)/max(dTheta)]
    Cols <- gray(1:1000/1000)[ind]
    plot(j/length(Refined[[i]])*cos(Theta),j/length(Refined[[i]])*sin(Theta),
         col=Cols,type="b",pch=16,ylab="",xlab="",cex=CEX,
         ylim=c(-1.2,1.2),xlim=c(-1.2,1.2),xaxt="n",yaxt="n")
    if(j < length(Refined[[i]])) par(new=TRUE)    
  }
  text(x=c(1.2,1.1*cos(pi/4),1.1*cos(3*pi/4),1.1*cos(pi/4),1.1*cos(3*pi/4)),
       y=c(0,1.3*sin(pi/4),1.3*sin(3*pi/4),1.3*sin(-pi/4),1.3*sin(-3*pi/4)),
       labels=c(0,expression(frac(pi,4)),expression(frac(3*pi,4)),
                expression(-frac(pi,4)),expression(-frac(3*pi,4))))
  iName = panel$cname[match(dimnames(Y)[[1]][index[i]], panel$ccode)]
  jName = panel$cname[match(dimnames(Y)[[1]][Refined[[i]][j]], panel$ccode)]  
  text(x=sum(-1.2,1.2)/2, y=1.2, iName, col='black')
  text(x=sum(-1.2,1.2)/2, y=-1.2, jName, col='blue')
  dev.off()
  jpeg(filename=paste0(graphicsPath, "nodeInfluenceDensity/NodInflAggregate_Actor",index[i],".jpg"))
  par(mar=c(4, 4, 2, 2))
  for(j in 1:length(Refined[[i]])){
    curve(dvonmises(x,mu=0,kappa=conc[[i]][order(conc[[i]],decreasing=TRUE)][j]),
          xlim=c(-pi,pi),col=1,
          xlab="",ylab="",xaxt="n",
          ylim=c(0.9*dvonmises(pi,mu=0,kappa=range(conc[[i]])[1]),
                 1.1*dvonmises(0,mu=0,kappa=range(conc[[i]])[2])))
    abline(h=0,v=0,lty=3)
    iName = panel$cname[match(dimnames(Y)[[1]][index[i]], panel$ccode)]
    jName = panel$cname[match(dimnames(Y)[[1]][Refined[[i]][j]], panel$ccode)]  
    text(x=sum(-pi,pi)/2, y=1.1*dvonmises(0,mu=0,kappa=range(conc[[i]])[2]), iName, col='black')
    text(x=sum(-pi,pi)/2, y=0.9*dvonmises(pi,mu=0,kappa=range(conc[[i]])[1]), jName, col='blue')    
    if(j < length(Refined[[i]])) par(new=TRUE)
  }
  axis(1,at=c(-pi,-pi/2,0,pi/2,pi),labels=c(expression(-pi),expression(-pi/2),
                                            0,expression(pi/2),expression(pi)))
  dev.off()
}