###
###Initial \beta_{IN} and \beta_{OUT}; scale latent positions
###
initialize.wrap <- function(x){
  -c.initialize1(X[[1]],c(n,p,TT),Y,XSCALE=1/n,
                 BETAIN=x[1],BETAOUT=x[2],w[,1])
}
initialize.grad.wrap <- function(x){
  -c.initialize1.grad(X[[1]],c(n,p,TT),Y,XSCALE=1/n,
                      BETAIN=x[1],BETAOUT=x[2],w[,1])[2:3]
}
Optim <- optim(par=c(1,1),fn=initialize.wrap,
               gr=initialize.grad.wrap,method="BFGS")
X[[1]] <- X[[1]]/n
Bin[1] <- max(Optim$par[1],1e-4)
Bout[1] <- max(Optim$par[2],1e-4)

Xit0 <- t(X[[1]][,1,])
for(tt in 2:TT)Xit0 <- rbind(Xit0,t(X[[1]][,tt,]))
Xit0 <- Xit0 -
  kronecker(rep(1,n*TT),matrix(apply(Xit0,2,mean),1,p))