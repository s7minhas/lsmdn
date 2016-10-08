####  
# Organize lat space array
EX = array(0,dim=c(p,TT,n))
for(it in (burnin+1):N){ EX = EX + X[[it]] }
EX = EX/(N-burnin)
dimnames(EX) = list(NULL, dimnames(Y)[[3]], dimnames(Y)[[1]])

if(quarterly){
  save(EX, file=paste0(outPath,'EX_Qly.rda')) } else {
  save(EX, file=paste0(outPath,'EX_Yly.rda')) }
####  

####  
# Identify a few major countries
iCntries = c('2', '710', '365') # US, China, Russia
pShapes = rep(1,length(cntries))
pShapes[match(iCntries, cntries)] = 15:17
pSize = rep(.75, length(cntries))
pSize[match(iCntries, cntries)] = 1.25
####  

####  
# static plot
multMarg = 1
xl = c(multMarg*min(EX[1,,]),multMarg*max(EX[1,,]))
yl = c(multMarg*min(EX[2,,]),multMarg*max(EX[2,,]))
lims = range(c(xl,yl))

pdf(file=paste0(graphicsPath, "latentPositions_inColor.pdf"))
# par(mar=c(2,2,4,2))
for(tt in 1:TT){	
  plot(0,xlim=lims,ylim=lims,xlab="",ylab="", xaxt="n",yaxt="n",
       main="Posterior mean latent positions showing temporal dynamics")
  if(tt>1){
    arrows(EX[1,tt-1,],EX[2,tt-1,],EX[1,tt,],EX[2,tt,],length=0.05,col='grey', lwd=.5)
  }
  points(t(EX[,tt,]), pch=pShapes, cex=pSize, col=ccols)
  points(t(EX[,tt,iCntries]), pch=pShapes, cex=pSize, col=ccols)
#  if(tt==1) textxy(EX[1,tt,],EX[2,tt,],labs=c(1:26)[-21]) #load "calibrate" package
  if(tt<TT) par(new=TRUE)
}
dev.off()
####  