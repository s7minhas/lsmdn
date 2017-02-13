//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' update lat positions, beta params, and accrate for non-negative continuous data
//' @param Xitm1 an n x p x T array of the current latent coordinates, where the second dimension is the number dimensions of the latent space, and the third is time 
//' @param dims vector of dimensions of Xitm1
//' @param tunex Variance of the normal random walk proposal for the latent space
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param tuneBIO Variance of the normal random walk proposal for betaIn and betaOut
//' @param ww vector of radius/weights
//' @param t2 variance of initial latent positions
//' @param s2 variance of change in latent positions
//' @param g2 variance of eta
//' @param xiBin mean of prior for betaIn
//' @param xiBout mean of prior for betaOut
//' @param nuBin variance of prior for betaIn
//' @param nuBout variance of prior for betaOut
//' @param Cauchy use a cauchy proposal or not
//' @param rnormsVec Vector of noise for random walks for X
//' @param rnormsBIO Vector of noise for random walks for betaIn/betaOut
//' @return returns list of:
//' \item{Xnew}{New array of latent positions}
//' \item{BinNew}{New values for $\beta_{IN}$}
//' \item{BoutNew}{New values for $\beta_{OUT}}
//' \item{AccRate}{Updated acceptance rate}
//' @export
// [[Rcpp::export]]

List updateNonNegNorm(
  arma::cube Xitm1, Rcpp::IntegerVector dims, double tunex, arma::cube Y,
  double BIN, double BOUT, double tuneBIO,
  arma::colvec ww, double t2, double s2, double g2,
  double xiBin, double xiBout, double nuBin,
  double nuBout, int Cauchy,
  Rcpp::NumericVector rnormsVec, arma::colvec rnormsBIO
  ) {
  
  arma::cube Xold(Xitm1);
  arma::cube Xnew = arma::zeros(dims[1],dims[2],dims[0]);
  
  for(int i=0;i<dims(0);i++) {
    Xnew.slice(i)=Xold.slice(i);
  }
  
  arma::mat Xconc = arma::zeros(dims(0)*dims(2),dims(1)); 
  arma::colvec centerings = arma::zeros(dims(1),1);
  
  arma::cube rnorms(rnormsVec.begin(),dims(1),dims(2),dims(0));
  
  double BinNew =0, BoutNew =0;
  
  double AccProb =0;
  double dz=0, dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  arma::colvec AccRate = arma::zeros(dims(0)*dims(2)+3,1);
  
  /*-------------------- Latent Positions-------------------*/
  
  for(int tt=0;tt < dims[2]; tt++)
{
  for(int i=0;i<dims[0];i++)
{
  AccProb=0;
  if(Cauchy<0.5)
{
  /*---Normal Random Walk---*/
  //  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*arma::randn(dims(1),1);
  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*rnorms.slice(i).col(tt);
}else{
  /*---Cauchy Random Walk---*/
  for(int ell=0;ell<dims(1);ell++)
{
  uu = arma::randu();
  Xnew.slice(i)(ell,tt) = Xold.slice(i)(ell,tt) + tunex*tan(PI*(uu-0.5));
}
}
  
  
  for(int j=0;j<dims[0];j++)
{
  if(j != i){
  dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
  dx = arma::norm(Xold.slice(i).col(tt)-Xold.slice(j).col(tt),2);
  if(Y.slice(tt)(i,j) == 0){
  AccProb += log(1 - 0.5 * erfc(-1*(BIN*(1 - dz/ww(i)) + BOUT*(1 - dz/ww(j)))/sqrt(g2)*M_SQRT1_2)); 
  AccProb += - log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)))/sqrt(g2)*M_SQRT1_2));

  }
  if(Y.slice(tt)(j,i) == 0){
  AccProb += log(1 - 0.5 * erfc(-(BIN*(1 - dz/ww(j)) + BOUT*(1 - dz/ww(i)))/sqrt(g2)*M_SQRT1_2)) ;
  AccProb += - log(1 - 0.5 * erfc(-(BIN*(1 - dx/ww(j)) + BOUT*(1 - dx/ww(i)))/sqrt(g2)*M_SQRT1_2));
  }
  if(Y.slice(tt)(i,j) > 0){
  AccProb += -1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dz/ww(i)) + BOUT*(1 - dz/ww(j)) ),2)/g2 ;
  AccProb += 1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)) ),2)/g2 ;
  }
  if(Y.slice(tt)(j,i) > 0){
  AccProb += -1/2*pow(Y.slice(tt)(j,i) - (BIN*( 1 - dz/ww(j)) + BOUT*(1 - dz/ww(i)) ),2)/g2 ;
  AccProb += 1/2*pow(Y.slice(tt)(j,i) - (BIN*( 1 - dx/ww(j)) + BOUT*(1 - dx/ww(i)) ),2)/g2  ;
}
}
}
  if(tt==0)
{
  insides = trans(Xnew.slice(i).col(tt))*(Xnew.slice(i).col(tt))/t2;
  AccProb += -0.5*insides(0,0);
  insides = trans(Xold.slice(i).col(tt))*(Xold.slice(i).col(tt))/t2;
  AccProb -= -0.5*insides(0,0);
}
  if(tt>0)
{
  muit = Xnew.slice(i).col(tt-1);
  insides = trans(Xnew.slice(i).col(tt)-muit)*(Xnew.slice(i).col(tt)-muit)/s2;
  AccProb += -0.5*insides(0,0);
  insides = trans(Xold.slice(i).col(tt)-muit)*(Xold.slice(i).col(tt)-muit)/s2;
  AccProb -= -0.5*insides(0,0);  
}
  if(tt <dims[2]-1)
{
  muit = Xnew.slice(i).col(tt);
  insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
  AccProb += -0.5*insides(0,0);
  muit = Xold.slice(i).col(tt);
  insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
  AccProb -= -0.5*insides(0,0);  
}
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate(3+tt*dims(0)+i) = 1;
}else
{
  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt);
}
  
}
}
  
  /*---Centering---*/
  for(int i=0;i<dims(0);i++)
{
  for(int tt=0;tt<dims(2);tt++)
{
  Xconc.row(i*dims(2)+tt) = trans(Xnew.slice(i).col(tt));
}
}
  for(int ell=0;ell<dims(1);ell++)
{
  centerings(ell) = sum(Xconc.col(ell))/(dims(0)*dims(2));
}
  for(int i=0;i<dims(0);i++)
{
  for(int tt=0;tt<dims(2);tt++)
{
  for(int ell=0;ell<dims(1);ell++)
{
  Xnew.slice(i)(ell,tt) = Xnew.slice(i)(ell,tt) - centerings(ell);
}
}
}
  
  
  /*-------------------- BetaIn and BetaOut-------------------*/
  AccProb=0;
  if(Cauchy<0.5)
{
  //  BinNew = BIN + tuneBIO*arma::randn();
  BinNew = BIN + tuneBIO*rnormsBIO(0);
}else{
  uu = arma::randu();
  BinNew = BIN + tuneBIO*tan(PI*(uu-0.5));
}
  
  for(int tt=0;tt<dims(2);tt++)
{
  for(int i = 0; i < dims(0); i++)
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
  if(Y.slice(tt)(i,j) == 0){
  AccProb += log(1 - 0.5 * erfc(-1*(BinNew*(1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)))/sqrt(g2)*M_SQRT1_2)); 
  AccProb += - log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)))/sqrt(g2)*M_SQRT1_2));

  }
  if(Y.slice(tt)(i,j) > 0){
  AccProb += -1/2*pow(Y.slice(tt)(i,j) - (BinNew*( 1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)) ),2)/g2 ;
  AccProb += 1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)) ),2)/g2 ;
  }
}
}
}
}
  
  AccProb += -0.5*(BinNew-nuBin)*(BinNew-nuBin)/xiBin;
  AccProb -= -0.5*(BIN-nuBin)*(BIN-nuBin)/xiBin;
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate(0) = 1;
}else
{
  BinNew = BIN;
}
  
  AccProb=0;
  if(Cauchy<0.5)
{
  //  BoutNew = BOUT + tuneBIO*arma::randn();
  BoutNew = BOUT + tuneBIO*rnormsBIO(1);
}else{
  uu = arma::randu();
  BoutNew = BOUT + tuneBIO*tan(PI*(uu-0.5));
}
  
  
  for(int tt=0;tt<dims(2);tt++)
{
  for(int i = 0; i < dims(0); i++)
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
  if(Y.slice(tt)(i,j) == 0){
  AccProb += log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/ww(i)) + BoutNew*(1 - dx/ww(j)))/sqrt(g2)*M_SQRT1_2)); 
  AccProb += - log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)))/sqrt(g2)*M_SQRT1_2));

  }
  if(Y.slice(tt)(i,j) > 0){
  AccProb += -1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/ww(i)) + BoutNew*(1 - dx/ww(j)) ),2)/g2 ;
  AccProb += 1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)) ),2)/g2 ;
  }
}
}
}
}
  
  AccProb += -0.5*(BoutNew-nuBout)*(BoutNew-nuBout)/xiBout;
  AccProb -= -0.5*(BOUT-nuBout)*(BOUT-nuBout)/xiBout;
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate(1) = 1;
}else
{
  BoutNew = BOUT;
}
  
return(Rcpp::List::create(
  Rcpp::Named("X")=Xnew,
  Rcpp::Named("betaIn")=BinNew,
  Rcpp::Named("betaOut")=BoutNew,
  // Rcpp::Named("lambda")=lambNew, # max comment these out once incorp exo into rest of update fns
  // Rcpp::Named("alpha")=alpha, # max comment these out once incorp exo into rest of update fns
  Rcpp::Named("accRate")=AccRate
  )); 
}