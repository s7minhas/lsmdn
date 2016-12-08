//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' update lat positions, beta params, and accrate for binomial data using log-likelihood approximation
//' @param Xitm1 data cube
//' @param dims vector of dims
//' @param tunex
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param tuneBIO
//' @param ww vector of weights
//' @param t2
//' @param s2
//' @param xiBin
//' @param xiBout
//' @param nuBin
//' @param nuBout
//' @param Cauchy
//' @param rnormsVec vector
//' @param rnormsBIO vector
//' @param ELout array
//' @param ELin array
//' @param subseq matrix
//' @param degr array
//' @return returns list of:
//' \item{Xnew}
//' \item{BinNew}
//' \item{BoutNew}
//' \item{AccRate}
//' @export lsmdn
// [[Rcpp::export]]

List update_llApprox(
  arma::cube Xitm1, arma::vec dims, double tunex, arma::cube Y,
  double BIN, double BOUT, double tuneBIO,
  arma::colvec ww, double t2, double s2,
  double xiBin, double xiBout, double nuBin,
  double nuBout, int Cauchy,
  Rcpp::NumericVector rnormsVec, arma::colvec rnormsBIO,
  arma::cube ELout, arma::cube ELin, 
  Rcpp::IntegerMatrix subseq, arma::cube degr
  ) {

  arma::cube Xold(Xitm1);
  arma::cube Xnew = arma::zeros(dims[1],dims[2],dims[0]);
  for(int i=0;i<dims(0);i++) {
	  Xnew.slice(i)=Xold.slice(i);
	}
  arma::mat Xconc = arma::zeros(dims(0)*dims(2),dims(1)); 
  arma::colvec centerings = arma::zeros(dims(1),1);

  arma::cube rnorms(rnormsVec.begin(),dims(1),dims(2),dims(0));
  
  int n0 = subseq.ncol(); 
  
  double BinNew =0, BoutNew =0;
  
  double AccProb =0, contrOld=0, contrNew=0;
  double dz=0, dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  arma::colvec AccRate = arma::zeros(dims(0)*dims(2)+3,1);

//-------------------- Latent Positions-------------------

  for(int tt=0;tt < dims[2]; tt++)
{
  for(int i=0;i<dims[0];i++)
{
  AccProb=0;
if(Cauchy<0.5)
{
//---Normal Random Walk---
//  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*arma::randn(dims(1),1);
  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*rnorms.slice(i).col(tt);
}else{
//---Cauchy Random Walk---
  for(int ell=0;ell<dims(1);ell++)
  {
    uu = arma::randu();
    Xnew.slice(i)(ell,tt) = Xold.slice(i)(ell,tt) + tunex*tan(PI*(uu-0.5));
  }
}

//IN edges
if(degr(i,0,tt)>0)
{
for(int j=0;j<degr(i,0,tt);j++)
{
  dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELin(i,j,tt)-1).col(tt),2);
  dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(ELin(i,j,tt)-1).col(tt),2);
  AccProb += (dx-dz)*(BIN/ww(i)+BOUT/ww(ELin(i,j,tt)-1));
}
}
//OUT edges
if(degr(i,1,tt)>0)
{
for(int j=0;j<degr(i,1,tt);j++)
{
  dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
  dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
  AccProb += (dx-dz)*(BIN/ww(ELout(i,j,tt)-1)+BOUT/ww(i));
}
}
//Control estimate

contrOld=0;contrNew=0;
for(int j=0;j<n0;j++)
{
  dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
  dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
  contrNew += log(1+exp(BIN*(1-dz/ww(i))+BOUT*(1-dz/ww(subseq(i,j)-1))))+
              log(1+exp(BIN*(1-dz/ww(subseq(i,j)-1))+BOUT*(1-dz/ww(i))));
  contrOld += log(1+exp(BIN*(1-dx/ww(i))+BOUT*(1-dx/ww(subseq(i,j)-1))))+
              log(1+exp(BIN*(1-dx/ww(subseq(i,j)-1))+BOUT*(1-dx/ww(i))));
}
  AccProb += (dims(0)-1)/n0*(contrOld-contrNew);


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

//OUT edges
if(degr(i,1,tt)>0)
{
for(int j=0;j<degr(i,1,tt);j++)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
  AccProb += (BinNew-BIN)*(1-dx/ww(ELout(i,j,tt)-1));
}
}
//Control estimate
contrOld=0;contrNew=0;
for(int j=0;j<n0;j++)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
  contrNew += log(1+exp(BinNew*(1-dx/ww(subseq(i,j)-1))+BOUT*(1-dx/ww(i))));
  contrOld += log(1+exp(BIN*(1-dx/ww(subseq(i,j)-1))   +BOUT*(1-dx/ww(i))));
}
  AccProb +=(dims(0)-1)/n0*(contrOld-contrNew);

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

//OUT edges
if(degr(i,1,tt)>0)
{
for(int j=0;j<degr(i,1,tt);j++)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
  AccProb += (BoutNew-BOUT)*(1-dx/ww(i));
}
}
//Control estimate
contrOld=0;contrNew=0;
for(int j=0;j<n0;j++)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
  contrNew += log(1+exp(BinNew*(1-dx/ww(subseq(i,j)-1))+BoutNew*(1-dx/ww(i))));
  contrOld += log(1+exp(BinNew*(1-dx/ww(subseq(i,j)-1))+BOUT*(1-dx/ww(i))));
}
  AccProb +=(dims(0)-1)/n0*(contrOld-contrNew);

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

  return(Rcpp::List::create(Xnew,BinNew,BoutNew,AccRate));
}