//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' update weights and accprob for binomial data using log-likelihood approximation
//' @param X data cube
//' @param dims vector of dims
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param tuneW
//' @param wwOld old vector of weights
//' @param wwNew new vector of weights
//' @param ELout array
//' @param ELin array
//' @param subseq matrix
//' @param degr array
//' @return returns list of:
//' \item{wwNew}
//' \item{AccRate}
//' @export lsmdn
// [[Rcpp::export]]

List wAccProb_llApprox(
	arma::cube X, arma::vec dims, arma::cube Y, 
	double BIN, double BOUT, double tuneW,
	arma::colvec wwOld, arma::colvec wwNew,
	arma::cube ELout, arma::cube ELin, arma::mat subseq, 
	arma::cube degr
	) {
  
  double AccProb =0, contrOld=0, contrNew=0, dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  int AccRate = 0;
  
  int n0 = subseq.ncol(); 
  
  /*---------------------------------------*/
  
for(int tt=0;tt<dims(2);tt++)
{
for(int i = 0; i < dims(0); i++)
{

//OUT edges
if(degr(i,1,tt)>0)
{
for(int j=0;j<degr(i,1,tt);j++)
{
  dx = arma::norm(X.slice(i).col(tt)-X.slice(ELout(i,j,tt)-1).col(tt),2);
  AccProb += dx*(BIN*(1/wwOld(ELout(i,j,tt)-1)-1/wwNew(ELout(i,j,tt)-1)) + 
             BOUT*(1/wwOld(i)-1/wwNew(i)));
}
}
//Control estimate

contrOld=0;contrNew=0;
for(int j=0;j<n0;j++)
{
  dx = arma::norm(X.slice(i).col(tt)-X.slice(subseq(i,j)-1).col(tt),2);
  contrNew += log(1+exp(BIN*(1-dx/wwNew(subseq(i,j)-1))+BOUT*(1-dx/wwNew(i))));
  contrOld += log(1+exp(BIN*(1-dx/wwOld(subseq(i,j)-1))+BOUT*(1-dx/wwOld(i))));
}
  AccProb +=(dims(0)-1)/n0*(contrOld-contrNew);

}
}

  
  for(int i =0;i<dims(0);i++)
{
  AccProb += (tuneW*wwNew(i)-1)*log(wwOld(i))- (tuneW*wwOld(i)-1)*log(wwNew(i)) -
  (tuneW*wwNew(i)-0.5)*log(tuneW*wwNew(i))-tuneW*wwNew(i) +
  (tuneW*wwOld(i)-0.5)*log(tuneW*wwOld(i))+tuneW*wwOld(i);
}
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate = 1;
}else
{
  wwNew = wwOld;
}
   
  return(Rcpp::List::create(wwNew,AccRate));
}