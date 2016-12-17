//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 


//' update weights and accprob for non-negative continuous normal data
//' @param X  an n x p x T array of latent coordinates, where the second dimension is the number dimensions of the latent space, and the third is time 
//' @param dims vector of dimensions of X
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param tuneW variance of the proposal distribution for w
//' @param wwOld old vector of weights
//' @param wwNew new vector of weights
//' @param g2 Variance of eta
//' @return returns list of:
//' \item{wwNew}{New weights/radius}
//' \item{AccRate}{Updated acceptance probability}
//' @export
// [[Rcpp::export]]

List wAccProbNonNegNormal(
	arma::cube X, Rcpp::IntegerVector dims, arma::cube Y, 
	double BIN, double BOUT, double tuneW,
	arma::colvec wwOld, arma::colvec wwNew, double g2
	) {
  
  double AccProb =0,dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  int AccRate = 0;
    
  /*---------------------------------------*/
  
  for(int tt=0;tt<dims(2);tt++)
{
  for(int i = 0; i < dims(0); i++)
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
  if(Y.slice(tt)(i,j) == 0){
  AccProb += log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/wwNew(i)) + BOUT*(1 - dx/wwNew(j)))/sqrt(g2)*M_SQRT1_2)); 
  AccProb += - log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/wwOld(i)) + BOUT*(1 - dx/wwOld(j)))/sqrt(g2)*M_SQRT1_2));

  }
  if(Y.slice(tt)(i,j) > 0){
  AccProb += -1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/wwNew(i)) + BOUT*(1 - dx/wwNew(j)) ),2)/g2 ;
  AccProb += 1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/wwOld(i)) + BOUT*(1 - dx/wwOld(j)) ),2)/g2 ;
  }
}
}
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
