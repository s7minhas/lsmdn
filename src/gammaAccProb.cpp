//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Draw from inv gamma and adjust acceptance rate, for non-negative normal data
//' @param X data cube
//' @param dims vector of dims
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param alph
//' @param bta
//' @param ww vector of weights
//' @param g2 current value
//' @param g2new new value
//' @return returns list of starting values:
//' \item{g2New}{new g2 value}
//' \item{AccRate}{updated accRate}
//' @export
// [[Rcpp::export]]

List gammaAccProb(
  arma::cube X, Rcpp::IntegerVector dims, arma::cube Y,
  double BIN, double BOUT, double alph,
  double bta, arma::colvec ww, double g2, double g2new
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
  AccProb += log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)))/sqrt(g2new)*M_SQRT1_2)); 
  AccProb += - log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)))/sqrt(g2)*M_SQRT1_2));

  }
  if(Y.slice(tt)(i,j) > 0){
  AccProb += -1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)) ),2)/g2new ;
  AccProb += 1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)) ),2)/g2 ;
  }
}
}
}
}
  
  for(int i =0;i<dims(0);i++)
{
  AccProb += - (alph + 1)*(log(g2new/g2)) - bta*(1/g2new - 1/g2);
}
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate = 1;
}else
{
  g2new = g2;
}
   
return(Rcpp::List::create(g2new,AccRate));

}