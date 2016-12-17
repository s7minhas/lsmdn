//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' initialize beta values for binomial data
//' @param X n x p x T array of initial values for the latent space
//' @param dims vector of dimentions for X
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param Xscale 1/n where n is the number of actors
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param ww vector of radiuses/weights
//' @export
// [[Rcpp::export]]

double initializeBinom(
  arma::cube X, Rcpp::IntegerVector dims, arma::cube Y, double Xscale, 
  double BIN, double BOUT, arma::colvec ww
  ) {
  
  double ret =0,dx=0, eta=0;
  
  /*---------------------------------------*/
  
  for(int tt=0;tt<dims(2);tt++)
{
  for(int i = 0; i < dims(0); i++)
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = Xscale*arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
  eta = (BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i)));
  ret += Y.slice(tt)(i,j)*eta-
              log(1+exp(eta));
}
}
}
}
  
return(ret);
}