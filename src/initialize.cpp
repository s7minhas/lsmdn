//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' initialize
//' @param X data cube
//' @param dims vector of dims
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param Xscale
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param ww vector of weights
//' @export initBetaInOut
// [[Rcpp::export]]

double initialize(
  arma::cube X, arma::vec dims, arma::cube Y, double Xscale, 
  double BIN, double BOUT, arma::colvec ww
  ) {
  
  double ret =0,dx=0, eta=0;
  
  /*---------------------------------------*/
  
  for(int tt=0; tt<dims(2); tt++) {
    for(int i = 0; i < dims(0); i++) {
      for(int j = 0; j < dims(0); j++) {        
        if(i != j) {
          dx = Xscale*arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
          eta = (BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i)));
          ret += Y.slice(tt)(i,j)*eta-log(1+exp(eta));
        }
      }
    }
  }
  
return(ret);
}