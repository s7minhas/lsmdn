//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' initialize using gradient
//' @param X data cube
//' @param dims vector of dims
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param Xscale
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param ww vector of weights
//' @export initBetaInOut
// [[Rcpp::export]]

arma::vec initializeGrad(
  arma::cube X, arma::vec dims, arma::cube Y, double Xscale, 
  double BIN, double BOUT, arma::colvec ww  
  ) {

  double dx=0, eta=0;
  Rcpp::NumericVector ret(3);
  
  /*---------------------------------------*/
  
  for(int tt=0; tt<dims(2); tt++) {
    for(int i = 0; i < dims(0); i++) {
      for(int j = 0; j < dims(0); j++) {        
        if(i != j) {
          dx = arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
          eta = BIN*(1-Xscale*dx/ww(j))+BOUT*(1-Xscale*dx/ww(i));
          ret(0) = ret(0)+ dx*(BIN/ww(j)+BOUT/ww(i))*(1/(1+exp(-eta))-Y(i,j,tt));
          ret(1) = ret(1)+ (1-Xscale*dx/ww(j))*(Y(i,j,tt)-1/(1+exp(-eta)));
          ret(2) = ret(2)+ (1-Xscale*dx/ww(i))*(Y(i,j,tt)-1/(1+exp(-eta)));
        }
      }
    }
  }
  
return(ret);
}