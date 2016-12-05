#include <RcppArmadillo.h>
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double cInitialize1(
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