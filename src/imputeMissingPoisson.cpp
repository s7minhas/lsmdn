//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Impute missingness in Y for Count Data
//' @param X  an n x p x T array of latent coordinates, where the second dimension is the number dimensions of the latent space, and the third is time 
//' @param dims vector of dimensions of X
//' @param MM n x n x T array with 1s where Y is missing
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param Ttt time period value
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param ww vector of weights/radius
//' @return Y with imputed values
//' @export
// [[Rcpp::export]]

arma::cube imputeMissingPoisson(
  arma::cube X, Rcpp::IntegerVector dims, Rcpp::IntegerVector MM, 
  arma::cube Y, int Ttt, 
  double BIN, double alpha, arma::colvec ww, arma::cube WL
  ) {

  int ttt = Ttt-1;  
  double dx=0, uu=0, Prob=0;
  double Mean;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  double BOUT = fabs(alpha) - BIN;

  
  /*---------------------------------------*/
  
  for(int i = 0; i < dims(0); i++) {
    if(std::find(MM.begin(),MM.end(),i) !=MM.end()) {
      for(int j = 0; j < dims(0); j++) {
        if(i != j) {
          dx = arma::norm(X.slice(i).col(ttt)-X.slice(j).col(ttt),2);
          Prob = alpha + WL.slice(ttt)(i,j) + BIN*(-dx/ww(j))+BOUT*(-dx/ww(i));
          Mean = exp(Prob);
          Y(i,j,ttt) = Rcpp::rpois(1, Mean)[0];
        }
      }
    }
  }
  
  return(Y);
}