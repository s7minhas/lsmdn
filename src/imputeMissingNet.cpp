//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Impute missingness in Y
//' @param X data cube
//' @param dims vector of dims
//' @param MM 
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param Ttt
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param ww vector of weights
//' @return Y with imputed values
//' @export
// [[Rcpp::export]]

arma::cube imputeMissingNet(
  arma::cube X, Rcpp::IntegerVector dims, Rcpp::IntegerVector MM, 
  arma::cube Y, int Ttt, 
  double BIN, double BOUT, arma::colvec ww
  ) {

  int ttt = Ttt-1;  
  double dx=0, uu=0, Prob=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  
  /*---------------------------------------*/
  
  for(int i = 0; i < dims(0); i++) {
		if(std::find(MM.begin(),MM.end(),i) !=MM.end()) {
			for(int j = 0; j < dims(0); j++) {
				if(i != j) {
					dx = arma::norm(X.slice(i).col(ttt)-X.slice(j).col(ttt),2);
					Prob = BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i));
					Prob = 1/(1+exp(-Prob));
					uu= arma::randu();
					if(uu<Prob) {
						Y(i,j,ttt) = 1;
					} else {
						Y(i,j,ttt) = 0;
					}
  /*  if(std::find(MM.begin(),MM.end(),j) ==MM.end())
{
  Prob = BIN*(1-dx/ww(i))+BOUT*(1-dx/ww(j));
  Prob = 1/(1+exp(-Prob));
  uu= arma::randu();
  if(uu<Prob)
{
  Y(j,i,ttt) = 1;
}else{
  Y(j,i,ttt) = 0;
}
}*/
  
				}
			}
		}
	}
  
  return(Y);
}
