#include <RcppArmadillo.h>
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::cube cMissing(
  arma::cube X, arma::vec dims, arma::vec MM, 
  arma::cube Y, int ttt, 
  double BIN, double BOUT, arma::colvec ww
  ) {

  int ttt = ttt-1;  
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
