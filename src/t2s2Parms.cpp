#include <RcppArmadillo.h>
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List cT2s2Parms(
	arma::cube X, arma::vec dims, double thetaT, 
	double thetaS, double phiT, double phiS 
	) {
	
  double shapeT=0, scaleT=0, shapeS=0, scaleS=0;  
  arma::mat insides = arma::zeros(1,1);
  
  //---------------------------------------
  
  shapeT = thetaT + 0.5*dims(0)*dims(1);
  scaleT = phiT;
  shapeS = thetaS + 0.5*dims(0)*dims(1)*(dims(2)-1);
  scaleS = phiS;
  for(int i =0;i<dims(0);i++) {
	  insides = 0.5*trans(X.slice(i).col(0))*(X.slice(i).col(0));
	  scaleT += insides(0,0);
	  for(int tt=1;tt<dims(2);tt++) {
		  insides = 0.5*trans(X.slice(i).col(tt)-X.slice(i).col(tt-1))*
		  (X.slice(i).col(tt)-X.slice(i).col(tt-1));
		  scaleS += insides(0,0);
		}
	}
  
  
 return(Rcpp::List::create(shapeT,scaleT,shapeS,scaleS));
 }
 