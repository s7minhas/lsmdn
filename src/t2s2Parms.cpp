//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' generate shape and scale parameters to update t2 and s2 for all types of data
//' @param X  an n x p x T array of latent coordinates, where the second dimension is the number dimensions of the latent space, and the third is time 
//' @param dims vector of dimensions of X
//' @param thetaT shape parameter for t
//' @param thetaS shape parameter for s
//' @param phiT scale parameter for t
//' @param phiS scale parameter for s
//' @return returns list of:
//' \item{shapeT}{Shape Parameter for tau2}
//' \item{scaleT}{Scale Parameter for tau2}
//' \item{shapeS}{Shape Parameter for sigma2}
//' \item{scaleS}{Scale Parameter for sigma2}
//' @export
// [[Rcpp::export]]

List t2s2Parms(
	arma::cube X, Rcpp::IntegerVector dims, double thetaT, 
	double thetaS, double phiT, double phiS 
	) {
	
  double shapeT=0, scaleT=0, shapeS=0, scaleS=0;  
  arma::mat insides = arma::zeros(1,1);
  
  //---------------------------------------
  
  shapeT = thetaT + 0.5*dims(0)*dims(1);
  scaleT = phiT;
  shapeS = thetaS + 0.5*dims(0)*dims(1)*(dims(2)-1);
  scaleS = phiS;
  for(int i =0;i<dims(0);i++)
{
  insides = 0.5*trans(X.slice(i).col(0))*(X.slice(i).col(0));
  scaleT += insides(0,0);
  for(int tt=1;tt<dims(2);tt++)
{
  insides = 0.5*trans(X.slice(i).col(tt)-X.slice(i).col(tt-1))*
  (X.slice(i).col(tt)-X.slice(i).col(tt-1));
  scaleS += insides(0,0);
}
}
   
 return(Rcpp::List::create(
  Rcpp::Named("shapeT")=shapeT,
  Rcpp::Named("scaleT")=scaleT,
  Rcpp::Named("shapeS")=shapeS,
  Rcpp::Named("scaleS")=scaleS
  ));
 }