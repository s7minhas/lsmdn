//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' get get predicted values for Y with count data
//' @param Ex matrix of expected positions for X at time T
//' @param sig2 vector of values for sigma^2 from different iterations
//' @param x1 matrix of values for the first dimension of X from different iterations
//' @param x2 matrix of values for the second dimention of X from different iterations
//' @param Bin  vector of values for betaIn from different iterations
//' @param Bout  vector of values for betaOut from different iterations
//' @param ww  matrix of weights/radii from different iterations
//' @return Matrix with expected values of Y
//' @export
// [[Rcpp::export]]


Rcpp::NumericMatrix predictionPoisson(
	Rcpp::NumericMatrix Ex, Rcpp::NumericVector sig2, 
  Rcpp::NumericMatrix x1, Rcpp::NumericMatrix x2,
	Rcpp::NumericVector Bin, Rcpp::NumericVector Bout, 
  Rcpp::NumericMatrix ww
	) {
  
  int nrows = Ex.nrow();
  int nIter = x1.ncol();
  Rcpp::NumericMatrix yhat(nrows,nrows);
  double dx, tempyhat, sumX, tempX1, tempX2;
  
  
  for(int i = 0; i < nrows-1; i++)
{
  for( int j = i+1; j < nrows; j++)
{
  sumX=0;
  dx = sqrt( pow(Ex(i,0)-Ex(j,0),2) + pow(Ex(i,1)-Ex(j,1),2) );
  for(int l = 0; l < nIter; l++)
{
  
  tempX1 = 1/(2*PI*sig2(l))*exp(-1/(2*sig2(l))*(pow(Ex(i,0)-x1(i,l),2)+pow(Ex(i,1)-x2(i,l),2)) );
  tempX2 = 1/(2*PI*sig2(l))*exp(-1/(2*sig2(l))*(pow(Ex(j,0)-x1(j,l),2)+pow(Ex(j,1)-x2(j,l),2)) );
  sumX = sumX + tempX1*tempX2/nIter;
  
  tempyhat = tempX1*tempX2/nIter*exp(Bin(l)*(dx/ww(j,l)-1)+Bout(l)*(dx/ww(i,l)-1));
  yhat(i,j) = yhat(i,j) + tempyhat;
  tempyhat = tempX1*tempX2/nIter*exp(Bin(l)*(dx/ww(i,l)-1)+Bout(l)*(dx/ww(j,l)-1));
  yhat(j,i) = yhat(j,i) + tempyhat;
}
  yhat(i,j) = yhat(i,j)/sumX;
  yhat(j,i) = yhat(j,i)/sumX;
}
}
  
  return(yhat);
}