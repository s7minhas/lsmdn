#include <RcppArmadillo.h>
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec cPrediction(
	arma::mat Ex, arma::vec sig2, arma::mat x1, arma::mat x2,
	arma::vec Bin, arma::vec Bout, arma::mat ww
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
  sumX = sumX + tempX1*tempX2/10000;
  
  tempyhat = tempX1*tempX2/10000/(1+exp(Bin(l)*(dx/ww(j,l)-1)+Bout(l)*(dx/ww(i,l)-1)));
  yhat(i,j) = yhat(i,j) + tempyhat;
  tempyhat = tempX1*tempX2/10000/(1+exp(Bin(l)*(dx/ww(i,l)-1)+Bout(l)*(dx/ww(j,l)-1)));
  yhat(j,i) = yhat(j,i) + tempyhat;
}
  yhat(i,j) = yhat(i,j)/sumX;
  yhat(j,i) = yhat(j,i)/sumX;
}
}
  
  return(yhat);
}