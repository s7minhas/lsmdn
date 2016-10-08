#include <RcppArmadillo.h>
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::depends(RcppArmadillo)]]

//' Regular Log-logistic likelihood
//' 
//' @param theta Current parameter values. 
//' @param y Dependent variables. 
//' @param X Duration equation covariates.
//' 
//' @keywords internal
// [[Rcpp::export]]

c.initialize1 <-  cxxfunction(
  signature(Data="numeric",DIMS="integer",Yy="numeric",
            XSCALE="numeric",BETAIN="numeric",
            BETAOUT="numeric",WW="numeric"),
  body='
/*Dims is c(n,p,TT)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);

  double Xscale= Rcpp::as<double>(XSCALE);

  double ret =0,dx=0, eta=0;
  
  /*---------------------------------------*/
  
  for(int tt=0;tt<dims(2);tt++)
{
  for(int i = 0; i < dims(0); i++)
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = Xscale*arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
  eta = (BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i)));
  ret += Y.slice(tt)(i,j)*eta-
              log(1+exp(eta));
}
}
}
}
  
  return wrap(ret);
  
  ', plugin="RcppArmadillo")