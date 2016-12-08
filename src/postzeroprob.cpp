//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Posterior zero prob
//' @param xi1 matrix
//' @param xi2 matrix
//' @param xj1 matrix
//' @param xj2 matrix
//' @param ss2 vector
//' @param lam
//' @param pp0
//' @return stuff
//' @export
// [[Rcpp::export]]

arma::vec postZeroProb(
	arma::mat xi1, arma::mat xi2, arma::mat xj1, arma::mat xj2,
	arma::vec ss2, double lam, double pp0
	) {

  Rcpp::NumericVector normCDF(1);
  Rcpp::NumericVector temp1(1);
  Rcpp::NumericVector ret(1);
  
  double angle = 0, temp = 0, x = 0, y = 0;
  int nrows = xi1.nrow();
  int ncols = xi1.ncol();
  
  for(int iter = 0; iter < ncols; iter++)
{
  temp = 0;
  for(int time = 1; time < nrows; time++)
{
  x = xj1(time,iter)-xi1(time-1,iter);
  y = xj2(time,iter)-xi2(time-1,iter);
  angle = atan2(y,x);
  temp += (xi1(time,iter)-xi1(time-1,iter))*cos(angle)+(xi2(time,iter)-xi2(time-1,iter))*sin(angle);
}
  temp1(0) = (lam*temp-ss2(iter))/(lam*sqrt(ss2(iter)*(nrows-1)));
  normCDF = pnorm(temp1);
  ret = ret + normCDF*sqrt(2*PI*ss2(iter)/(nrows-1))*exp(pow(lam*temp-ss2(iter),2)/(2*lam*lam*ss2(iter)*(nrows-1)));
}
  ret = ret*(1-pp0)/(ncols*pp0*lam);
  
  return(ret);
}