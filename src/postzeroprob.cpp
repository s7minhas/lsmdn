c.postzeroprob <- cxxfunction(
  signature(Xi1="numeric",Xi2="numeric",Xj1="numeric",Xj2="numeric",SS2="numeric",LAM="numeric",PP0="numeric"),
  body='
  
  
  double lam = Rcpp::as<double>(LAM);
  double pp0 = Rcpp::as<double>(PP0);
  Rcpp::NumericMatrix xi1(Xi1);
  Rcpp::NumericMatrix xi2(Xi2);
  Rcpp::NumericMatrix xj1(Xj1);
  Rcpp::NumericMatrix xj2(Xj2);
  Rcpp::NumericVector ss2(SS2);
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
  
  return Rcpp::wrap(ret);
  ', 
  plugin="Rcpp")