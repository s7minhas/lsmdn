c.prediction <- cxxfunction(
  signature(EX="numeric",SIG2="numeric",X1T="numeric",X2T="numeric",BIN="numeric",BOUT="numeric",WW="numeric"),
  body='

  Rcpp::NumericMatrix Ex(EX);
  Rcpp::NumericMatrix x1(X1T);
  Rcpp::NumericMatrix x2(X2T);
  Rcpp::NumericMatrix ww(WW);
  Rcpp::NumericVector sig2(SIG2);
  Rcpp::NumericVector Bin(BIN);
  Rcpp::NumericVector Bout(BOUT);
  
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
  
  return yhat;
  ',plugin="Rcpp")