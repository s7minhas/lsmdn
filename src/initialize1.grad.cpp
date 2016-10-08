c.initialize1.grad <-  cxxfunction(
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

  double dx=0, eta=0;
  Rcpp::NumericVector ret(3);
  
  /*---------------------------------------*/
  
  for(int tt=0;tt<dims(2);tt++)
{
  for(int i = 0; i < dims(0); i++)
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
  eta = BIN*(1-Xscale*dx/ww(j))+BOUT*(1-Xscale*dx/ww(i));
  ret(0) = ret(0)+ dx*(BIN/ww(j)+BOUT/ww(i))*(1/(1+exp(-eta))-Y(i,j,tt));
  ret(1) = ret(1)+ (1-Xscale*dx/ww(j))*(Y(i,j,tt)-1/(1+exp(-eta)));
  ret(2) = ret(2)+ (1-Xscale*dx/ww(i))*(Y(i,j,tt)-1/(1+exp(-eta)));
}
}
}
}
  
  return wrap(ret);
  
  ', plugin="RcppArmadillo")