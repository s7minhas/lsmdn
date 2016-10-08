c.WAccProb1 <-  cxxfunction(
  signature(Data="numeric",DIMS="integer",Yy="numeric",
            BETAIN="numeric",BETAOUT="numeric",TUNEW="numeric",
            WWOld="numeric",WWNew="numeric"),
  body='
/*Dims is c(n,p,TT,K)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  
  arma::colvec wwOld = Rcpp::as<arma::colvec>(WWOld);
  arma::colvec wwNew = Rcpp::as<arma::colvec>(WWNew);
  double tuneW = Rcpp::as<double>(TUNEW);
  
  double AccProb =0,dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  int AccRate = 0;
  
  
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
  AccProb += Y.slice(tt)(i,j)*dx*(BIN*(1/wwOld(j)-1/wwNew(j)) + 
  BOUT*(1/wwOld(i)-1/wwNew(i))) +
  log(1+exp(BIN*(1-dx/wwOld(j))+BOUT*(1-dx/wwOld(i)))) -
  log(1+exp(BIN*(1-dx/wwNew(j))+BOUT*(1-dx/wwNew(i))));
}
}
}
}
  
  for(int i =0;i<dims(0);i++)
{
  AccProb += (tuneW*wwNew(i)-1)*log(wwOld(i))- (tuneW*wwOld(i)-1)*log(wwNew(i)) -
  (tuneW*wwNew(i)-0.5)*log(tuneW*wwNew(i))-tuneW*wwNew(i) +
  (tuneW*wwOld(i)-0.5)*log(tuneW*wwOld(i))+tuneW*wwOld(i);
}
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate = 1;
}else
{
  wwNew = wwOld;
}
   
  return Rcpp::List::create(wwNew,AccRate);
  
  ', plugin="RcppArmadillo")