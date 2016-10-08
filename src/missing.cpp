c.missing <-  cxxfunction(
  signature(Data="numeric",DIMS="integer",MMM="integer",Yy="numeric",Ttt="integer",
            BETAIN="numeric",BETAOUT="numeric",WW="numeric"),
  body='
/*Dims is c(n,p,TT,K); 
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  int ttt = Rcpp::as<int>(Ttt)-1;
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::IntegerVector MM(MMM);
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  
  double dx=0, uu=0, Prob=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  
  /*---------------------------------------*/
  
  for(int i = 0; i < dims(0); i++)
{
  if(std::find(MM.begin(),MM.end(),i) !=MM.end())
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = arma::norm(X.slice(i).col(ttt)-X.slice(j).col(ttt),2);
  Prob = BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i));
  Prob = 1/(1+exp(-Prob));
  uu= arma::randu();
  if(uu<Prob)
{
  Y(i,j,ttt) = 1;
}else{
  Y(i,j,ttt) = 0;
}
  /*  if(std::find(MM.begin(),MM.end(),j) ==MM.end())
{
  Prob = BIN*(1-dx/ww(i))+BOUT*(1-dx/ww(j));
  Prob = 1/(1+exp(-Prob));
  uu= arma::randu();
  if(uu<Prob)
{
  Y(j,i,ttt) = 1;
}else{
  Y(j,i,ttt) = 0;
}
}*/
  
}
}
}
}
  
  return Rcpp::wrap(Y);
  
  ', plugin="RcppArmadillo")