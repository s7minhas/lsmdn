c.WAccProb2 <-  cxxfunction(
  signature(Data="numeric",DIMS="integer",Yy="numeric",
            BETAIN="numeric",BETAOUT="numeric",TUNEW="numeric",
            WWOld="numeric",WWNew="numeric",
            ELOUT="integer",ELIN="integer",SUBSEQ="numeric",DEG="integer"),
  body='
/*Dims is c(n,p,TT,dinmax,doutmax)
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
  
    double AccProb =0, contrOld=0, contrNew=0, dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  int AccRate = 0;
  
  Rcpp::NumericVector ELinVec(ELIN);
  arma::cube ELin(ELinVec.begin(),dims(0),dims(3),dims(2));
  Rcpp::NumericVector ELoutVec(ELOUT);
  arma::cube ELout(ELoutVec.begin(),dims(0),dims(4),dims(2));
  Rcpp::IntegerMatrix subseq(SUBSEQ);
  int n0 = subseq.ncol(); 
  Rcpp::NumericVector DegrVec(DEG);
  arma::cube degr(DegrVec.begin(),dims(0),2,dims(2));

  
  /*---------------------------------------*/
  

for(int tt=0;tt<dims(2);tt++)
{
for(int i = 0; i < dims(0); i++)
{

//OUT edges
if(degr(i,1,tt)>0)
{
for(int j=0;j<degr(i,1,tt);j++)
{
  dx = arma::norm(X.slice(i).col(tt)-X.slice(ELout(i,j,tt)-1).col(tt),2);
  AccProb += dx*(BIN*(1/wwOld(ELout(i,j,tt)-1)-1/wwNew(ELout(i,j,tt)-1)) + 
             BOUT*(1/wwOld(i)-1/wwNew(i)));
}
}
//Control estimate

contrOld=0;contrNew=0;
for(int j=0;j<n0;j++)
{
  dx = arma::norm(X.slice(i).col(tt)-X.slice(subseq(i,j)-1).col(tt),2);
  contrNew += log(1+exp(BIN*(1-dx/wwNew(subseq(i,j)-1))+BOUT*(1-dx/wwNew(i))));
  contrOld += log(1+exp(BIN*(1-dx/wwOld(subseq(i,j)-1))+BOUT*(1-dx/wwOld(i))));
}
  AccProb +=(dims(0)-1)/n0*(contrOld-contrNew);

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