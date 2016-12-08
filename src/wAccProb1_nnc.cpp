#include <RcppArmadillo.h>
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List cWAccProb1_nnc(
	arma::cube X, arma::vec dims, arma::cube Y, 
	double BIN, double BOUT, double tuneW,
	arma::colvec wwOld, arma::colvec wwNew, double g2
	) {
  
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
  if(Y.slice(tt)(i,j) == 0){
  AccProb += log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/wwNew(i)) + BOUT*(1 - dx/wwNew(j)))/sqrt(g2)*M_SQRT1_2)); 
  AccProb += - log(1 - 0.5 * erfc(-1*(BIN*(1 - dx/wwOld(i)) + BOUT*(1 - dx/wwOld(j)))/sqrt(g2)*M_SQRT1_2));

  }
  if(Y.slice(tt)(j,i) == 0){
  AccProb += log(1 - 0.5 * erfc(-(BIN*(1 - dx/wwNew(j)) + BOUT*(1 - dx/wwNew(i)))/sqrt(g2)*M_SQRT1_2)) ;
  AccProb += - log(1 - 0.5 * erfc(-(BIN*(1 - dx/wwOld(j)) + BOUT*(1 - dx/wwOld(i)))/sqrt(g2)*M_SQRT1_2));
  }
  if(Y.slice(tt)(i,j) > 0){
  AccProb += -1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/wwNew(i)) + BOUT*(1 - dx/wwNew(j)) ),2)/g2 ;
  AccProb += 1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/wwOld(i)) + BOUT*(1 - dx/wwOld(j)) ),2)/g2 ;
  }
  if(Y.slice(tt)(j,i) > 0){
  AccProb += -1/2*pow(Y.slice(tt)(j,i) - (BIN*( 1 - dx/wwNew(j)) + BOUT*(1 - dx/wwNew(i)) ),2)/g2 ;
  AccProb += 1/2*pow(Y.slice(tt)(j,i) - (BIN*( 1 - dx/wwOld(j)) + BOUT*(1 - dx/wwOld(i)) ),2)/g2  ;
}}
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
   
return(Rcpp::List::create(wwNew,AccRate));
}
