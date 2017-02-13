//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Draw from Metropolis Hasting acceptance probability for Gamma^2, for continuous data
//' @param X  an n x p x T array of latent coordinates, where the second dimension is the number dimensions of the latent space, and the third is time 
//' @param dims vector of dimensions of X
//' @param Y an n x n x T array of relational matrices, where the third dimension corresponds to different time periods
//' @param BIN betaIn value
//' @param BOUT betaOut value
//' @param alph Shape parameter for Gamma^2 value
//' @param bta Scale parameter for Gamma^2 value
//' @param ww vector of weights/radius
//' @param g2 current value of Gamma^2
//' @param g2new proposed new value of Gamma^2
//' @return returns new value for Gamma^2 and updated acceptance rate:
//' \item{g2New}{new g2 value}
//' \item{AccRate}{updated accRate}
//' @export
// [[Rcpp::export]]

List gammaAccProbGaussian(
  arma::cube X, Rcpp::IntegerVector dims, arma::cube Y,
  double BIN, double BOUT, double alph,
  double bta, arma::colvec ww, double g2, double g2new
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
  AccProb += -1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)) ),2)/g2new ;
  AccProb += 1/2*pow(Y.slice(tt)(i,j) - (BIN*( 1 - dx/ww(i)) + BOUT*(1 - dx/ww(j)) ),2)/g2 ;
}
}
}
}
  
  for(int i =0;i<dims(0);i++)
{
  AccProb += - (alph + 1)*(log(g2new/g2)) - bta*(1/g2new - 1/g2);
}
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate = 1;
}else
{
  g2new = g2;
}
   
return(Rcpp::List::create(
  Rcpp::Named("g2")=g2new,
  Rcpp::Named("accRate")=AccRate
  ));

}