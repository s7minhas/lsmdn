c.update2 <-  cxxfunction(
  signature(Xitm1="numeric",DIMS="integer",TUNEX="numeric",Yy="numeric",
            BETAIN="numeric",BETAOUT="numeric",TUNEBIO="numeric",
            WW="numeric",t2X="numeric",s2X="numeric",
            xiBIN="numeric",xiBOUT="numeric",nuBIN="numeric",
            nuBOUT="numeric",CAUCHY="integer",
            RNORMS="numeric",RNORMSBIO="numeric",
            ELOUT="integer",ELIN="integer",SUBSEQ="numeric",DEG="integer"),
  body='
/*Dims is c(n,p,TT,dinmax,doutmax)
Z nxT; X pxTxn; Y pxpxT; 
*/
  Rcpp::IntegerVector dims(DIMS);
  int Cauchy = Rcpp::as<int>(CAUCHY);

  Rcpp::NumericVector XOLD(Xitm1);
  arma::cube Xold(XOLD.begin(),dims[1],dims[2],dims[0]);
  arma::cube Xnew = arma::zeros(dims[1],dims[2],dims[0]);
  for(int i=0;i<dims(0);i++)
{
  Xnew.slice(i)=Xold.slice(i);
}
  arma::mat Xconc = arma::zeros(dims(0)*dims(2),dims(1)); 
  arma::colvec centerings = arma::zeros(dims(1),1);

  Rcpp::NumericVector rnormsVec(RNORMS);
  arma::cube rnorms(rnormsVec.begin(),dims(1),dims(2),dims(0));
  arma::colvec rnormsBIO = Rcpp::as<arma::colvec>(RNORMSBIO);

  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  Rcpp::NumericVector ELinVec(ELIN);
  arma::cube ELin(ELinVec.begin(),dims(0),dims(3),dims(2));
  Rcpp::NumericVector ELoutVec(ELOUT);
  arma::cube ELout(ELoutVec.begin(),dims(0),dims(4),dims(2));
  Rcpp::IntegerMatrix subseq(SUBSEQ);
  int n0 = subseq.ncol(); 
  Rcpp::NumericVector DegrVec(DEG);
  arma::cube degr(DegrVec.begin(),dims(0),2,dims(2));
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  double xiBin = Rcpp::as<double>(xiBIN);
  double xiBout = Rcpp::as<double>(xiBOUT);
  double nuBin = Rcpp::as<double>(nuBIN);
  double nuBout = Rcpp::as<double>(nuBOUT);
  double BinNew =0, BoutNew =0;
  double tunex = Rcpp::as<double>(TUNEX);
  double tuneBIO = Rcpp::as<double>(TUNEBIO);
  
  double t2 = Rcpp::as<double>(t2X);
  double s2 = Rcpp::as<double>(s2X);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  
  double AccProb =0, contrOld=0, contrNew=0;
  double dz=0, dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  arma::colvec AccRate = arma::zeros(dims(0)*dims(2)+3,1);

  

//-------------------- Latent Positions-------------------

  for(int tt=0;tt < dims[2]; tt++)
{
  for(int i=0;i<dims[0];i++)
{
  AccProb=0;
if(Cauchy<0.5)
{
//---Normal Random Walk---
//  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*arma::randn(dims(1),1);
  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*rnorms.slice(i).col(tt);
}else{
//---Cauchy Random Walk---
  for(int ell=0;ell<dims(1);ell++)
  {
    uu = arma::randu();
    Xnew.slice(i)(ell,tt) = Xold.slice(i)(ell,tt) + tunex*tan(PI*(uu-0.5));
  }
}

//IN edges
if(degr(i,0,tt)>0)
{
for(int j=0;j<degr(i,0,tt);j++)
{
  dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELin(i,j,tt)-1).col(tt),2);
  dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(ELin(i,j,tt)-1).col(tt),2);
  AccProb += (dx-dz)*(BIN/ww(i)+BOUT/ww(ELin(i,j,tt)-1));
}
}
//OUT edges
if(degr(i,1,tt)>0)
{
for(int j=0;j<degr(i,1,tt);j++)
{
  dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
  dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
  AccProb += (dx-dz)*(BIN/ww(ELout(i,j,tt)-1)+BOUT/ww(i));
}
}
//Control estimate

contrOld=0;contrNew=0;
for(int j=0;j<n0;j++)
{
  dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
  dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
  contrNew += log(1+exp(BIN*(1-dz/ww(i))+BOUT*(1-dz/ww(subseq(i,j)-1))))+
              log(1+exp(BIN*(1-dz/ww(subseq(i,j)-1))+BOUT*(1-dz/ww(i))));
  contrOld += log(1+exp(BIN*(1-dx/ww(i))+BOUT*(1-dx/ww(subseq(i,j)-1))))+
              log(1+exp(BIN*(1-dx/ww(subseq(i,j)-1))+BOUT*(1-dx/ww(i))));
}
  AccProb += (dims(0)-1)/n0*(contrOld-contrNew);


  if(tt==0)
{
  insides = trans(Xnew.slice(i).col(tt))*(Xnew.slice(i).col(tt))/t2;
  AccProb += -0.5*insides(0,0);
  insides = trans(Xold.slice(i).col(tt))*(Xold.slice(i).col(tt))/t2;
  AccProb -= -0.5*insides(0,0);
}
  if(tt>0)
{
  muit = Xnew.slice(i).col(tt-1);
  insides = trans(Xnew.slice(i).col(tt)-muit)*(Xnew.slice(i).col(tt)-muit)/s2;
  AccProb += -0.5*insides(0,0);
  insides = trans(Xold.slice(i).col(tt)-muit)*(Xold.slice(i).col(tt)-muit)/s2;
  AccProb -= -0.5*insides(0,0);  
}
  if(tt <dims[2]-1)
{
  muit = Xnew.slice(i).col(tt);
  insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
  AccProb += -0.5*insides(0,0);
  muit = Xold.slice(i).col(tt);
  insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
  AccProb -= -0.5*insides(0,0);  
}
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate(3+tt*dims(0)+i) = 1;
}else
{
  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt);
}

}
}

/*---Centering---*/
for(int i=0;i<dims(0);i++)
{
for(int tt=0;tt<dims(2);tt++)
{
  Xconc.row(i*dims(2)+tt) = trans(Xnew.slice(i).col(tt));
}
}
for(int ell=0;ell<dims(1);ell++)
{
  centerings(ell) = sum(Xconc.col(ell))/(dims(0)*dims(2));
}
for(int i=0;i<dims(0);i++)
{
for(int tt=0;tt<dims(2);tt++)
{
for(int ell=0;ell<dims(1);ell++)
{
  Xnew.slice(i)(ell,tt) = Xnew.slice(i)(ell,tt) - centerings(ell);
}
}
}


/*-------------------- BetaIn and BetaOut-------------------*/

  AccProb=0;
if(Cauchy<0.5)
{
//  BinNew = BIN + tuneBIO*arma::randn();
  BinNew = BIN + tuneBIO*rnormsBIO(0);
}else{
  uu = arma::randu();
  BinNew = BIN + tuneBIO*tan(PI*(uu-0.5));
}

for(int tt=0;tt<dims(2);tt++)
{
for(int i = 0; i < dims(0); i++)
{

//OUT edges
if(degr(i,1,tt)>0)
{
for(int j=0;j<degr(i,1,tt);j++)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
  AccProb += (BinNew-BIN)*(1-dx/ww(ELout(i,j,tt)-1));
}
}
//Control estimate
contrOld=0;contrNew=0;
for(int j=0;j<n0;j++)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
  contrNew += log(1+exp(BinNew*(1-dx/ww(subseq(i,j)-1))+BOUT*(1-dx/ww(i))));
  contrOld += log(1+exp(BIN*(1-dx/ww(subseq(i,j)-1))   +BOUT*(1-dx/ww(i))));
}
  AccProb +=(dims(0)-1)/n0*(contrOld-contrNew);

}
}

  AccProb += -0.5*(BinNew-nuBin)*(BinNew-nuBin)/xiBin;
  AccProb -= -0.5*(BIN-nuBin)*(BIN-nuBin)/xiBin;
  
  uu= arma::randu();
if(uu<exp(AccProb))
{
  AccRate(0) = 1;
}else
{
  BinNew = BIN;
}
  

  AccProb=0;
if(Cauchy<0.5)
{
//  BoutNew = BOUT + tuneBIO*arma::randn();
  BoutNew = BOUT + tuneBIO*rnormsBIO(1);
}else{
  uu = arma::randu();
  BoutNew = BOUT + tuneBIO*tan(PI*(uu-0.5));
}


for(int tt=0;tt<dims(2);tt++)
{
for(int i = 0; i < dims(0); i++)
{

//OUT edges
if(degr(i,1,tt)>0)
{
for(int j=0;j<degr(i,1,tt);j++)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(ELout(i,j,tt)-1).col(tt),2);
  AccProb += (BoutNew-BOUT)*(1-dx/ww(i));
}
}
//Control estimate
contrOld=0;contrNew=0;
for(int j=0;j<n0;j++)
{
  dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(subseq(i,j)-1).col(tt),2);
  contrNew += log(1+exp(BinNew*(1-dx/ww(subseq(i,j)-1))+BoutNew*(1-dx/ww(i))));
  contrOld += log(1+exp(BinNew*(1-dx/ww(subseq(i,j)-1))+BOUT*(1-dx/ww(i))));
}
  AccProb +=(dims(0)-1)/n0*(contrOld-contrNew);

}
}

  AccProb += -0.5*(BoutNew-nuBout)*(BoutNew-nuBout)/xiBout;
  AccProb -= -0.5*(BOUT-nuBout)*(BOUT-nuBout)/xiBout;
  
  uu= arma::randu();
if(uu<exp(AccProb))
{
  AccRate(1) = 1;
}else
{
  BoutNew = BOUT;
}


  return Rcpp::List::create(Xnew,BinNew,BoutNew,AccRate);

  ', plugin="RcppArmadillo")