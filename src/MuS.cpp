#include <iostream>
#include <cmath>
#include "MuS.h"
#include "Lib.h"

#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

// public methods -------------------------------------------------
void MuS::drawPost()   
{
  // Rprintf("MuS::drawPost\n");
   mu = post_m + post_s*norm_rand();
}
double MuS::getLogILik()   //probably ok
{
  //   Rprintf("MuS::getLogILik\n");
   double res = .5*log(a/(a+b));
   res -= s2/(2.0*sigma2);
   res -= .5*(a*b*ybar*ybar)/(a+b);
   return res;  
}
//LOOP
void MuS::updatepost()  //change index in loop (0 to <nob ??) DONE
{
  // Rprintf("MuS::updatepost\n");
   int i;
   double d;
   if(nob) {
     /* ybar = mean(y(indices));
     s2 = var(y(indices)) * indices.size(); //check
     */
     ybar=0.0;
     for(i=1;i<=nob;i++) {
       //Rprintf("index: %d\n",indices(i)-1);
       //Rcout<<indices<<std::endl;
       //Rcout<<y<<std::endl;
       ybar += y(indices(i)-1);
     }
      ybar /= nob;
      s2=0.0;
      for(i=1;i<=nob;i++) {d=y(indices(i)-1)-ybar; s2 += d*d;}
      b = nob/sigma2; 
      post_m = (b*ybar)/(a+b);
      post_s = 1.0/std::sqrt(a+b);
   }
   else {
      post_m=0.0;
      post_s= 1.0/std::sqrt(a);
      b=0.0;
   }
}
void MuS::setData(int nob, mat x, vec y, //probably ok
                  ivec indices, vec w)
{
  //Rprintf("MuS::setData\n");
    //Rprintf("\t nob = %d\n",nob);
    //Rcout << trans(w) << std::endl;
    //Rcout << trans(indices) << std::endl;
    //Rcout <<trans(y) <<std::endl;
   this->nob=nob;
   this->y=y;
   this->indices=indices;
   //Rcout << nob << std::endl;
   //Rcout << trans(y) << std::endl;
   updatepost();
}

//LOOP
vec MuS::getFits(int np, mat xpred, ivec indpred)  //change index on loop?
{
  // Rprintf("MuS::getFits\n");
  //double *rv = new double[np+1];
  vec rv(np+1);
  for(int i=1;i<=np;i++) rv(i)=post_m;
   return rv;
}

