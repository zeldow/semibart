#ifndef GUARD_Sdev
#define GUARD_Sdev

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


class Sdev {
public:
  //Sdev();
  //~Sdev();
  double getS() {return this->s;}
  void setS(double s) {this->s = fabs(s);}
  void setPrior(int nu, double lambda) {
    //Rprintf("setPrior\n");
    this->nu = nu;
    this->lambda = lambda;
  }
  void setData(int nob, vec e) {
    //Rprintf("setData\n");
    this->nob = nob;
    this->e = e;
  }
    void drawPost()
  {
    //  Rprintf("Sdev::drawPost\n");
    //double ss = 0.0;
    //for(int i = 0;i<nob;i++) ss += e[i]*e[i];
    double ss = dot(e,e);
  int nupost(nu+nob);
  double nlpost(nu*lambda + ss);
  s = sqrt(nlpost/R::rchisq((double)nupost)); //works with R::rchisq;
  //not with Rcpp::rchisq
  }

 private:
  //parameter
  double s;
  //prior
  int nu;
  double lambda;
  //data
  int nob;
  vec e; 
};
#endif
