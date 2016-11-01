#ifndef GUARD_Beta
#define GUARD_Beta

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

class Beta
{
 public:
  //do i need constructor/destructors?
  vec getBeta() {return this->draw;}
  void setBeta(vec draw) {this->draw = draw;}
  void setPrior(double b, vec p) {
    //Rprintf("setPrior\n");
    this->sig_beta = b*b;
    this->mean_beta = p;
  }
  void setAA(mat AA) {this->AA = AA; }
  void setAy(mat Ay) {this->Ay = Ay; }
  void setSig(double a) {this->sig_reg = a*a;}
  void drawPost() {
    //Rprintf("Beta::drawPost\n");
    mat I_A(NumA,NumA,fill::eye);
    //mat lamN = sig_reg*sig_beta*(AA/sig_reg + I_A/sig_beta);
    mat lamN = (AA*sig_beta + I_A*sig_reg)/(sig_beta*sig_reg);
    vec muN = lamN.i()*(Ay*sig_beta+mean_beta*sig_reg)/(sig_reg*sig_beta);
    //mat Ydraw = randn(1,NumA);
    mat Ydraw(1,NumA);
    for(int i = 0; i < NumA; i++) Ydraw(0,i) = R::rnorm(0.0,1.0);
    mat finalsig = lamN.i();
    draw = muN + trans(Ydraw*chol(finalsig));
    }

 private:
  //prior on Betas
    double sig_beta; //this times identity
  vec mean_beta; //needs to be vector
  //parameter
  vec draw; //needs to be vector
  //data
  double sig_reg;
  mat AA;
  mat Ay;
};
#endif
