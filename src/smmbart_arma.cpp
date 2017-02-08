// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <Rcpp.h>
// #include <Rcpp/RNGScope.h>
// #include <Rcpp/routines.h>
#include <Rmath.h>

#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <search.h>
#include <vector>
#include <valarray>
#include <algorithm>

#include "global.h"
#include "Node.h"
#include "Funs.h"
#include "Prior.h"
#include "MuS.h"
#include "Sdev.h"
#include "Beta.h"
#include "List1.h"
#include "EndNodeModel.h"
#include "Lib.h"


using namespace Rcpp;
using namespace arma;

//set global vars

double pBD, pChange, pSwap;  //probabilities of each step

double musig;

double **RuleMat;
int *Ivec;

int NumX; // number of x variables
int NumY; //number of y variables
int NumA;

double kfac;

ivec VarType;
int NumObs; // number of observations
mat X; // x data, note: cats are double
mat trt;
vec y;
vec ydat;
vec ydat1;

vec logitweights; //weights for logistic regression

double nu = 8; //for approximation to logistic regression 


vec weights;

ivec RuleNum; // integer vec of length NumX, ith is number of split
				// points for ORD var and number of CATs for CAT var




// typedefs
typedef double *dp;
typedef Node *NodeP;

//make sure these are necessary
CPriParams PriParams;
EndNodeModel* endNodeModel=0;


// [[Rcpp::export]]
List semibartcpp(arma::mat iX, arma::mat itrt, arma::vec iy,
	     double sigma, int sigdf, double sigquant, double kfac,
	     double power, double base, 
	     arma::vec meanb, double sigb,
	     int ntree, int ndpost, arma::ivec inumcut,
	     int iusequants,
	     double binary_offset, int probitlink,
	     int verbose, int printevery)
{

  X = iX; trt = itrt;
  bool binary = (binary_offset > -1000.0);
  if(binary)
    sigma = 1.0;
  
  if(verbose) {
    if(binary)
      Rprintf("\n\nRunning BART with binary y\n\n");
    else
      Rprintf("\n\nRunning BART with numeric y\n\n");
  }
  mat AtA = trt.t()*trt; 

 
  NumX = X.n_cols;
  NumA = trt.n_cols;
  NumObs = X.n_rows;

  //define globally
  VarType.zeros(NumX+1); 
  RuleNum.zeros(NumX+1); //# split pts for ord and # cats for cat
  for(int i=1;i<=NumX;i++) {
    VarType(i) = ORD;
      }

  //initialize logitweights to ones
  logitweights.ones(NumObs);
  
  bool usequants = true;
   if(!(iusequants)) usequants=false;
  
  RuleMat = new dp [NumX + 1];
  int cnq,cfac,cnc,cind,coffset;
  vec xcol, sortedx; //Vector to store specific columns of X??

  double maxx,minx,xinc;
  if(usequants) {
    for(int i=1;i<=NumX;i++) {
      xcol = X.col(i-1);
      sortedx = unique(xcol);
      cnq = sortedx.n_elem;
      if(cnq<=(inumcut(i-1)+1)) {
	cfac = 1;
	cnc = cnq-1;
	coffset=0;
      } else {
	cnc = inumcut(i-1);
	cfac = cnq/cnc;
	coffset = (cfac/2);
      }
      RuleNum(i) = cnc; 
      RuleMat[i] = new double[cnc+1];
      for(int j=0;j<cnc;j++) {
	cind = std::min(j*cfac+coffset,cnq-2);
	RuleMat[i][j+1] = (sortedx(cind)+sortedx(cind+1))/2.0;
      }
    } 
  } else {
    for(int i=1;i<=NumX;i++) {
      xcol = X.col(i-1);
      maxx = xcol.max();
      minx = xcol.min();
      RuleNum(i) = inumcut(i-1);
      xinc = ( maxx - minx ) / ( RuleNum(i) + 1 );
      RuleMat[i] = new double[RuleNum(i)];
      for(int j=1;j<RuleNum(i);j++) RuleMat[i][j] = minx + j*xinc;
    }
  }
  
  Ivec = new int[NumObs+1];
   for(int i=1;i<=NumObs;i++) Ivec[i]=i;

   if(binary) {
     ydat = iy;
     y = 2.0*(iy+0.5) - 1.0 - binary_offset;
   } else {
     y = iy;
   }
  
  if(binary)
    musig = 3.0/(kfac*sqrt(ntree));
  else 
    musig = 0.5/(kfac*sqrt(ntree));


  PriParams.base = base;
  PriParams.power = power;

  double lambda;
  double qchi = R::qchisq(1.0-sigquant,sigdf,1,0);
  lambda = (sigma*sigma*qchi)/sigdf;

  //need to define
  // sets up priors on the end node parameters
  MuS mu;
  mu.setSigma(sigma);
  mu.setPriorS(musig);
  endNodeModel = &mu;
  
  // sets up prior on regression sigma
  Sdev sd;
  sd.setPrior(sigdf, lambda); //define these

  // sets up the prior on the blip parameters
  Beta bet;
  bet.setPrior(sigb, meanb);
  bet.setAA(AtA);
  bet.setSig(sigma);
  
  // set probabilities of each tree proposal
  pBD = 0.5;
  pChange = 0.4;
  pSwap = 0.1;
  
  std::vector<Node*> theTrees(ntree);  //vector of trees -- change double* to Node*
  typedef std::vector<Node*>::size_type nvs;
  for(nvs i=0;i<theTrees.size();i++)
    {theTrees[i] = new Node; theTrees[i]->SetData();} //change double to Node
  
  mat mtrainFits(ntree,NumObs); mtrainFits.zeros();
  //Check what these are for
  int Done = 0;
  double alpha = 0.0;
  int step = 0;
  //end check
  
  //this is the total fit for the trees
  vec mtotalfit(NumObs); mtotalfit.zeros();
  vec mfits(NumObs); mfits.zeros();

  vec eps(NumObs); eps.zeros();//residuals for backfitting

  // Storage for MCMC parameters of interest
  // stores beta MCMC parameters
  mat betaReps(ndpost,NumA); betaReps.zeros();
  // stores sigma MCMC valyes
  vec sigmaReps(ndpost); sigmaReps.zeros();
  
  // stores current values of Beta
  vec curr_beta(NumA); curr_beta.zeros();
  // double curr_sig;
  

  mat Aty; //for storing trt.t()*y (or eps)
  
  //MCMC loop
  for(int k=0;k<ndpost;k++) { // loop through MCMC ndpost times
    
 //   if ( Progress::check_abort() ) {
//      return List::create(_["sigmaReps"] = sigmaReps,
//                          _["betaReps"] = betaReps);
//    }

    for(nvs i=0;i<theTrees.size();i++) {    // loop through each tree
      
      ydat1 = y; //set eps to be y (temp storage)
      ydat1 = ydat1 - mtotalfit; // subtract current total fit from trees
      for(int j=0;j<NumA;j++) // subtract current fit from blip
	    {
	      ydat1 = ydat1 - trt.col(j) * curr_beta(j); // likely better way for this
	    }
      ydat1 = ydat1 + trans(mtrainFits.row(i)); //add back in fit for tree of interest
      //will eventually refit this
      
      // run metrop step
      alpha = Metrop(&theTrees[i],&Done,&step); //unclear how this works
      theTrees[i]->currentFits(&mu,NumObs,X,ydat1,weights,mfits); //need other params?

      // get new fits
      mtotalfit = mtotalfit - trans(mtrainFits.row(i)); //subtract old fits
      mtotalfit = mtotalfit + mfits; //add new fits
      mtrainFits.row(i) = trans(mfits);
    }  //end loop through trees

    
    //update Beta
    //subtract mtotalfit from y;
    eps = y;
    eps = eps - mtotalfit;
  
    if(binary && probitlink==1) {  //logit link
      Aty = trt.t()*diagmat(logitweights)*eps;
      bet.setAy(Aty);
      AtA = trt.t()*diagmat(logitweights)*trt;
      bet.setAA(AtA);
      bet.drawPost();
      curr_beta = bet.getBeta();  //divide to get back on scale of logit link
      betaReps.row(k) = trans(curr_beta);
    } else {  //continuous or probit link
      Aty = trt.t()*eps;
      bet.setAy(Aty);
      bet.drawPost();
      curr_beta = bet.getBeta();
      betaReps.row(k) = trans(curr_beta);
    }
    
    if(!binary) {
      //update Sigma
      eps = y;
      eps = eps - mtotalfit;
      for(int j=0;j<NumA;j++) //subtract off blip fits
	    {
	      eps = eps - trt.col(j) * curr_beta(j); //maybe better way to code this
	    }
      sd.setData(NumObs,eps); 
      sd.drawPost();
      mu.setSigma(sd.getS());
      bet.setSig(sd.getS());
      sigmaReps(k) = sd.getS();
    }

    if(binary) {
      if(probitlink==0) {
        // double u, Z, blip;
	double Z, blip;
        for(int i = 1; i <= NumObs; i++) {
	        blip = 0.0;
	        for(int j=0;j<NumA;j++)
	          {
	            blip = blip+trt(i-1,j)*curr_beta(j); 
	          }
	        //u = unif_rand();
	        // u = R::runif(0.0,1.0);
	        if(ydat(i-1) > 0) {
	          //Z = R::qnorm((1.0-u)*R::pnorm(-mtotalfit(i-1)-blip-binary_offset,0.0,1.0,1,0)+u,0.0,1.0,1,0);
	          Z = R::qnorm(R::runif(R::pnorm(0,mtotalfit(i-1)+blip+binary_offset,1.0,1,0),1),mtotalfit(i-1)+blip+binary_offset,1.0,1,0);
	        } else {
	          //Z = -R::qnorm((1.0-u)*R::pnorm(mtotalfit(i-1)+blip+binary_offset,0.0,1.0,1,0)+u,0.0,1.0,1,0);
	          Z = R::qnorm(R::runif(0,R::pnorm(0,mtotalfit(i-1)+blip+binary_offset,1.0,1,0)),mtotalfit(i-1)+blip+binary_offset,1.0,1,0);
	        }
	        //y(i-1) = mtotalfit(i-1) + blip + Z;
	        y(i-1) = Z;
        } 
      } else if(probitlink==1) {  //t8 approx for logit
        double Z, blip;
        for(int i = 0; i < NumObs; i++) {
          blip = 0.0;
          for(int j=0;j<NumA;j++)
          {
            blip = blip+trt(i,j)*curr_beta(j); 
          }
          
          if(ydat(i) > 0) {
            Z = R::qnorm(R::runif(R::pnorm(0,mtotalfit(i)+blip+binary_offset,sqrt(1.0/logitweights(i)),1,0),1),mtotalfit(i)+blip+binary_offset,sqrt(1.0/logitweights(i)),1,0);
          } else {
            Z = R::qnorm(R::runif(0,R::pnorm(0,mtotalfit(i)+blip+binary_offset,sqrt(1.0/logitweights(i)),1,0)),mtotalfit(i)+blip+binary_offset,sqrt(1.0/logitweights(i)),1,0);
          }
          y(i) = Z;
          
          logitweights(i) = R::rgamma( (nu + 1.0 ) / 2 , 2 / ( nu + pow(y(i) - mtotalfit(i) - blip - binary_offset, 2) ) );
        }
      }
    }
    
    if(verbose) {
      if(((k+1)%printevery==0)) Rprintf("iteration: %d (of %d)\n",k+1,ndpost);
    }

  }
    
  
  delete [] Ivec;
  for(int i=1;i<=NumX;i++) delete [] RuleMat[i];
  delete [] RuleMat;
  for(nvs i=0;i<theTrees.size();i++)
    theTrees[i]->deall();

  if(verbose) {
    Rprintf("Function call finished.\n");
  }

   return List::create(_["sigmaReps"] = sigmaReps,
		      _["betaReps"] = betaReps);
  
}
