#ifndef GUARD_global
#define GUARD_global

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

#define CAT 1  //  1 is flag for categorical type of variable
#define ORD 2 //   2 is flag for ordered variable

#define BIRTH 1
#define DEATH 2
#define SWAP 3
#define CHANGE 4

extern double pBD;
extern double pSwap;
extern double pChange;

extern int NumX;
extern int NumA;

extern int NumObs;
extern int *Ivec;

extern ivec RuleNum; 
extern double **RuleMat; 
extern ivec VarType;

extern mat X;
extern mat trt;
extern vec y;
extern vec ydat;
extern vec ydat1;
extern vec weights;




#include "CPriParams.h"
#include "EndNodeModel.h"

extern CPriParams PriParams;
extern EndNodeModel* endNodeModel;
#endif

