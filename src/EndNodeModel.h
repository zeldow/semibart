#ifndef GUARD_EndNodeModel_h
#define GUARD_EndNodeModel_h

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// probably will need to change these to armadillo format
class EndNodeModel
{
	public:
		virtual double getLogILik()  = 0;
		//virtual void setData(int nob,  double** x, double* y, int* indices, double* w)  =0;
		virtual void setData(int nob,  mat x, vec y, ivec indices, vec w)  =0;
		//virtual double* getFits(int np,  double** xpred, int* indpred) =0;
		virtual vec getFits(int np,  mat xpred, ivec indpred) =0;
		virtual int getEstimateDim() const  =0;
		//virtual double* getParameterEstimate() =0;
		virtual ~EndNodeModel() {}
};

#endif


