#ifndef GUARD_Lib
#define GUARD_Lib

#include <iostream>
#include <string>
#include <vector>
#include <cstring>

#include <RcppArmadillo.h>
#include <Rmath.h>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cfloat>

typedef double *dp;
typedef int *ip;
typedef std::vector<double> Vec;
typedef std::vector<int> IVec;

using namespace Rcpp;
using namespace arma;

class Lib
{
public:
        //constructors, destructor-----------------------------------
        Lib() {}
        ~Lib() {}
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        // static public methods-----------------------------------------
        //static double** almat(long n, long m);
        //static void dealmat(double** mat);
        //static int** ialmat(long n, long m);
        //static void idealmat(int** mat);

     
        //static void acov(Vec& x,int nl, Vec& acov, bool cor=true);
        //static void acov(vec& x,int nl, vec& acov, bool cor=true);
        //static double tssd(vec& x,int n, int nl);
        //static void batchMeans(vec& x, int bsize, vec& means);

       	//static int Disc(Vec& p);
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

};


#endif // GUARD_Lib
