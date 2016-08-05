#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <fstream>

#include "Lib.h"


#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;


// NOTE: Vec = vector of doubles
// NOTE: IVec = vector of integers

//public methods-----------------------------------------------------
//this makes some matrix ??
// deleted mean and sdev function -- bring back if needed


//CHANGED SEE ORIGINAL IF NEEDED
//void Lib::acov(Vec& x,int nl, Vec& acov, bool cor)
/*void Lib::acov(vec& x,int nl, vec& acov, bool cor)
{
        double c;
        int n = x.size();
        double m = mean(x);
        acov.reset(); //removes all elements leaving vector with 0 elements
	acov.insert_cols(0,nl+1);
        for(int i=0;i<=nl;i++) //change to armadillo
        {
                c=0.0;
                for(int j=0;j<(n-i);j++) c += (x(j)-m)*(x(j+i)-m);
                acov(i) = c;
        }
        if(cor)
        {
	  double c0 = acov(0);
	  acov = acov / c0;
        }
        else
	   acov =  acov/n;

}

double Lib::tssd(vec& x, int n, int nl)
{
        vec gamma;
        acov(x,nl,gamma,false);
        double v=gamma(0);
        for(int i=1;i<=nl;i++) v += 2*(1.0-((double)i/n))*gamma(i);
        return sqrt(v/n);
}
*/

/*
// Don't think this is used
// The function in Funs.cpp seems to be though
int Lib::Disc(vec& p)
// draw from discrete distributin given by p, return index
{
    double sum;
    double u = unif_rand();

    int i=0;
    //sum=p[0];
    sum=p(0);
    while(sum<u) {
        i += 1;
        //sum += p[i];
	sum += p(i);
    }
    return i;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
