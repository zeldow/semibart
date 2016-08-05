#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "Rlob.h"


#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

typedef double *dp;
typedef int *ip;



//deleted bunch of print functions

int Bern(double p)
//Bernouill random number generator
{
  //Rprintf("Bern\n");
        if(unif_rand() < p) {
                return 1;
        } else {
                return 0;
        }
}

//LOOP
//this is like multinomial function from Jason MCMC code
int Disc(double *p)
// draw from discrete distributin given by p, return index
{
  //Rprintf("Disc\n");
   double sum;
   double u = unif_rand();

   int i=1;
   sum=p[1];
   while(sum<u) {
      i += 1;
      sum += p[i];
   }
   return i;
}

double min(double a,double b) // used in ChangeRule.cpp I think
{
  //Rprintf("min\n");
	if(a<b) {
		return a;
	} else {
		return b;
	}
}

double  max(double a, double b)
{
  //Rprintf("max\n");
	if(a>b) {
		return a;
	} else {
		return b;
	}
}


//LOOP
void indtd(int k,int ind,int *d)  //used in ChangeRule and Prior
{
  //Rprintf("indtd\n");
int j,i,nind;

nind=ind;

for(i=0;i<k;i++) {
	j=k-i-1;
	d[j+1] = (int)(((double) nind)/(pow(2.0,(double)j)));
	nind=nind-d[j+1]*(int)pow(2.0,(double)j);
}
}


//LOOP
int ISum(int n, int *Iv) //Used in Funs and Swap
{
  //Rprintf("ISum\n");
	int sum=0;
	int i;
	for(i=1;i<=n;i++) sum += Iv[i];
	return sum;
}




