#ifndef GUARD_Rlob
#define GUARD_Rlob

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>


#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

typedef double *dp;
typedef int *ip;

// these are probably not necessary but need to figure out how to change to arma
//double **almat(long n, long m);
//void dealmat(double **m);
//int **ialmat(long n, long m);
//void idealmat(int **m);

//deleted bunch of print functions

//double ran1(long *idum);
//double gasdev(long *idum);

int Bern(double p);
int Disc(double *p);
double min(double a,double b); //min of a and b
//double myDoubleAbs(double a); //unused

void indtd(int k,int ind,int *d);
//int dtind(int k,int *d); //unused

int ISum(int n, int *Iv); // sums a Vector of integers
double max(double a, double b); // max of a and b

//void choldc(double **a, int n, double p[]);
//void sym_chol_inv(int n,double **a,double **li);
//double sym_inv_det(int n,double **a, double **ai);
//void mul_ltl(int n,double **l,double **a);
//void solve_rtxb(int n,double **r,double *x,double *b);
//void solve_rxb(int n, double **r,double *x,double *b);
//double gammln(double xx);

//int compare( const void *arg1, const void *arg2 ); //unused
//void stanAndSortForCart(int n, int k, double **raw, double **stan, int *numUnique, double **uniqueVals, double* meanV, double* scaleV); //unused

#endif



