#include <RcppArmadillo.h>
#include "global.h"
#include "List1.h"
#include "Node.h"
#include "Rlob.h"
#include "Funs.h"
#include "Likelihood.h"
#include "Prior.h"
#include "BirthDeath.h"
#include "ChangeRule.h"
#include "Swap.h"

using namespace Rcpp;
using namespace arma;

double LogLNode(Node *n)
//returns of integrated likelihood of data at node n
{	
  //Rprintf("LogLNode\n");
	int nob;
	//int *divec;
	ivec divec;
	double lp=0.0;

	//Rprintf("\tnob in LogLNode: %d\n",nob);
	
	//MakeIntVec(&(n->DataList),&divec,&nob);
	MakeIntVec(&(n->DataList),divec,&nob);
	//Rprintf("\tnob in LogLNode: %d\n",nob);

	//Rcout << trans(divec) << std::endl;

	if(nob==0) lp = -10000000.0;
	else
	{
	  //Rcout << "TEST1" << nob << std::endl;
	  //Rcout << "TEST2" << trans(y) << std::endl;
	  
	  // inputs data
	  //endNodeModel->setData(nob,XDatR,YDat1,divec,weights); 
	  endNodeModel->setData(nob,X,ydat1,divec,weights); 
	  // outputs Log integrated likelihood
	  lp = endNodeModel->getLogILik(); //from MuS.cpp??
	}

	//delete [] divec;

	return lp;
}

//LOOP
double LogLT(Node *branch,Node *top)
{
  //Rprintf("LogLT\n");
	int nbot;
	NodeP *botvec;
	//think this is from Funs.cpp
	MakeBotVec(branch,&botvec,&nbot);
	//Rprintf("\tnbot in LogLT: %d\n",nbot);
	//Rcout << botvec << std::endl;
	int i;
	double LP=0.0;

	for(i=1;i<=nbot;i++) {

		LP += LogLNode(botvec[i]);	
	}

	delete [] botvec;
	return LP;
}
