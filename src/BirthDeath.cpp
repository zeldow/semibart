#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include <RcppArmadillo.h>
#include <Rmath.h>


#include "global.h"
#include "List1.h"
#include "Node.h"

#include "Rlob.h"
#include "Funs.h"
#include "Likelihood.h"
#include "Prior.h"
#include "BirthDeath.h"

using namespace Rcpp;
using namespace arma;

double PBirth(Node *top,Node **n,double *prob)
//returns prob of a birth step, used in BirthDeath
//top: top node of tree
//n: the drawn bottom node
//prob: prob of drawn bottom node
{
  //Rprintf("PBirth\n");
	double PB;
	int can = DrNode(top,n,prob); //checks if birth is possible

	if (!can){
	  PB=0; // if birth is not possible, prob(birth) = 0
	} else if (top->Bot) {
	  PB=1; //if top node is also bottom, prob(birth) =1 (no death possible)
	} else {
	  PB = .5; //else 50/50
	}
	
	//Rprintf("PB: %.1f\n",PB);
	return PB;
}


double DrNogNode(Node *top,Node **n)
//returns prob of a nog node
//top: top node of tree
//n:  on exit , drawn nog node
{

  //Rprintf("DrNogNode\n");
	int nnog,NodeI;
	NodeP *nogvec;
	MakeNogVec(top,&nogvec,&nnog); 

	double u=unif_rand();
	NodeI =  (int)floor(u*nnog)+1;
	*n = nogvec[NodeI];

	delete [] nogvec;

	return 1.0/(double)nnog;

}


//LOOP
int DrNode(Node *top,Node **node, double *nprob)
//draw from from the set of bottom nodes such that a birth is possible
//returns 0 is no bots can grow, 1 otherwise
//top: top of tree
//node: pointer to drawn bottom node
//nprob: prob of drawn bottom node
{
  //Rprintf("DrNode\n");
	int can;
	int nbot;
	int index;
	NodeP *botvec;
	MakeBotVec(top,&botvec,&nbot); // get vector of bottom nodes

	double *probs = new double [nbot+1];
	double sum=0;
	int i;
	for(i=1;i<=nbot;i++) {
		probs[i] = NodeProb(botvec[i]);
		sum += probs[i];
	}

	if(sum>0) {
		can=1;
		for(i=1;i<=nbot;i++) probs[i]=probs[i]/sum;
		index=Disc(probs);
		*node = botvec[index];
		*nprob = probs[index];
	} else {
		can=0;
		*nprob=0;
	}

	delete [] probs;
	delete [] botvec;

	//Rprintf("can: %d\n",can);

	return can;


}

//LOOP
double PrBotNode(Node *top,Node *node)
// called in BirthDeath
// returns:
//	prob of drawing the bottom node as a birth node
{
  //Rprintf("PrBotNode\n");
	double PrNode=-1;
	
	int nbot;
	
	NodeP *botvec;
	MakeBotVec(top,&botvec,&nbot); // get vector of bottom nodes

	double *probs = new double [nbot+1];
	double sum=0;
	int i;
	for(i=1;i<=nbot;i++) {
		probs[i] = NodeProb(botvec[i]);
		sum += probs[i];
	}

	for(i=1;i<=nbot;i++) probs[i]=probs[i]/sum;
	for(i=1;i<=nbot;i++) {
		if(botvec[i] == node) {
			PrNode = probs[i];
			i=nbot+1;
		}
	}

	if(PrNode==-1) Rprintf("error in PrBotNode: node not found in botnode list\n");

	delete [] probs;
	delete [] botvec;

	return PrNode;
}






double NodeProb(Node *n)
//called in DrNode and PrBotNode
//returns:
//	 up to proportionality prob of drawing the node for birthstep
{
  //Rprintf("NodeProb\n");
	int sum = SumGoodVar(n);
	if(sum>0) {
		return 1;
	} else {
		return 0;
	}
}


void KillChildren(Node *nog)
//just like it says
{
  //Rprintf("KillChildren\n");
	delete nog->LeftC;
	delete nog->RightC;
	nog->Bot=1;  //now bottom node
	nog->Nog=0;

	(nog->rule).deall();

	if (!(nog->Top)){
		Node *P = nog->Parent;

		if(nog==P->LeftC) {
			if((P->RightC)->Bot) P->Nog=1;
		} else {
			if((P->LeftC)->Bot) P->Nog=1;
		}
	}
}

double BirthDeath(Node *top,int *BD,int *Done)
//does either a birth or death step
//top: top of tree
//BD: on exit, 1 if birth , 0 if death
//Done: on exit, 1 if step taken , 0 otherwise
{
  //Rprintf("BirthDeath\n");
  //Rcout << "OH SHIT 2    " << y2 << std::endl;
	double PGn,Pbot,PBx,PGl,PGr,PDy,Pnog;
	double PDx,PBy;
	double temprob;


	int VarI;
	int LeftEx,RightEx;
	double alpha1,alpha2,alpha;


	double Ly,Lx;

	Rule *rule=new Rule;

	Node *n,*tempnode;

	PBx=PBirth(top,&n,&Pbot);

	
	
	
	if(Bern(PBx)) {
	
		*BD=1;

		
		PGn = PGrow(n);	
	
		Lx = LogLT(n,top);

		VarI = DrPriVar(n);//draw variable
		DrPriRule(VarI,n,&LeftEx,&RightEx);//draw rule
		SpawnChildren(n,LeftEx,RightEx); //create children

		// this omitted because we will implement it at the metrop level
		/*if ((((n->LeftC)->DataList).length < 5) || (((n->RightC)->DataList).length < 5)){
			// back out if we have less than 5 obs in either new kid; clean up first
			KillChildren(n);
			return -1;
		}*/

		PGl = PGrow(n->LeftC);
		PGr = PGrow(n->RightC);
		Ly = LogLT(n,top);
		
		Pnog = 1.0/((double)(top->NumNogNodes()));
	
		PDy = 1.0-PBirth(top,&tempnode,&temprob);

		alpha1 = (PGn*(1.0-PGl)*(1.0-PGr)*PDy*Pnog)/((1.0-PGn)*PBx*Pbot);
		alpha2 = alpha1*exp(Ly-Lx);
		alpha = min(1.0,alpha2);
			
		if(Bern(alpha)) {		
			*Done=1;
		} else {
			
			KillChildren(n);
			*Done=0;
		}
	} else {
			
			*BD=0;
			PDx=1-PBx;
			Pnog = DrNogNode(top,&n);
			PGl = PGrow(n->LeftC);
			PGr = PGrow(n->RightC);
			
			Lx = LogLT(n,top);

			CopyRule(&(n->rule),rule);
			LeftEx=1-((n->LeftC)->VarAvail[(n->rule).Var]);
			RightEx=1-((n->RightC)->VarAvail[(n->rule).Var]);

			KillChildren(n);
		
			Ly = LogLT(n,top);
			PBy = PBirth(top,&tempnode,&temprob);
			PGn = PGrow(n);
			Pbot=PrBotNode(top,n);
			alpha1 =((1.0-PGn)*PBy*Pbot)/(PGn*(1.0-PGl)*(1.0-PGr)*PDx*Pnog);
			alpha2 = alpha1*exp(Ly-Lx);
			alpha = min(1,alpha2);
			
			if(Bern(alpha)) {
				
				*Done=1;
			} else {
				//put back rule and children
				
				CopyRule(rule,&(n->rule));
				SpawnChildren(n,LeftEx,RightEx);
				
				*Done=0;

			}
			
	}
	delete rule;

	return alpha;

	
}



