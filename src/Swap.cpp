#include <stdio.h>
#include <iostream>
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
#include "Swap.h"

//LOOP
int CheckCatRule(Node *n,int VarI, int *cats)
{
  //Rprintf("CheckCatRule\n");

	int i;

	int RN = RuleNum(VarI);

	if(n->Bot) {
		delete [] cats;
		return 1;
	} else {

		int *catsl = new int [RN+1];
		int *catsr = new int [RN+1];
		for(i=1;i<=RN;i++) {
			catsl[i]=cats[i];
			catsr[i]=cats[i];
		}

		if(((n->rule).Var)==VarI) {
			for(i=1;i<=RN;i++) {
				if(cats[i]) {
					if((n->rule).CatRule[i]) {
						catsl[i]=0;
					} else {
						catsr[i]=0;
					}
				}
			}
		}
		
		delete [] cats;
		if((!ISum(RN,catsl)) || (!ISum(RN,catsr))) {
			delete [] catsl;
			delete [] catsr;
			
			return 0;
		} else {
			
			if(!CheckCatRule(n->LeftC,VarI,catsl)) {
				return 0;
			} else if(!CheckCatRule(n->RightC,VarI,catsr)) {
				return 0;
			} else {
				return 1;
			}
			
		}
	}
}

					
int CheckOrdRule(Node *n, int VarI, int left, int right)
{
  //Rprintf("CheckOrdRule\n");
	if(n->Bot) {
		return 1;
	} else {

		int rv = (n->rule).Var;
		if(rv==VarI) {
			int or1 = (n->rule).OrdRule;
			if((or1>=left) && (or1<=right)) {
				if(!CheckOrdRule(n->LeftC,VarI,left,or1-1)) {
					return 0;
				} else if (!CheckOrdRule(n->RightC,VarI,or1+1,right)) {
					return 0;
				} else {
					return 1;
				}
			} else {
				return 0;
			}

		} else {
			if(!CheckOrdRule(n->LeftC,VarI,left,right)) {
				return 0;
			} else if (!CheckOrdRule(n->RightC,VarI,left,right)) {
				return 0;
			} else {
				return 1;
			}
		}
	}


}

int CheckRule(Node *n,int VarI)
//starting at node n, check rules using VarI to see if they make sense
{
  //Rprintf("CheckRule\n");
  if(VarType(VarI)==CAT) {
	  int *cats = new int [RuleNum(VarI)+1];
		GetSetCats(n,VarI,cats);
		return CheckCatRule(n,VarI,cats);
	} else {
		int LeftI,RightI;
		GetSplitInterval(&LeftI,&RightI,n,VarI);
		return CheckOrdRule(n,VarI,LeftI,RightI);
	}
}

//LOOP
int AreRulesEqual(Rule *r1,Rule *r2)
{
  //Rprintf("AreRulesEqual\n");
	int i;
	
	if((r1->Var)==(r2->Var)) {
	  if(VarType(r1->Var)==CAT) {
		  for(i=1;i<=RuleNum(r1->Var);i++) if(r1->CatRule[i]!=r2->CatRule[i]) return 0;
			return 1;
		} else {
			if((r1->OrdRule)!=(r2->OrdRule)) {
				return 0;
			} else {
				return 1;
			}
		}
	} else {
		return 0;
	}
}

double SwapRule(Node *top,int *Done)
// step which tries swapping rules
{
  //Rprintf("SwapRule\n");

	double alpha=0.0; // note backout = (alpha = -1)
	double u; //uniform (0,1)

	//first get the list of nodes with swappable rules
	int Nswap,dadVarI,kidVarI,checkrule;
	NodeP *swapvec;
	MakeSwapVec(top,&swapvec,&Nswap);


	//if there are no swappable rule back out
	if(!Nswap) {
		alpha =  -1.0;
	} else {


	// randomly choose a node with a swappable rule = dad
	//u=ran1(&idum);
        u = unif_rand();
	int NodeI =  (int)floor(u*Nswap)+1;
	Node *dad = swapvec[NodeI];
	
	//check whether children have the same rule

	int SameRule = AreRulesEqual(&((dad->LeftC)->rule),&((dad->RightC)->rule));

	if(!SameRule) {

		//find out which children have rules and pick one
		int lI=0,rI=0;
		Node *kid;
		if((dad->LeftC)->rule.Var) lI=1;
		if((dad->RightC)->rule.Var) rI=1;
		if(!(lI+rI)) Rprintf("error in SwapRule: neither child of dad has a rule\n");
	
		if((lI+rI)==2) {
			//u=ran1(&idum);
			u=unif_rand();
			if(u<.5) {
				kid = dad->LeftC;
			} else {
				kid = dad->RightC;
			}
		} else if(lI) {
			kid = dad->LeftC;
		} else {
			kid = dad->RightC;
		}

	
		Rule KidRule(kid->rule); //store kid's rule
		CopyRule(&(dad->rule),&(kid->rule)); //copy dadrule onto kidrule
		CopyRule(&KidRule,&(dad->rule)); //copy kidrule onto dadrule
	
		//check if rule make sense after given the swap+++++++++++++++++++++++++++
		dadVarI = (dad->rule).Var;
		kidVarI = (kid->rule).Var;

		checkrule = CheckRule(dad,dadVarI);
		if((dadVarI!=kidVarI) && checkrule) checkrule = CheckRule(dad,kidVarI);
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		//put the rules back
		CopyRule(&(kid->rule),&(dad->rule));
		CopyRule(&KidRule,&(kid->rule));
	
		
	
		//if the swap was ok (rules made sense)
		if(checkrule) {

				//get logpri and logL from current tree (X)
				double XLogPi = LogPriT(top);
				double XLogL = LogLT(dad,top);

				CopyRule(&(dad->rule),&(kid->rule)); //copy dadrule onto kidrule
				CopyRule(&KidRule,&(dad->rule)); //copy kidrule onto dadrule

				//fix data at nodes below cnode given new rule
				FixDataBelow(dad);
			
				//  fix VarAvail
				dadVarI = (dad->rule).Var;
				kidVarI = (kid->rule).Var;
				UpDateVarAvail(dad,dadVarI);
				if(!(dadVarI==kidVarI)) UpDateVarAvail(dad,kidVarI);

				//get logpri and logL from current tree (X)
				double YLogPi = LogPriT(top);
				double YLogL = LogLT(dad,top);

				alpha = min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));
				if(Bern(alpha)) {
					*Done = 1;
				} else {
					CopyRule(&(kid->rule),&(dad->rule));
					CopyRule(&KidRule,&(kid->rule));

					//fix data at nodes below cnode given new rule
					FixDataBelow(dad);
			
					//  fix VarAvail
					dadVarI = (dad->rule).Var;
					kidVarI = (kid->rule).Var;
					UpDateVarAvail(dad,dadVarI);
					if(!(dadVarI==kidVarI)) UpDateVarAvail(dad,kidVarI);
					*Done = 0;
				}

			} else {
				alpha =  -1; //not a legal swap	
			}
	} else {

		//get prior of current tree
		double XLogPi = LogPriT(top);
		//next line is new!
		double XLogL = LogLT(dad,top);
	

		//swap the rules (for the case where the two children are the same
		CopyRule(&(dad->rule),&((dad->RightC)->rule));
		CopyRule(&((dad->LeftC)->rule),&(dad->rule));
		CopyRule(&((dad->RightC)->rule),&((dad->LeftC)->rule));

		//check if rule is ok
		dadVarI = (dad->rule).Var;
		kidVarI = ((dad->LeftC)->rule).Var;

		checkrule = CheckRule(dad,dadVarI);
		if((dadVarI!=kidVarI) && checkrule) checkrule = CheckRule(dad,kidVarI);
		if(checkrule) {


			//fix data at nodes below cnode given new rule
			FixDataBelow(dad);
		
			//  fix VarAvail
			kidVarI = ((dad->LeftC)->rule).Var;
			dadVarI = (dad->rule).Var;
			UpDateVarAvail(dad,dadVarI);
			if(!(dadVarI==kidVarI)) UpDateVarAvail(dad,kidVarI);

			//get prior at new tree
			double YLogPi = LogPriT(top);
			//next line is new!
			double YLogL = LogLT(dad,top);
			
			// next line was old version of the code, didn't check log likelihood
			// and I think it should because grandkids will be affected...
			//alpha = min(1.0,exp(YLogPi-XLogPi));

			// next line is the new version.
			alpha = min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));
				
	
			if(Bern(alpha)) {
				*Done = 1;
			} else {
				CopyRule(&(dad->rule),&((dad->RightC)->rule));
				CopyRule(&((dad->LeftC)->rule),&(dad->rule));
				CopyRule(&((dad->RightC)->rule),&((dad->LeftC)->rule));

				//fix data at nodes below cnode given new rule
				FixDataBelow(dad);
			
				//  fix VarAvail
				kidVarI = ((dad->LeftC)->rule).Var;
				dadVarI = (dad->rule).Var;
				UpDateVarAvail(dad,dadVarI);
				if(!(dadVarI==kidVarI)) UpDateVarAvail(dad,kidVarI);
				*Done = 0;
			}
		} else {
			//checkrule failed, swap back
			CopyRule(&(dad->rule),&((dad->RightC)->rule));
			CopyRule(&((dad->LeftC)->rule),&(dad->rule));
			CopyRule(&((dad->RightC)->rule),&((dad->LeftC)->rule));
			alpha=-1;
			*Done=0;
		}
	}
	}
	delete [] swapvec;
	return alpha;
}

