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
#include "ChangeRule.h"






void CatFindBots(Node *curr,int VarI,int cat,NodeP *botvec,int *fcount)
//adds 1 to fcount i if category cat associated with variable VarI
// can get from node curr to the ith bottom node
// 
{
  //Rprintf("CatFindBots\n");
	if(curr->Bot) {
		int ind=1;
		while(!(curr==botvec[ind])) ind+=1;
		fcount[ind] += 1;
	} else {
		if((curr->rule).Var==VarI) {
			if((curr->rule).CatRule[cat]) {
				CatFindBots(curr->RightC,VarI,cat,botvec,fcount);
			} else {
				CatFindBots(curr->LeftC,VarI,cat,botvec,fcount);
			}
		} else {
			CatFindBots(curr->RightC,VarI,cat,botvec,fcount);
			CatFindBots(curr->LeftC,VarI,cat,botvec,fcount);
		}
	}
}

void FindGoodOrdRules(Node *n,int VarI, int &l, int &u)
//good rule have splits in [l,u]
{
  //Rprintf("FindGoodOrdRules\n");
	int LeftI,RightI; 
	LeftI = 1; // right value if you top out
	RightI = RuleNum(VarI); // right value if you top out

	int lmin,lmax,rmin,rmax;

	GetSplitInterval(&LeftI,&RightI,n,VarI);

	lmin=RightI+1;
	rmin=RightI+1;
	lmax=LeftI-1;
	rmax=LeftI-1;

	OrdFindMinMax(n->LeftC,VarI,&lmin,&lmax);
	OrdFindMinMax(n->RightC,VarI,&rmin,&rmax);

	l = (int)max(LeftI,lmax+1);
	u = (int)min(RightI,rmin-1);

}

void OrdFindMinMax(Node *n,int VarI, int *min, int *max)
//used to find good ord rule, go down tree adjusting min and max whenever VarI is used
{
  //Rprintf("OrdFindMinMax\n");
  if(VarType(VarI)==CAT) Rprintf("error in OrdFindMinMax, CAT var\n");

	if(!(n->Bot)) {
		if(VarI == ((n->rule).Var)) {
			if(((n->rule).OrdRule)<(*min)) (*min) = (n->rule).OrdRule;
			if(((n->rule).OrdRule)>(*max)) (*max) = (n->rule).OrdRule;
		}
		OrdFindMinMax(n->LeftC,VarI,min,max);
		OrdFindMinMax(n->RightC,VarI,min,max);
	}
}


//LOOP
void FindGoodCatRules(Node *n,int VarI, int *RuleInd,int &firstone)
// finds out which categorical rule using VarI are good.
// a good rule is one that does not result in logically empty bottom nodes
//n: the node at which the rule is to be set
//VarI: the variable ~ the rule
//RuleInd: integer vector whose length is the number of possible rules 2^(NR-1) - 1
//	on exit 1 if rule ok 0 otherwise, already allocated
//firstone: first category still "alive" at node n, depends on tree above n
{
  //Rprintf("FindGoodCatRules\n");

	int i,j;

	int NR = RuleNum(VarI);
	int *sel = new int [NR+1];
	int numrul =  (int)pow(2.0,NR-1)-1;
	for(i=1;i<=numrul;i++) RuleInd[i]=0;

	int *cats = new int [NR+1];
	GetSetCats(n,VarI,cats);
	firstone = FirstOne(NR,cats);
	if(!firstone) {
		Rprintf("error in FindGoodCatRule: no availble cats\n");
	} else {
		sel[firstone]=1;
	}

	int *sel1 = new int [NR-1+1];

	NodeP *lbotvec;
	int lnbot;
	MakeBotVec(n->LeftC,&lbotvec,&lnbot); //note lbotvec allocated here
	int *lfcount=new int [lnbot+1]; //allocation

	NodeP *rbotvec;
	int rnbot;
	MakeBotVec(n->RightC,&rbotvec,&rnbot); //note rbotvec allocated here
	int *rfcount=new int [rnbot+1]; //allocation

	
	// for eachrule see if it is any good
	// you do this be sending the categories left or right according to the rule
	// and checking to make sure that no bottom nodes are empty
	for(i=0;i<numrul;i++) {

		indtd(NR-1,i,sel1);
		
		for(j=1;j<firstone;j++) sel[j]=sel1[j];
		for(j=(firstone+1);j<=NR;j++) sel[j] = sel1[j-1];
		
		for(j=1;j<=lnbot;j++) lfcount[j]=0;
		for(j=1;j<=rnbot;j++) rfcount[j]=0;

		for(j=1;j<=NR;j++) {
			if(cats[j]) {
				if(sel[j]) {
					CatFindBots(n->RightC,VarI,j,rbotvec,rfcount);
				} else {
					CatFindBots(n->LeftC,VarI,j,lbotvec,lfcount);
				}
			}
			if((NoZero(lnbot,lfcount)) && (NoZero(rnbot,rfcount))) {
				RuleInd[i+1]=1;
				break;
			}
		}		
		
	}

	delete [] sel;
	delete [] sel1;
	delete [] cats;
	delete [] lbotvec;
	delete [] rbotvec;
	delete [] lfcount;
	delete [] rfcount;
}

//LOOP
int NoZero(int n,int *v)
{
  //Rprintf("NoZero\n");
	int retval = 1;
	int i;
	for(i=1;i<=n;i++) {
		if(v[i]==0) {
			retval=0;
			break;
		}
	}

	return retval;
}

//LOOP
int FirstOne(int n,int *v)
{
  //Rprintf("FirstOne\n");
	int i;
	for(i=1;i<=n;i++) {
		if(v[i]==1) return i;
	}

	return 0;
}

//LOOP
double ChangeRule(Node *top,int *Done)
// step which tries changing the rule 
{
  //Rprintf("ChangeRule\n");

	int i,j;
	double XLogPi,XLogL,YLogPi,YLogL;
	int ruleI;
	
	
	double alpha;
	double u;
	int Nnotbot;
	NodeP *notbotvec;
	
	// get list of nodes with rule = nodes which are not bottom
	MakeNotBotVec(top,&notbotvec,&Nnotbot);
	if(Nnotbot==0) {
		delete [] notbotvec;
		return -1;
	}
	
	// randomly choose a notbot node = cnode
	//u=ran1(&idum);
	u= unif_rand();
	int NodeI =  (int)floor(u*Nnotbot)+1;
	Node *cnode = notbotvec[NodeI];

	//given the node, choose a new variable for the new rule
	int YVarI = DrPriVar(cnode);

	// if new var is CAT do one thing, if ORD another
	if(VarType(YVarI)==CAT) {
		
		// get the list of good cat rules given var choice
		int firstone;
		int NR = RuleNum(YVarI);
		int numr = (int)pow(2.0,NR-1)-1;
		int *RuleInd = new int [numr+1];
		FindGoodCatRules(cnode,YVarI,RuleInd,firstone);
		int sum = 0;
		for(i=1;i<=numr;i++) sum += RuleInd[i];
		
		//if there are any good cat rules
		if(sum) {
			
			// draw the rule from list of good ones
			//u=ran1(&idum);
                        u = unif_rand();
			ruleI = (int)floor(u*sum)+1;
			ruleI = GetSkipBadInd(numr,RuleInd,ruleI);

			//get logpri and logL from current tree (X)
			XLogPi = LogPriT(top);
			XLogL = LogLT(cnode,top);

			// copy old rule
			Rule rule;
			CopyRule(&(cnode->rule),&rule);
			

			// change rule at cnode to the new one
			int *sel = new int [NR-1+1];
			indtd(NR-1,ruleI-1,sel);
			(cnode->rule).Var = YVarI;
			delete [] (cnode->rule).CatRule;
			(cnode->rule).CatRule = new int [NR+1];
			for(j=1;j<firstone;j++) (cnode->rule).CatRule[j]=sel[j];
			(cnode->rule).CatRule[firstone]=1;
			for(j=(firstone+1);j<=NR;j++) (cnode->rule).CatRule[j] = sel[j-1];
			
			//fix data at nodes below cnode given new rule
			FixDataBelow(cnode);
			
			//  fix VarAvail
			UpDateVarAvail(cnode,YVarI);
			if(!(YVarI==rule.Var)) UpDateVarAvail(cnode,rule.Var);
			
			

			//get logpri and logL from candidate tree (Y)
			YLogPi = LogPriT(top);
			YLogL = LogLT(cnode,top);
			
			//draw go nogo
			alpha = min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));
			if(Bern(alpha)) {
				
				
				*Done=1;
				


			} else {

				// if nogo put rule, data, and VarAvail back
				CopyRule(&rule,&(cnode->rule));
				FixDataBelow(cnode);

				//  fix VarAvail
				UpDateVarAvail(cnode,YVarI);
				if(!(YVarI==rule.Var)) UpDateVarAvail(cnode,rule.Var);
				
				*Done=0;
			}

			
			delete [] sel;

		
		} else {

			// if no rules for that var abort step
			alpha = -1;
		}

		delete [] RuleInd;


	} else {

		//ORD variable
		
		// get the set of good rules = [l,r]
		int l,r;
		FindGoodOrdRules(cnode,YVarI,l,r);
		int numsplit = r-l+1;

		// if there are any rules
		if(numsplit>0) {

			//draw the rule
			//u=ran1(&idum);
                        u = unif_rand();
			ruleI = l+(int)floor(u*numsplit);

			//get logpri and logL from current tree (X)
			XLogPi = LogPriT(top);
			XLogL = LogLT(cnode,top);
			
			// copy old rule
			int XVarI = (cnode->rule).Var;
			int XOrdRule = (cnode->rule).OrdRule;
			
			// change rule at cnode to the new one
			(cnode->rule).Var = YVarI;
			(cnode->rule).OrdRule = ruleI;
			
			//fix data at nodes below cnode given new rule
			FixDataBelow(cnode);

			UpDateVarAvail(cnode,YVarI);
			if(!(YVarI==XVarI)) UpDateVarAvail(cnode,XVarI);

			//get logpri and logL from candidate tree (Y)
			YLogPi = LogPriT(top);
			YLogL = LogLT(cnode,top);
			
			//draw go nogo
			alpha = min(1.0,exp(YLogPi+YLogL-XLogPi-XLogL));
			if(Bern(alpha)) {	
				// if go fix VarAvail
				*Done=1;
				
			} else {
				// if nogo put rule and data back
				(cnode->rule).Var = XVarI;
				(cnode->rule).OrdRule = XOrdRule;

				FixDataBelow(cnode);

				UpDateVarAvail(cnode,YVarI);
				if(!(YVarI==XVarI)) UpDateVarAvail(cnode,XVarI);

				*Done=0;
			}

		} else {
			// if no rules for that var abort step
			alpha=-1;
		}

	}
	delete [] notbotvec;
	return alpha; // note -1 means backed out
}




