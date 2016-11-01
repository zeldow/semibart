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
#include "ChangeRule.h"
#include "Swap.h"

typedef Node *NodeP;

//LOOP
void MakeBotVec(Node *top,NodeP **botvec,int *NBot)
// allocates and defines array of bot nodes for tree top, sets NBot
{
  //Rprintf("MakeBotVec\n");
	int i;

	List1 *bots;
	top->GetBotList(&bots);  //gets bottom nodes?
	*NBot = bots->length;    //number of bottom nodes?
	*botvec = new NodeP [*NBot+1]; //vector of bottom nodes?
	
	// rest puts each bottom node into slot of vector???	
	Cell *cell = bots->first; 
	(*botvec)[1]=(Node *)cell->contents;
	
	for(i=2;i<=(*NBot);i++) {
		cell = cell->after;
		(*botvec)[i]=(Node *)cell->contents;	
	}
	
	bots->deall();
	delete bots;
}

//LOOP
void MakeNogVec(Node *top,NodeP **nogvec,int *NNog)
// allocates and defines list of nog nodes for tree top, sets NNog
{
  //Rprintf("MakeNogVec\n");
	int i;

	List1 *bots;
	top->GetNogList(&bots);
	*NNog = bots->length;
	*nogvec = new NodeP [*NNog+1];
	
	if (*NNog){
	Cell *cell = bots->first;
	(*nogvec)[1]=(Node *)cell->contents;
	
	for(i=2;i<=(*NNog);i++) {
		cell = cell->after;
		(*nogvec)[i]=(Node *)cell->contents;	
	}
	}
	bots->deall();
	delete bots;
}

//LOOP
void MakeSwapVec(Node *top,NodeP **swapvec,int *Nswap)
// 
{
  //Rprintf("MakeSwapVec\n");
	int i;

	List1 *swaps;
	top->GetSwapsList(&swaps);
	*Nswap = swaps->length;
	*swapvec = new NodeP [*Nswap+1];
	
	if (*Nswap){
	Cell *cell = swaps->first;
	(*swapvec)[1]=(Node *)cell->contents;
	
	for(i=2;i<=(*Nswap);i++) {
		cell = cell->after;
		(*swapvec)[i]=(Node *)cell->contents;	
	}
	}
	swaps->deall();
	delete swaps;
}

//LOOP
void MakeNotBotVec(Node *top,NodeP **notbotvec,int *Nnotbot)
// allocates and defines list of nog nodes for tree top, sets NNog
{
  //Rprintf("MakeNotBotVec\n");
	int i;

	List1 *bots;
	top->GetNotBotList(&bots);
	*Nnotbot = bots->length;
	*notbotvec = new NodeP [*Nnotbot+1];
	
	if (*Nnotbot){
	Cell *cell = bots->first;
	(*notbotvec)[1]=(Node *)cell->contents;
	
	for(i=2;i<=(*Nnotbot);i++) {
		cell = cell->after;
		(*notbotvec)[i]=(Node *)cell->contents;	
	}
	}
	bots->deall();
	delete bots;
}



//LOOP
void MakeIntVec(List1 *intlist, ivec &ivec_, int *n)
//allocate and define vec of integers corresponding to list of int pointers
{
  //Rprintf("MakeIntVec\n");
	int i;

	*n = intlist->length;
	//*ivec = new int [*n + 1];
	ivec_.set_size(*n+1);

	Cell *cell = intlist->first;
	if(*n>0) ivec_(1)=*((int *)(cell->contents));
	
	for(i=2;i<=(*n);i++) {
		cell = cell->after;
		ivec_(i)=*((int *)(cell->contents));
	}
}


//Need to make sure divec is pointing to the right person

//LOOP
void AddDatChildren(Node *n)
{
  //Rprintf("AddDatChildren\n");
  if(!(n->rule).Var) 
		Rprintf("error in AddDatChildren: rule not set\n");
	if(((n->LeftC)->DataList.length!=0) || ((n->RightC)->DataList.length!=0))
		Rprintf("error in AddDatChildren: data already set\n");

	//int *divec;
	ivec divec;
	int nob;
	//MakeIntVec(&(n->DataList),&divec,&nob);
	MakeIntVec(&(n->DataList),divec,&nob);

	int i;
	for(i=1;i<=nob;i++) {
	  if((n->rule).Right(X.row(divec(i)-1))) {
	    (n->RightC)->SetData(divec(i));
			
	  }
	  else {
	    (n->LeftC)->SetData(divec(i));
	  }
	}
		
	
	//delete [] divec;
}

//LOOP
void FixDataBelow(Node *cnode)
{
  //Rprintf("FixDataBelow\n");
  //int *divec;
  ivec divec;
	int nobs;
	int i;
	
	(cnode->LeftC)->ClearData();
	(cnode->RightC)->ClearData();

	//MakeIntVec(&(cnode->DataList), &divec, &nobs);
	MakeIntVec(&(cnode->DataList), divec, &nobs);
	for(i=1;i<=nobs;i++) {
	  if ((cnode->rule).Right(X.row(divec(i)-1))) {
	    (cnode->RightC)->SetData(divec(i));
		} else {
	    (cnode->LeftC)->SetData(divec(i));
		}
	}
	//delete [] divec;

}



void UpDateOrdVarAvail(Node *n, int VarI, int left, int right)
{
  //Rprintf("updateordvaravail\n");
	int numsplit = right-left+1;
	if(numsplit<1) {
		(n->VarAvail)[VarI]=0;
	} else {
		(n->VarAvail)[VarI]=1;
	}

	if(!(n->Bot)) {
		int lleft,lright,rleft,rright;
		lleft=left;
		rleft=left;
		lright=right;
		rright=right;

		if(((n->rule).Var)==VarI) {
			lright = (n->rule).OrdRule - 1;
			rleft = (n->rule).OrdRule + 1;
		}

		UpDateOrdVarAvail(n->LeftC,VarI,lleft,lright);
		UpDateOrdVarAvail(n->RightC,VarI,rleft,rright);
	}
}

//LOOP
void UpDateCatVarAvail(Node *n, int VarI, int *cats)
{
  //Rprintf("updatecatvaravail\n");

	int i;
	int RN = RuleNum(VarI);

	if(ISum(RN,cats)<2) {
		(n->VarAvail)[VarI]=0;
	} else {
		(n->VarAvail)[VarI]=1;
	}

	if(!(n->Bot)) {
	
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

		UpDateCatVarAvail(n->LeftC,VarI,catsl);
		UpDateCatVarAvail(n->RightC,VarI,catsr);
	}

	delete [] cats;
}

//LOOP
void UpDateVarAvail(Node *n,int VarI)
{
  //Rprintf("updatevaravail\n");
        if(VarType(VarI)==CAT) {
	  int *cats = new int [RuleNum(VarI)+1];
		GetSetCats(n,VarI,cats);
		UpDateCatVarAvail(n,VarI,cats);// note cats is deleted in here
	} else {
		int LeftI,RightI;
		GetSplitInterval(&LeftI,&RightI, n,VarI);
		UpDateOrdVarAvail(n,VarI,LeftI,RightI);
	}


}


double Metrop(Node **top,int *Done,int *step)
{
  //Rprintf("Metrop\n");
	double alpha;
	int BD;

	//double u = unif_rand();
	double u = R::runif(0.0,1.0);
	//Rprintf("u: %.1f\n",u);
	if(u<pBD) {
		alpha = BirthDeath(*top,&BD,Done);
		if(BD) {
			*step = BIRTH;
		} else {
			*step = DEATH;
		}
	} else if(u<pBD+pSwap) { 
		alpha = SwapRule(*top,Done);
		*step=SWAP;
	} else {
		alpha = ChangeRule(*top,Done);
		*step = CHANGE;
	}

	return alpha;
}

