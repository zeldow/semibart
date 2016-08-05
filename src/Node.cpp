#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "MuS.h"

#include "global.h"
#include "List1.h"


#include "Rlob.h"
#include "Lib.h"

#include "Node.h"
#include "Funs.h"

#include "Prior.h"

#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

typedef double *dp;

double Rule::SplitVal()
{
  //Rprintf("Rule::SplitVal\n");
	if(Var==0) return -1000.0;
	if(VarType[Var]!=ORD) return -2000.0;
	return RuleMat[Var][OrdRule];
}


//LOOP
Rule::Rule(Rule &rule) {
  //Rprintf("Rule::Rule(rule)\n");
        Var=rule.Var; //variable of rule -- integer
	if(VarType(Var)==ORD) { 
		OrdRule = rule.OrdRule;
	} else {
	  CatRule = new int [RuleNum(Var)+1]; //trouble
		int i;
		for(i=1;i<=RuleNum(Var);i++) CatRule[i] = rule.CatRule[i];
	}
}

Rule::Rule()
{
  //Rprintf("Rule::Rule()\n");
	Var=0;
	OrdRule=0;
	CatRule=0;
	
}

Rule::~Rule()
{
  //Rprintf("Rule::~Rule\n");
	
  if(Var && (VarType(Var) == CAT)) {
		
		delete [] CatRule;
		
	
	}
}

void Rule::deall()
{
  //Rprintf("Rule::deall()\n");
  if(Var && (VarType(Var) == CAT)) {
		
		delete [] CatRule;
	
	}

	Var=0;
	OrdRule=0;
	CatRule=0;
	
}

//LOOP
int Rule::Right(rowvec x)
{
  int i;
  
  if(VarType(Var)==CAT) {
    
    for(i=1;i<=RuleNum(Var);i++) {
      
      if(x(Var-1) == RuleMat[Var][i]) {
	  
	if(CatRule[i]){
	  return 1;  
	} else {
	  return 0;
	}
	
      }
    }
    
    return 0;
  } else {
    if(x(Var-1)>RuleMat[Var][OrdRule]) {
      return 1;
    } else {
      return 0;
    }
  }
}
	

void CopyRule(Rule *r1,Rule *r2)
//rule of r1 copied to r2
{
  //Rprintf("CopyRule\n");
	int i;
	int RN;

	if(r2->Var) r2->deall(); //deletes rule of r2
	if(r1->Var) {  //if a Var is associated with this rule
	  r2->Var = r1->Var; //copy variable
	  if(VarType(r1->Var)==ORD) { //check if ordinal
			r2->OrdRule = r1->OrdRule;
	  } else {  //if categorical
	                //RN = RuleNum[r1->Var]; 
	                RN = RuleNum(r1->Var);
			r2->CatRule = new int[RN + 1];
			for(i=1;i<=RN;i++) (r2->CatRule)[i] = (r1->CatRule)[i];
		}
	}
}



//LOOP
void Node::deall() 
{
  //Rprintf("Node::deall()\n");
	if(!Bot) {
		LeftC->deall();
		RightC->deall();
		delete LeftC;
		delete RightC;
		rule.deall();
		Bot=1;
		Nog=0;
		if(!Top) {
			Node *brother = Brother(this);
			if(brother->Bot) (Parent->Nog)=1;
		}

	}
	if(Top) {
		Bot=1;
		Nog=0;
		if(DataList.length) DataList.deall();
		rule.deall();
		int i;
		for (i=1;i<=NumX;i++)
		  VarAvail[i]=1;
	}
 		
}



void Node::ClearData()
{
  //Rprintf("Node::ClearData\n");
	DataList.deall();
	if(!Bot) {
		LeftC->ClearData();
		RightC->ClearData();
	}
}

//LOOP
Node::Node()
{
  //Rprintf("Node::Node\n");
  Top = 1;
  	Bot = 1;
  	Nog = 0;

	VarAvail = new int [NumX+1];
	int i;
	for(i=1;i<=NumX;i++) VarAvail[i]=1;
}

Node::~Node()
{
  //Rprintf("Node::~Node\n");

	delete [] VarAvail;
	if(DataList.length) DataList.deall();

	
}


int Node::NumBotNodes()
// returns number of bottom nodes
{
  //Rprintf("Node::NumBotNodes\n");
	if(Bot) {
		return 1;
	} else {
		return (LeftC->NumBotNodes() + RightC->NumBotNodes()); 
	}
}

int Node::NumNogNodes()
{
  //Rprintf("Node::NumNogNodes\n");
	if(Bot) return 0;
	if(Nog) return 1;
	return (LeftC->NumNogNodes() + RightC->NumNogNodes());
}

void Node::GetBotList(List1 **list)
// gets list of pointers to bottom nodes
{
  //Rprintf("Node::GetBotList\n");
	if(Bot) {	
		*list=new List1;
		(**list).length=1;
		Cell *tempcell;
		tempcell=new Cell;
		(*tempcell).contents=this;
		(*tempcell).Beg=1;
		(*tempcell).End=1;
		(**list).first=tempcell;
		(**list).last=tempcell;
	} else {
		List1 *llist,*rlist;
		LeftC->GetBotList(&llist);
		RightC->GetBotList(&rlist);
		CombineLists(llist,rlist,list);
	}
}

void Node::GetNogList(List1 **list)
// gets list of pointers to nog nodes
{
  //Rprintf("Node::GetNogList\n");
	if(Bot) {
		*list=new List1;
		(**list).length=0;
	} else {

		if(Nog) {	
			*list=new List1;
			(**list).length=1;
			Cell *tempcell;
			tempcell=new Cell;
			(*tempcell).contents=this;
			(*tempcell).Beg=1;
			(*tempcell).End=1;
			(**list).first=tempcell;
			(**list).last=tempcell;
		} else {
			List1 *llist,*rlist;
			LeftC->GetNogList(&llist);
			RightC->GetNogList(&rlist);
			CombineLists(llist,rlist,list);
		}
	}
}



void Node::GetNotBotList(List1 **list)
// gets list of pointers to nog nodes
{
  //Rprintf("Node::GetNotBotList\n");

	if(Bot) {
		*list=new List1;
		(**list).length=0;
	} else if(Nog) {
			
		*list=new List1;
		(*list)->length=0;
		AddCellToEnd(*list,(void *)this);

	} else {
		List1 *llist,*rlist;
		LeftC->GetNotBotList(&llist);
		RightC->GetNotBotList(&rlist);
		CombineLists(llist,rlist,list);
		AddCellToEnd(*list,(void *)this);

	}
	
}

void Node::GetSwapsList(List1 **list)
{
  //Rprintf("Node::GetSwapsList\n");
	if(Nog || Bot) {
		*list=new List1;
		(**list).length=0;
	} else if((LeftC->Bot || LeftC->Nog) && (RightC->Bot || RightC->Nog)) {
		*list=new List1;
		(*list)->length=0;
		AddCellToEnd(*list,(void *)this);
	} else {

		List1 *llist,*rlist;
		LeftC->GetSwapsList(&llist);
		RightC->GetSwapsList(&rlist);
		CombineLists(llist,rlist,list);
		AddCellToEnd(*list,(void *)this);
	}
}




void Node::FindNode(vec x,Node **n)
// gets pointer to bottom node for x
{
  //Rprintf("Node::FindNode\n");
	
	if(Bot) {
		*n=this;
	}
	else {
		if(rule.Right(x)) {
			RightC->FindNode(x,n);
			
		}
		else {
			LeftC->FindNode(x,n);
		}
	}

}

//LOOP
voidP* Node::GetBotArray()
{
  //Rprintf("Node::GetBotArray\n");
	voidP* botArray=0;

	int nbot = NumBotNodes();

	botArray = new voidP[nbot+1];

	int i;

	List1 *bots;
	GetBotList(&bots);
		
	Cell *cell = bots->first;
	botArray[1] = (void*)cell->contents;
	
	for(i=2;i<=nbot;i++) {
		cell = cell->after;
		botArray[i] = (void*)cell->contents;
	}
	
	bots->deall();
	delete bots;

	return botArray;
}

//LOOP
ivec Node::GetIndPart(int numObs, mat xx)
{
  //Rprintf("Node::GetIndPart\n");
  ivec indPart(numObs+1);
  int i,j;
  voidP* botvec = GetBotArray();
  Node* nn;
  for(i=1;i<=numObs;i++)
    {
      FindNode(xx.row(i-1),&nn);
      for(j=1;((void*)nn) != botvec[j];j++);
      indPart(i) = j;
    }
  return indPart;
}

//LOOP
//changed output from double**
vec Node::GetFits(void* model,int  nTrain, mat xTrain, mat xTrainR, vec yTrain,  vec w)
{
  //Rprintf("Node::GetFits\n");
	int i,j;
	EndNodeModel* mod = (EndNodeModel*)model;
	vec fits(nTrain);

	ivec indPartTrain = GetIndPart(nTrain,xTrain);
	int nbot = NumBotNodes();
	
	ivec indObsTrain;
	
	int nobTrain=0;
	vec tempFits;
	int count;

	for(i=1;i<=nbot;i++)
	{
		nobTrain=0;
		for(j=1;j<=nTrain;j++) {if(indPartTrain(j)==i) nobTrain += 1;}

		indObsTrain.set_size(nobTrain+1);
		count=0;
		for(j=1;j<=nTrain;j++) {if(indPartTrain(j)==i) {count +=1; indObsTrain(count)=j;}}
		
		mod->setData(nobTrain,xTrainR,yTrain,indObsTrain,w);
		tempFits = mod->getFits(nobTrain,xTrainR,indObsTrain);	     
		for(j=1;j<=nobTrain;j++) fits(indObsTrain(j)-1) = tempFits(j);
		tempFits.zeros();
		tempFits.reset();
		indObsTrain.reset();
		

	}

	return fits;
}

//LOOP
void Node::currentFits(MuS* mod,int nTrain,mat xTrain,vec yTrain,vec w, vec &fits)
{
  //Rprintf("Node::currentFits\n");
        double ybar,postmu,postsd,b,a; //posterior of mu in a bottom node
        double nodeMu; //draw of mu, for a bottom node

        voidP* botvec = GetBotArray(); //bottom nodes
	int* indPartTest;
	
	int nbot = NumBotNodes();
	int nobTrain=0;
        int *itr;


	for(int i=1;i<=nbot;i++) { // loop over bottom nodes-------------
                //data is list of indices of train obs in the bottom node
	  //	  Rprintf("check1\n");
                List1& data = ((Node *)botvec[i])->DataList;
                nobTrain = data.length;
		//Rcout << nobTrain << std::endl;
                itr = new int[nobTrain+1]; //copy list contents to itr
		//	Rprintf("check2\n");
                Cell *cell = data.first;
                if(nobTrain>0) itr[1]=*((int *)(cell->contents));
		//Rcout << itr[1] << std::endl;
                ybar = yTrain(itr[1]-1);
                for(int j=2;j<=nobTrain;j++) {
                   cell = cell->after;
                   itr[j]=*((int *)(cell->contents));
                   ybar += yTrain(itr[j]-1);
                }
		//	Rprintf("check3\n");
                ybar /= nobTrain;

                b=nobTrain/mod->getSigma2();a=mod->getA();
                postmu = b*ybar/(a+b); postsd = 1.0/sqrt(a+b);
                nodeMu = postmu + postsd*norm_rand();
		//Rprintf("nodeMu: %.5f\n",nodeMu);
		
		for(int j=1;j<=nobTrain;j++) fits(itr[j]-1) = nodeMu; //altered
		//Rprintf("check4\n");
                delete [] itr;
	} //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        delete [] botvec;
}

//LOOP
// needs to be mindful of indices for mat
void Node::SetData()
{
  //Rprintf("Node::SetData()\n");
	int i;
	for(i=1;i<=NumObs;i++) SetData(i);
}

//LOOP
void Node::SetData(int i)
{
  //Rprintf("Node::SetData(%d)\n",i);
	
	Cell *cell;
	
	cell=new Cell;
	cell->contents= (&Ivec[i]); 
	cell->End=1;
	
	

	if(DataList.length==0) {
	
		DataList.first=cell;
		DataList.last=cell;
					
		cell->Beg=1;
		
			
		
	} else {
			
			
			DataList.last->End=0;
			DataList.last->after=cell;
			cell->before=DataList.last;
			DataList.last=cell;
			
			cell->Beg=0;
			
	}
	
	DataList.length +=1;
			
	if (!Bot){

	       if(rule.Right(X.row(i-1))) { //changed
			RightC->SetData(i);
			
		}
		else {
			LeftC->SetData(i);
		}
		
	
	}
	
	
}
	



// returns the depth of the node
int Depth(Node *nn)
{
  //Rprintf("Depth\n");
	int d=0;
	while(!(nn->Top)) {
		d += 1;
		nn = nn->Parent;
	}
	//Rprintf("\tDepth: %d\n",d);
	return d;
}


// returns the brother of the node (same depth)
Node *Brother(Node *n)
{
  //Rprintf("brother\n");
	if(n->Top) return 0;

	if(n==((n->Parent)->LeftC)) {
		return (n->Parent)->RightC;
	} else {
		return (n->Parent)->LeftC;
	}
}









