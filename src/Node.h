#include <RcppArmadillo.h>
#include <stdio.h>
#include <vector>
#include <cmath>

#include "List1.h"
#include "MuS.h"

using namespace Rcpp;
using namespace arma;

typedef void *voidP;

class Rule {
public:

	Rule();
	~Rule();
	Rule(Rule &rule);

	int Var; // index of variable associated with rule, 0 on creation
	int OrdRule; // if ordered this will be the index of the split value
	int *CatRule;
	//int Right(double *x); // returns 1 if vector x "goes" to right node, 0 else
	int Right(rowvec x);
	void deall();
	double SplitVal();
};

void CopyRule(Rule *r1,Rule *r2);

class Node {
public:

	Node();
	~Node();

	// next three are indicators for 
	// whether the node is a top, bottom, or nogrand
	int Top;
	int Bot;
	int Nog;
	
	// pointers for tree structure
	Node *Parent;
	Node *LeftC;
	Node *RightC;

	//the rule for splitting the observations
	Rule rule;
	int *VarAvail; // ith is 1 if variable i has rules left, 0 else
	
	//list of observations corresponding to node
	List1 DataList;

	//functions
	int NumBotNodes();  // returns number of bottom nodes
	int NumNogNodes(); // returns number of Nog nodes
	void GetBotList(List1 **list);  // gets list of pointers to bottom nodes
	void GetNogList(List1 **list);
	void GetNotBotList(List1 **list);
	void GetSwapsList(List1 **list);
	//rowvec???
	void FindNode(vec x,Node **n); // gets pointer to bottom node for x
	voidP* GetBotArray();
	//int* GetIndPart(int numObs, mat xx);  //maybe use later
	ivec GetIndPart(int numObs, mat xx);
	vec GetFits(void* model,int  nTrain,mat xTrain, mat xTrainR,vec yTrain, vec w);
	//double** GetEstimates(void* model,int  nTrain,double** dTrain, double** dTrainR,double* yTrain, double* w);
	//void  currentFits(MuS* mod,int nTrain,mat xTrain,vec yTrain, vec w, double** fits); //probably should switch in mat for fits
	void  currentFits(MuS* mod,int nTrain,mat xTrain,vec yTrain, vec w, vec &fits); 
       	void SetData(int i); //goes thru the tree and adds i to list of obs indices for each appropriate node
	void SetData();
	void ClearData(); // deallocate all the DataLists for all nodes in tree
	void deall();
	//void CopyTree(Node *copy);
	//int DepthBelow();
};

int Depth(Node *n);
Node *Brother(Node *n);




class VarUsage {
public:
   int depth;
   int nodeIndex;
   int varIndex;
};



