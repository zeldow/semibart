typedef Node *NodeP;


void GetSetCats(Node *curr,int VarI, int *cats);
void GetSplitInterval(int *LeftI, int *RightI, Node *curr,int VarI);


double PGrow(Node *n);
int SumGoodVar(Node *n);
//void DrawPrior(Node *n);
int DrPriVar(Node *n);
int GetSkipBadInd(int n,int *vec,int ind);
void DrPriRule(int VarI,Node *GNode,int *LeftEx,int *RightEx);
int SpawnChildren(Node *n,int LeftEx,int RightEx);

//int MaxDepth(Node *top);

double LogPriT(Node *n);

