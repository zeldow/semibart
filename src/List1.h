#ifndef GUARD_List
#define GUARD_List

class Cell {
public:
	int Beg; //flag for is it the first cell(1=yes 0=no)
	int End; //flag for is it the last cell (1=yes 0=no)
	Cell *before;
	Cell *after;
	void *contents;
};

class List1 {
public:
	Cell *first;
	Cell *last;
	int length;
	void deall();
	List1();
	~List1();
};

typedef void *voidP;

void CombineLists(List1 *list1,List1 *list2, List1 **list); //used in Node.cpp
//void PrintList(List1 *list);
//void DelCell(List1 *list,Cell *cell);
void AddCellToEnd(List1 *list, void *p);  //used in Node.cpp
//void AddCellToBeg(List1 *list, void *p);
//void AddCellAfter(List1 *list,Cell *oldcell,void *p);
//void AddCellBefore(List1 *list,Cell *oldcell,void *p);
//void ListToVector(List1 *list,voidP **p,int *n);

#endif



