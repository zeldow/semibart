#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

typedef void *voidP;
/**
extern "C" {
#include <R.h>
#include <Rmath.h>
};*/
#include <RcppArmadillo.h>
#include <Rmath.h>

/**
List1 contains data for each Node
List1 also used as for lists of Node (e.g. lists of all bottom nodes on tree
Cells are the individual elements of each List1
*/

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

List1::~List1()  //destructor?
{
  //deall();
}

List1::List1() //constructor
{
	length=0;
}

//LOOP
void List1::deall()
{
	if(length>0) {
		Cell *kill=0,*next=0;
		kill=first;
		for(int i=1;i<=length;i++) {
			if(i<length) next=kill->after;
			delete kill;
			kill=next;
		}
		length=0;
	}
	
}


void CombineLists(List1 *list1,List1 *list2, List1 **list)
{
	int n1,n2;

	n1=(*list1).length; //dereferences List and finds length
	n2=(*list2).length;
	
	if(n1==0) {
		(*list)=list2;
		delete list1;
		return;
	}
	if(n2==0) {
		(*list)=list1;
		delete list2;
		return;
	}
	if((n1>0) && (n2>0)) {			
	
		(*list)=new List1;
	
	
		(**list).length=n1+n2;
		(**list).first=(*list1).first;
		(**list).last=(*list2).last;
	
		((*list1).last)->after=(*list2).first;
		((*list1).last)->End=0;
		((*list2).first)->before=(*list1).last;
		((*list2).first)->Beg=0;
	
		delete list1;
		delete list2;
	}
	
}

//adds information p to end of list as Cell??
void AddCellToEnd(List1 *list, void *p)
{
	int len = list->length;

	Cell *cell;
	cell = new Cell;
	cell->contents = p;
	cell->End = 1;
	if(len) {
		(list->last)->End=0;
		(list->last)->after = cell;

		cell->before = list->last;
		cell->Beg=0;
	} else {
		list->first=cell;
		cell->Beg=1;
	}
	list->last = cell;
	list->length +=1;
}

















