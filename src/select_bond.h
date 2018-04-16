/*
 * select_bond.h
 *
 *  Created on: 15-Mar-2018
 *      Author: Birva Patel
 */

#ifndef LMP_SELECT_BOND_H_
#define LMP_SELECT_BOND_H_


#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%" PRId64

#include <stdio.h>
#include <map>
#include <list>
#include<iostream>
#include"lammps.h"
#include "pointers.h"
using namespace std;
using namespace LAMMPS_NS;

class SelectBond: protected Pointers{
public:
	typedef int smallint;
	typedef int imageint;
	typedef int tagint;
	typedef int64_t bigint;

	int nlocal;
	int totalpairs;
	int *bondpair;
	int branchsize;
	int **bondAtom;	//(1 based)
	int *bondNum;	//(0 based)
	int **branchSizeMatrix;	//(1 based)
	int *branchNumber;
	int *cacheMatrix;		// Contains atoms count
	int *visited; //(1 based)
	int *visited_branch;
	int *visitedgroup;
	int *rotateGroupFlag; //(1 based)
	int **bondatompair;
	int atom1, atom2;
	int sizeGroup;
	int me;

	FILE *fp;

	list<int> *adj;     		// Pointer to an array containing
	list<int> *groupCache;  	// Groupcache

	SelectBond(LAMMPS *,char*);
	void makeBranchSizeMatrix();
	void addEdgeUndirected(int,int);
	void addEdgeDirected(int,int);
	void array2List();
	void cacheGroups();
	void cacheGroups(int,int);
	void selectBondPairAtoms();
	void selectDirection();
	void makeGroup();
	void createArrays();
	int branchSize(int);
	int size(int** );
	void getBondData(char*);
	int selectBondsAndDirectionsPlusCreateGroups();
	void addToGroup(int);
	void open(char*);
	void readline(char *);
	void parse_keyword(int, char *, char *);
	void skip_lines(int, char *);
	int parse(char *, char **, int);

	~SelectBond();
};

#endif /* LMP_SELECT_BOND_H_ */
