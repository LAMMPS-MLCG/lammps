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

	list<int> *adj;     		// Pointer to an array containing
	list<int> *groupCache;  	// Groupcache

	SelectBond(LAMMPS *);
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
	void getBondData();
	int selectBondsAndDirectionsPlusCreateGroups();
	void addToGroup(int);
	~SelectBond();
};

#endif /* LMP_SELECT_BOND_H_ */
