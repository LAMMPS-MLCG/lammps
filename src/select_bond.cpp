#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "select_bond.h"
#include <limits>
#include <list>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include "atom.h"
#include "lammps.h"
#include "group.h"
#include "memory.h"
#include "force.h"
#include "comm.h"
#include "error.h"
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MAXLINE 256
int INF = std::numeric_limits<int>::max();

using namespace LAMMPS_NS;
using namespace std;

/*TODO
 *
 * 6. implement group caching
 *
 */


SelectBond::SelectBond(LAMMPS *lmp, char* filename) : Pointers(lmp){
	sizeGroup = 0;
	nlocal = atom->natoms;
	atom1 = 0;
	atom2=0;
	bondpair = NULL;
	bondAtom = NULL;
	bondNum = NULL;
	branchNumber = NULL;
	rotateGroupFlag = NULL;
	branchSizeMatrix = NULL;
	visited = NULL;
	visited_branch = NULL;
	visitedgroup = NULL;
	cacheMatrix = NULL;
	bondatompair = NULL;
	fp = NULL;
	me = comm->me;
	srand(time(NULL));

	createArrays();
	getBondData(filename);
	array2List();
	//	for(int atom=1;atom<=nlocal;atom++){
	//		cout<<atom;
	//		for(auto i = adj[atom].begin(); i!= adj[atom].end(); i++){
	//			int v = *i;
	//			cout<<"->"<<v;
	//
	//		}
	//		cout<<endl;
	//	}
	//	for(int i=1;i<=nlocal;i++)
	//		cout<<i<<"-->"<<atom->map(i)<<endl;
	makeBranchSizeMatrix();
}

void SelectBond::makeBranchSizeMatrix() {
	int sum=0, max=0;
	memset(branchNumber,0,(nlocal+1)*sizeof(int));
	for(int i=1; i<=nlocal; i++){
		sum = 0; max = 0;
		memset(visited,0,(nlocal+1)*sizeof(int)); // 1 based

		for(std::list<int>::const_iterator iterator = adj[i].begin(); iterator != adj[i].end(); iterator++){
			int vertex = *iterator;
			memset(visited_branch,0,(nlocal+1)*sizeof(int)); // 1 based
			branchsize=0;
			visited[i]=1;
			visited_branch[i]=1;
			branchSizeMatrix[i][branchNumber[i]]=branchSize(vertex)+1;
			max = MAX(branchSizeMatrix[i][branchNumber[i]],max);
			branchNumber[i]++;
		}
		for(int j=1;j<=nlocal;j++)
			if(visited[j]) sum+=1;
		branchSizeMatrix[i][branchNumber[i]] = sum;
	}
}

int SelectBond::branchSize(int vertex) {
	if (visited_branch[vertex])
		return branchsize;

	visited_branch[vertex]=1;
	visited[vertex]=1;
	branchsize++;
	for(std::list<int>::const_iterator iterator = adj[vertex].begin(); iterator != adj[vertex].end(); ++iterator){
		int v = *iterator;
		if(!visited_branch[v]){
			branchSize(v);
		}
	}
	return branchsize;
}

//cache every Ni%10 vertex
void SelectBond::cacheGroups() {

	for(int x=0;x<=nlocal;x++)
		visited[x]=0;

	for(int i=0;i<=nlocal;i++)
		visitedgroup[i]=0;

	for(int i=10; i<nlocal; i += 10){
		for(int x=0;x<=nlocal;x++)
			visitedgroup[x]=0;
		cacheGroups(i/10,i);
	}
}

void SelectBond::cacheGroups(int index, int atom) {
	if(visitedgroup[atom])
		return;

	groupCache[index].push_back(atom);
	visitedgroup[atom] = 1;
	visited[atom] = 1;

	for(std::list<int>::const_iterator iterator = adj[atom].begin(); iterator != adj[atom].end(); ++iterator){
		cout<<*iterator<<'\n';
		int v = *iterator;
		if(visitedgroup[v])
			continue;
		else if(v%10)
			cacheGroups(index,v);
		else{
			groupCache[index].push_back(v);
			visitedgroup[v] = 1;
			visited[v] = 1;
		}
	}
}

void SelectBond::array2List() {
	for(int i=1 ; i<=nlocal ; i++){
		for(int j=0 ; j<bondNum[i] ; j++){
			if(bondAtom[i][j]!=0){
				addEdgeUndirected(i,bondAtom[i][j]);
			}
		}
	}
}

void SelectBond::addEdgeUndirected(int v, int w) {
	adj[v].push_back(w);
	adj[w].push_back(v);
}

void SelectBond::addEdgeDirected(int v,int w) {
	groupCache[v].push_back(w);
}

int SelectBond::size(int** arr) {
	return sizeof(arr) / sizeof(arr[0]);
}

void SelectBond::getBondData(char* filename) {
	tagint *tag = atom->tag;
	bondNum = new int[nlocal+1];
	memset(bondNum,0,(nlocal+1)*sizeof(int)); // 1 based

	for (int a1 = 0; a1 < nlocal; a1++) {
		for (int a2 = 0; a2 < atom->num_bond[a1]; a2++) {
			bondAtom[tag[a1]][a2] = atom->bond_atom[a1][a2];
			bondNum[tag[a1]]++;
		}
	}

	char line[MAXLINE],keyword[MAXLINE];
	int index = 0;
	open(filename);
	readline(line);
	parse_keyword(0,line,keyword);
	totalpairs = atoi (keyword);

	memory->create(bondatompair,totalpairs,2,"selectbond:pairs");

	for (int i = 0; i < totalpairs; i++){
		 readline(line);
		 sscanf (line,"%d\t%d",&bondatompair[index][0],&bondatompair[index][1]);
		 index++;
	}

	for(int i=0;i<totalpairs;i++)
		cout<<bondatompair[i][0]<<'\t'<<bondatompair[i][1]<<endl;
}

void SelectBond::selectBondPairAtoms() {

	int randindex = rand() % totalpairs;
	atom1 = bondatompair[randindex][0];
	atom2 = bondatompair[randindex][1];
}

void SelectBond::selectDirection() {
	//for atom1
	int countPositionatom1=0,countPositionatom2=0,branchsizeatom1,branchsizeatom2,totalsizeatom1,totalsizeatom2,netsizeatom1,netsizeatom2;
	for(std::list<int>::const_iterator iterator = adj[atom1].begin(); iterator != adj[atom1].end(); ++iterator){
		if((*iterator)== atom2)
			break;
		countPositionatom2++;
	}
	//for atom2
	for(std::list<int>::const_iterator iterator = adj[atom2].begin(); iterator != adj[atom2].end(); ++iterator){
		if((*iterator)== atom1)
			break;
		countPositionatom1++;
	}

	totalsizeatom1 = branchSizeMatrix[atom1][branchNumber[atom1]];
	totalsizeatom2 = branchSizeMatrix[atom2][branchNumber[atom2]];
	branchsizeatom1 = branchSizeMatrix[atom1][countPositionatom2];
	branchsizeatom2 = branchSizeMatrix[atom2][countPositionatom1];
	netsizeatom1 = totalsizeatom1 - branchsizeatom1;
	netsizeatom2 = totalsizeatom2 - branchsizeatom2;

	int swapAtom=0;
	if(netsizeatom2<netsizeatom1)
	{
		swapAtom = atom1;
		atom1 = atom2;
		atom2 = swapAtom;
	}
	bondpair[0] = atom->map(atom1);
	bondpair[1] = atom->map(atom2);

}

void SelectBond::makeGroup() {
	memset(rotateGroupFlag,0,(nlocal+1)*sizeof(int)); //exact size 1 based
	memset(visited,0,(nlocal+1)*sizeof(int)); // 1 based

	visited[atom2] = 1;
	rotateGroupFlag[atom1]=1;
	addToGroup(atom1);

	int s=0;
	sizeGroup = 0;
	for(int i=0;i<nlocal;i++)
		if(rotateGroupFlag[i])
			s++;
	sizeGroup = s;

	int igroup = group->find_or_create("rotategroup");
	if(igroup == -1)
		cout<<"cannot create the group"<<endl;
	else
	{
		int bit = group->bitmask[igroup];
		int inversebit = group->inversemask[igroup];
		int gid = group->find("rotategroup");
		if (gid > 0) {
			char *cmd[2];
			cmd[0] = (char *)"rotategroup";
			cmd[1] = (char *)"clear";
			group->assign(2,cmd);
		}

		igroup = group->find("rotategroup");
		bit = group->bitmask[igroup];
		int *flags = (int *)calloc(nlocal,sizeof(int));
		for (bigint i=1; i <= nlocal; ++i) {
			if(rotateGroupFlag[i])
			{
				const int id = atom->map(i);
				flags[id] = 1;
			}
		}
		group->create((char *)"rotategroup",flags);
		delete [] flags;
	}

}
//
//void SelectBond::addToGroup(int atom) {
//	if(visited[atom])
//		return;
//	visited[atom] = 1;
//	rotateGroupFlag[atom]=1;
//
//	for(std::list<int>::const_iterator iterator = adj[atom].begin(); iterator != adj[atom].end(); ++iterator){
//		int v = *iterator;
//		if(!visited[v]) //we keep this line mainly for processing optimization
//			addToGroup(v);
//	}
//}

void SelectBond::addToGroup(int atom) {
	list<int> queue;
	list<int>::iterator i;

	visited[atom] = 1;
	rotateGroupFlag[atom1]=1;
	queue.push_back(atom);

	while(!queue.empty()){
		atom = queue.front();
		queue.pop_front();

		for (i = adj[atom].begin(); i != adj[atom].end(); ++i) {
			if (!visited[*i]) {
				visited[*i] = 1;
				rotateGroupFlag[*i]=1;
				queue.push_back(*i);
			}
		}
	}
}

void SelectBond::createArrays() {
	memory->create(bondAtom,nlocal+1,atom->bond_per_atom,"selectbond:bondAtom");
	memory->create(bondpair,2,"selectbond:bondpair");
	memory->create(branchNumber,nlocal+1,"selectbond:branchNumber");
	memory->create(rotateGroupFlag,nlocal+1,"selectbond:rotateGroupFlag");
	memory->create(visited,nlocal+1,"selectbond:visited");
	memory->create(visited_branch,nlocal+1,"selectbond:visited_branch");
	memory->create(branchSizeMatrix,nlocal+1,atom->bond_per_atom,"selectbond:branchSizeMatrix");
	memory->create(cacheMatrix,nlocal+1,"selectbond:cacheMatrix");
	memory->create(visitedgroup,nlocal+1,"selectbond:visitedgroup");

	adj = new list<int>[nlocal+1]; //adjecency list of the bond atoms
}


int SelectBond::selectBondsAndDirectionsPlusCreateGroups()
{
	selectBondPairAtoms();
	selectDirection();
	makeGroup();
	return 0 ;
}

/* ----------------------------------------------------------------------
   open molecule file
------------------------------------------------------------------------- */

void SelectBond::open(char* file){
	fp = fopen(file,"r");
	if (fp == NULL) {
		char str[128];
		sprintf(str,"Cannot open bonds file %s",file);
		error->one(FLERR,str);
	}
}

/* ----------------------------------------------------------------------
   read and bcast a line
------------------------------------------------------------------------- */

void SelectBond::readline(char *line){
	int n;
	if (me == 0) {
		if (fgets(line,MAXLINE,fp) == NULL) n = 0;
		else n = strlen(line) + 1;
	}
	MPI_Bcast(&n,1,MPI_INT,0,world);
	if (n == 0) error->all(FLERR,"Unexpected end of bonds file");
	MPI_Bcast(line,n,MPI_CHAR,0,world);
}

/* ----------------------------------------------------------------------
   extract keyword from line
   flag = 0, read and bcast line
   flag = 1, line has already been read
------------------------------------------------------------------------- */

void SelectBond::parse_keyword(int flag, char *line, char *keyword){
	if (flag) {

		// read upto non-blank line plus 1 following line
		// eof is set to 1 if any read hits end-of-file

		int eof = 0;
		if (me == 0) {
			if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
			while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
				if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
			}
			if (fgets(keyword,MAXLINE,fp) == NULL) eof = 1;
		}

		// if eof, set keyword empty and return

		MPI_Bcast(&eof,1,MPI_INT,0,world);
		if (eof) {
			keyword[0] = '\0';
			return;
		}

		// bcast keyword line to all procs

		int n;
		if (me == 0) n = strlen(line) + 1;
		MPI_Bcast(&n,1,MPI_INT,0,world);
		MPI_Bcast(line,n,MPI_CHAR,0,world);
	}

	// copy non-whitespace portion of line into keyword

	int start = strspn(line," \t\n\r");
	int stop = strlen(line) - 1;
	while (line[stop] == ' ' || line[stop] == '\t'
			|| line[stop] == '\n' || line[stop] == '\r') stop--;
	line[stop+1] = '\0';
	strcpy(keyword,&line[start]);
}

SelectBond::~SelectBond() {
	memory->destroy(bondAtom);
	memory->destroy(bondpair);
	memory->destroy(branchNumber);
	memory->destroy(rotateGroupFlag);
	memory->destroy(visited);
	memory->destroy(visited_branch);
	memory->destroy(branchSizeMatrix);
	memory->destroy(cacheMatrix);
	memory->destroy(visitedgroup);
}

