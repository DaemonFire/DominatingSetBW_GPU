#ifndef _DATAREPRESENTATION_
#define _DATAREPRESENTATION_
#include <stddef.h>

typedef struct pointset {
	int size;
	int* members;
} pointset;

typedef struct setwithinsets {
	pointset *set;
	int size;
} setwithinsets;

typedef struct cutdata {
	int *matrixrevisited;
	int na;
	int *a;
	int nacomp;
	int *acomp;

	int *tc;
	int nrep;
	int *complementtc;
	int nrepincomp;
	int *pointtorep;
	int *pointtorepincomp;

	int *lra;
	int lracard;
	int *lnra;
	int lnracard;

	int *lracomp;
	int lracompcard;
	int *lnracomp;
	int lnracompcard;

	int *m;
	int *mcomp;

	int* tab;

	int* box;
} cutdata;

typedef struct dectree {
	struct dectree *left;
	struct dectree *right;
	int label;
	cutdata c;
	int computed;
} dectree; 

typedef struct graph {
	int* matrix;
	int* pos;
	int size;
	int sizeset;
	int* domset;
} graph;

#endif


