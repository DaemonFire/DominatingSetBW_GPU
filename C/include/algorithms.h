#ifndef _ALGORITHMS_
#define _ALGORITHMS_
#include "../include/datarepresentation.h"

__global__
void computeTwins (int* mat, int* width, int* result);

__global__
void computeRepresentatives(int* result, int* ptr);
	
__global__
void computeNeighboorhoods(int* lc, int* rc, int* lrc, int* lrcard, int* mat, int* width, int* na, int* rep, int* repc);

__global__
void computeAlgorithm (int* tabg, int* lra, int* lrb, int* lrw, int* lracard, int* lrbcard, int* lrwcard, int* lnracard, int* lnrbcard, int* lnrwcard, int* mw, int* macomp, int* mbcomp, int* nrepa, int* nrepb, int* nrepw, int* repacomp, int* repbcomp, int* repw, int* taba, int* tabb, int* ptrac, int* ptrbc, int* ptrw, int* nacomp, int* nbcomp, int* nw);

cutdata cutThatTree (graph* g, dectree* t);

__global__ 	
void computeMatrix(int* lr, int* ln, int* lrcard, int* tc, int* mg, int* mat, int* na, int* rep, int* repc, int* width, int* r);

int firstpreprocess(graph* g,  cutdata* c);

int secondpreprocess (cutdata* c, graph* g);

int thirdpreprocess (cutdata* c, graph* g);

int* toplevelalgorithm (dectree* t, graph* g, int n);

int stepalgorithm (dectree* t, graph* g);

int* computeDS (dectree* t, int much, int aleft, int acleft);

int getBW (dectree t, graph g);

#endif
