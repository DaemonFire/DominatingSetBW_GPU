#ifndef _TREEPRIMITIVES_
#define _TREEPRIMITIVES_
#include "../include/datarepresentation.h"

int getallleaves(dectree t, int *list);

int getnumberofleaves (dectree t);

dectree *generateTree (pointset p, graph g, int verthor);

pointset getCandidates (pointset left, pointset right, graph g);

setwithinsets incrementun(graph g, pointset left, setwithinsets unleft, int i);

pointset incrementalUNheuristic (graph g, int init);

dectree *generateTreeBWstep (graph g, pointset dec, int i);

dectree *generateTreeBW (graph g, int z);

int printTree (dectree t);

#endif
