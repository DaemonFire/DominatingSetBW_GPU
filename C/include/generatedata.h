#ifndef _GENERATEDATA_
#define _GENERATEDATA_

#include "../include/datarepresentation.h"
#include <stdio.h>

int preprocessingsolopoints (graph* g, int* sol, int threshold);

int computeconnexcomposants (graph* g, graph** components, int threshold);

graph generategraph (int nbpoint, int xrange, int yrange, int threshold);

int getEdgeNumber(graph g);

int storegraph (graph g, char* name);

graph loadgraph (char* path, int threshold);

void generateEdges(graph g, int threshold);

dectree loadtree (char* path);

dectree lookTree (FILE *f, char *node);

int computeSizeFromFile (char* path);

int storetree (dectree t, FILE* f, char* node);

graph loadgraphformat2 (char* path);

#endif
