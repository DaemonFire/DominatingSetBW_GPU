#include "../include/datarepresentation.h"
#include "../include/generatedata.h"
#include "../include/treeprimitives.h"
#include "../include/algorithms.h"

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <sys/time.h>
#include <time.h>

int main (int argc, char** argv){
	if (argc!=3)
		return EXIT_FAILURE;
	graph* g = (graph*)malloc(sizeof(graph));
	int threshold = atoi(argv[2]);
	*g =loadgraph(argv[1], threshold);
	graph *h = (graph*)malloc(sizeof(graph));
	*h = *g;
	int* sol = (int*)malloc(g->size*sizeof(int));
	int size=0;
	size = preprocessingsolopoints (g, sol, threshold);
	int sizeinit=size;

	printf("g.size=%d\n", g->size);
	printf("number of edges=%d\n", getEdgeNumber(*g));
	//dectree t = loadtree("tiefighter.tree");
	//dectree *t = generateTree(p,g,0);
	graph** components = (graph**)malloc(g->size*sizeof(graph*));
	int ncomp = computeconnexcomposants (g, components, threshold);
	inital(ncomp);
	struct timeval stop, start;
	gettimeofday(&start, NULL);
	for (int i=0; i<ncomp; i++){
	//	for (int j=0; j< components[i]->size; j++)
	//		printf("(%d, %d)\n",components[i]->pos[2*j], components[i]->pos[2*j+1]);
		dectree *t=generateTreeBW (*components[i], i);
		int* x=(int*)malloc(components[i]->size*sizeof(int));
		int nx=0;
		if (components[i]->size>0){
			nx = toplevelalgorithm (t, components[i], components[i]->size, x, i);
		}	
		else
			nx=0;
	
		if (nx<0)
			nx=0;
		
		for (int j=0; j<nx; j++)
			sol[size+j]=x[j];
	/*	for (int j=0; j<nx; j++)
			printf("%d, ", x[j]);
		printf("\n");*/
		size+=nx;
		components[i]->sizeset=nx;
		components[i]->domset=x;
	}
	for (int i=0; i<sizeinit; i++){
		printf("(%d, %d),  ", h->pos[2*sol[i]], h->pos[2*sol[i]+1]);
	}
	printf("|");
	int cursor =sizeinit;
	for (int i=0; i<ncomp; i++){
		for (int j=0; j<components[i]->sizeset; j++){
			printf("(%d, %d), ", components[i]->pos[2*sol[cursor+j]], components[i]->pos[2*sol[cursor+j]+1]);
		}
		cursor+=components[i]->sizeset;
		printf("|");
	}
	printf("\n");
	printf("There are %d composants\n", ncomp);
	gettimeofday(&stop, NULL);
	printf("Minimum Dominating Set is of size %d\n",size);
	//printf("\n");
	//generatePlotFile (*t, g);


	long timeToSet=1000000*(stop.tv_sec-start.tv_sec)+stop.tv_usec - start.tv_usec;

	printf("Time elapsed for set=%ld\n",timeToSet);

	return EXIT_SUCCESS;
}
