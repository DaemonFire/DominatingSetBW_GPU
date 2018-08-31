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
	if (argc!=2)
		return EXIT_FAILURE;
	graph* g = (graph*)malloc(sizeof(graph));
	*g =loadgraph(argv[1], 180);
	graph *h = (graph*)malloc(sizeof(graph));
	*h = *g;
	int* sol = (int*)malloc(g->size*sizeof(int));
	int size=0;
	//graph g = loadgraphformat2(argv[1]);

	printf("g.size=%d\n", g->size);
	printf("number of edges=%d\n", getEdgeNumber(*g));
	//dectree t = loadtree("tiefighter.tree");
	//dectree *t = generateTree(p,g,0);
	size = preprocessingsolopoints (g, sol, 180);
	int sizeinit=size;

	graph** components = (graph**)malloc(g->size*sizeof(graph*));
	int ncomp = computeconnexcomposants (g, components, 180);
	struct timeval stop, start;
	gettimeofday(&start, NULL);

	for (int i=0; i<ncomp; i++){
		dectree *t=generateTreeBW (*components[i]);
		int* x=(int*)malloc(g->size*sizeof(int));
		int sizetmp=0;
		if (components[i]->size>0)
			sizetmp = toplevelalgorithm (t, components[i], components[i]->size, x);
	
		if (sizetmp<0)
			sizetmp=0;

		for (int j=0; j<sizetmp; j++)
			sol[size+j]=x[j];
		size+=sizetmp;
		components[i]->sizeset=sizetmp;
		components[i]->domset=x;
	}
	for (int i=0; i<sizeinit; i++){
		printf("(%d, %d),  ", h->pos[2*sol[i]], h->pos[2*sol[i]+1]);
	}
	for (int i=0; i<ncomp; i++){
		for (int j=0; j<components[i]->sizeset; j++){
			printf("(%d, %d), ", components[i]->pos[2*sol[j]], g->pos[2*sol[j]+1]);
		}
	}
	gettimeofday(&stop, NULL);
	printf("\n");
	printf("There are %d composants\n", ncomp);

	printf("Minimum Dominating Set is of size %d\n",size);



	int timeToSet=stop.tv_usec - start.tv_usec;

	printf("Time elapsed for set=%d\n",timeToSet);


	return EXIT_SUCCESS;
}
