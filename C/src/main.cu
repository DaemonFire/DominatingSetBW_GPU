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
	graph g = loadgraph(argv[1], 180);

	//graph g = loadgraphformat2(argv[1]);

	printf("g.size=%d\n", g.size);
	printf("number of edges=%d\n", getEdgeNumber(g));
	//dectree t = loadtree("tiefighter.tree");
	//dectree *t = generateTree(p,g,0);

	struct timeval stop, start;
	gettimeofday(&start, NULL);

	dectree *t=generateTreeBW (g);
	gettimeofday(&stop, NULL);
	int timeToTree=stop.tv_usec - start.tv_usec;
	printf("Time elapsed for tree=%d\n",timeToTree);
/*
	FILE *f;
	if ((f=fopen("arbre.tree","r"))==NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}


	storetree (*t, f, "0");*/
/*	int bwmax = getBW(*t,g);
	printf("Boolean-Width=%d\n",bwmax);*/

//	printTree (*t);
	graph* l;
	l=(graph*)malloc(sizeof(graph));
	*l=g;
	int* x = toplevelalgorithm (t, l);
	gettimeofday(&stop, NULL);


	//printf("Minimum Dominating Set is of size %d\n",x.size);
/*	for (int i=0;i<x.size;i++){
		printf("%d %d\n",g.pos[3*x.members[i]+1], g.pos[3*x.members[i]+2]);
	}*/
	//printf("\n");

	//generatePlotFile (*t, g);


	int timeToSet=stop.tv_usec - start.tv_usec;

	printf("Time elapsed for set=%d\n",timeToSet);


	return EXIT_SUCCESS;
}
