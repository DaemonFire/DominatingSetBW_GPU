#include "../include/generatedata.h"	

#include <time.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>


graph generategraph (int nbpoint, int xrange, int yrange, int threshold){

	srand (time(NULL));
	graph g;
	g.size = nbpoint;
	g.matrix = (int*)malloc(nbpoint*nbpoint*sizeof(int));
	g.pos = (int*)malloc(nbpoint*3*sizeof(int));

	for (int i = 0; i < nbpoint; i++){
		int x = rand() % xrange;
		int y = rand() % yrange;
		g.pos[3*i] = i;
		g.pos[3*i+1] = x;
		g.pos[3*i+2] = y;		 
	}
	generateEdges(g,threshold);

	return g;
}

int getEdgeNumber(graph g){
	int n=0;
	for (int i=0;i<g.size;i++){
		for (int j=i+1;j<g.size;j++){
			if (g.matrix[i*g.size+j]==1)
			n++;
		}
	}
	return n;
}
int storegraph (graph g, char* name){
	int fd;
	if(fd = open (name, O_CREAT|O_TRUNC|O_RDWR)==-1){
		fprintf(stderr,"could not open the file\n");
		return EXIT_FAILURE;
	}

	chmod (name, 777);
	char* str = (char*) malloc (100);

	for (int i=0; i<g.size; i++){
		sprintf(str, "%d %d\n", g.pos[3*i+1], g.pos[3*i+2]);
		write (fd, str, strlen(str));
	}
	return fd;
}


graph loadgraph (char* path, int threshold) {
	FILE *f;
	char buffer[21];
	char* subtoken;
	graph g;
	g.size = computeSizeFromFile(path);
	if ((f=fopen(path,"r")) == NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}
	g.pos = (int*)malloc(3*g.size*sizeof(int));
	g.matrix = (int*)malloc(g.size*g.size*sizeof(int));

	for (int i=0; i<g.size; i++){
		fgets(buffer,80,f);
		int a = atoi(subtoken=strtok(buffer," "));
		int b = atoi(subtoken=strtok(NULL," "));
		g.pos[3*i]=i;
		g.pos[3*i+1]=a;
		g.pos[3*i+2]=b;
	}
	fclose(f);
	generateEdges(g,threshold);
	return g;
}


void generateEdges (graph g, int threshold){

	for (int i = 0; i<g.size; i++) {
		for (int j= 0; j<i; j++) {
			if ((g.pos[3*i+1]-g.pos[3*j+1])*(g.pos[3*i+1]-g.pos[3*j+1])+(g.pos[3*i+2]-g.pos[3*j+2])*(g.pos[3*i+2]-g.pos[3*j+2])<=(threshold*threshold)){
				g.matrix[g.size*i+j]=1;
				g.matrix[g.size*j+i]=1;
			}
			else {
				g.matrix[g.size*i+j]=0;
				g.matrix[g.size*j+i]=0;
			}
		}
		g.matrix[g.size*i+i]=0;
	}
}


dectree loadtree (char* path) {
	FILE *f;
	char buffer[21];
	char* subtoken;
	dectree t;

	if ((f=fopen(path,"r")) == NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}
	(fgets(buffer,80,f));
	char* a = strtok(buffer," ");
	int b = atoi(subtoken=strtok(NULL," "));
	if (b==-1){
		char* c = strtok(NULL," ");
		char* e = strtok(NULL," ");
		char* d = (char*)malloc(strlen(e)-1);
		strncpy(d,e,strlen(e)-1);
		dectree q = lookTree(f, c);
		dectree r = lookTree(f,d);
		t.label=-1;
		t.left=(dectree*)malloc(sizeof(dectree));
		t.right=(dectree*)malloc(sizeof(dectree));
		*(t.left) = q;
		*(t.right) = r;
	}

	else {
		t.label=b;
		t.right=NULL;
		t.left=NULL;
	}

	fclose(f);
	return t;
}


dectree lookTree (FILE *f, char *node){
	int found=0;
	dectree t;
	char buffer[81];
	char* subtoken;
	fseek(f,0,SEEK_SET);

	while ((found==0)&&(fgets(buffer,80,f)!=NULL)){
		char* a = strtok(buffer," ");
		if (strcmp(a,node)==0){
			int b = atoi(subtoken=strtok(NULL," "));
			if (b==-1){
				char* c = strtok(NULL," ");
				char* e = strtok(NULL," ");
				char* d = (char*)malloc(strlen(e)-1);
				strncpy(d,e,strlen(e)-1);
				dectree q = lookTree(f, c);
				dectree r = lookTree(f,d);
				t.label=-1;
				t.left=(dectree*)malloc(sizeof(dectree));
				t.right=(dectree*)malloc(sizeof(dectree));
				*(t.left) = q;
				*(t.right) = r;
			}
			else {
				t.label=b;
				t.right=NULL;
				t.left=NULL;
			}
			found=1;
		}
	}

	return t;
}


int computeSizeFromFile(char* path){
	FILE *f;
	char buffer[80];

	if ((f=fopen(path,"r"))==NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}

	int i = 0;

	while (fgets(buffer,80,f)!=NULL){
		i++;
	}

	fclose(f);
	return i;
}


int storetree (dectree t, FILE* f, char* node){

	char* str = (char*) malloc (100);
	if (t.label!=-1){
		sprintf(str, "%s %d\n", node, t.label);		
		fwrite ( str, strlen(str), 1, f);	
	}
	else {
		char* tmp1 = (char*)malloc(strlen(node)+1);
		char* tmp2 = (char*)malloc(strlen(node)+1);
		for (int i=0;i<strlen(node);i++){
			tmp1[i]=node[i];
			tmp2[i]=node[i];
		}
		tmp1[strlen(node)]='1';
		tmp2[strlen(node)]='2';
		sprintf(str, "%s -1 %s %s\n", node, tmp1, tmp2);
		fwrite (str, strlen(str), 1, f);
		storetree(*(t.left), f, tmp1);
		storetree(*(t.right), f, tmp2);
	}

	return EXIT_SUCCESS;
}


graph loadgraphformat2 (char* path){
	FILE *f;
	char buffer[21];
	char* subtoken;
	graph g;

	if ((f=fopen(path,"r")) == NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}

	fgets(buffer,80,f);
	subtoken=strtok(buffer," ");
	subtoken=strtok(NULL," ");
	g.size = atoi(subtoken=strtok(NULL," "));
	int numberEdges = atoi(subtoken=strtok(NULL," "));
	g.matrix = (int*)malloc(g.size*g.size*sizeof(int));
	for (int i=0; i<g.size*g.size;i++)
		g.matrix[i]=0;
	for (int i=0; i<numberEdges; i++){
		fgets(buffer,80,f);
		subtoken=strtok(buffer," ");
		int a = atoi(subtoken=strtok(NULL," "))-1;
		int b = atoi(subtoken=strtok(NULL," "))-1;
		g.matrix[a*g.size+b]=1;
		g.matrix[b*g.size+a]=1;
	}

	fclose(f);

	return g;
}
