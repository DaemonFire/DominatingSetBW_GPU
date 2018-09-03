#include "../include/generatedata.h"	

#include <time.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

int preprocessingsolopoints (graph* g, int* sol, int threshold){
	int nsolo=0;	
	for (int i=0; i<g->size; i++){
		int solo = 1;
		for (int j=0; j<g->size; j++){
			if (g->matrix[i*g->size+j]==1){
				solo=0;
				break;
			}
		}
		if (solo==1){
			sol[nsolo]=i;
			nsolo++;
		}
	}
	graph* h= (graph*)malloc(sizeof(graph));
	h->size=g->size-nsolo;
	h->matrix=(int*)malloc((g->size-nsolo)*(g->size-nsolo)*sizeof(int));
	h->pos=(int*)malloc(2*(g->size-nsolo)*sizeof(int));
	int cursor= 0;
	for (int i=0; i<g->size; i++){
		if (cursor<nsolo){
			if (i==sol[cursor]){
				cursor++;
			}
			else{
				h->pos[2*(i-cursor)]=g->pos[2*i];
				h->pos[2*(i-cursor)+1]=g->pos[2*i+1];
			}
		}
		else {
			h->pos[2*(i-cursor)]=g->pos[2*i];
			h->pos[2*(i-cursor)+1]=g->pos[2*i+1];
		}
	}
	generateEdges (*h, threshold);
	*g=*h;
	return nsolo;
}

int computeconnexcomposants (graph* g, graph** components, int threshold){
	int which = 0;
	int* computed = (int*)malloc(g->size*sizeof(int));
	for (int i=0; i<g-> size; i++)
		computed[i]=-1;
	int i=0;
	while (i<g->size){
		if (computed[i]==-1){
			computed[i]=which;
			which++;
		}
		for (int j=0; j<g->size; j++){
			if (g->matrix[i*g->size+j]==1){
				if (computed[j]==-1){
					computed[j]=computed[i];
				}
				else{
					if (computed[j]>=computed[i]){
						int previous=computed[j];
						computed[j]=computed[i];
						for (int k=0; k<g->size; k++){
							if (computed[k]==previous)
								computed[k]=computed[i];
						}
					}
					else{
						int previous = computed[i];
						int a = computed[j];
						for (int k=0; k<g->size; k++){
							if (computed[k]==previous)
								computed[k]=a;
						}
					}
				}

			}

		}
		i++;
	}
	for (int i=0; i<which; i++){
		int number=0;
		for (int j=0; j<g->size; j++){
			if (computed[j]==i){
				number++;
				break;
			}
		}
		if (number==0){
			which--;
			for (int j=0; j<g->size; j++){
				if (computed[j]>i)
					computed[j]--;
			}
			i--;
		}
	}

	for (int i=0; i<which; i++){
		components[i]=(graph*)malloc(sizeof(graph));
		components[i]->size=0;
		for (int j=0; j<g->size; j++){
			if (computed[j]==i)
				components[i]->size++;
		}
		components[i]->matrix=(int*)malloc(components[i]->size*components[i]->size*sizeof(int));
		components[i]->pos=(int*)malloc(2*components[i]->size*sizeof(int));
		int cursor=0;
		for (int j=0; j<g->size; j++){
			if (cursor>=components[i]->size)
				break;
			if (computed[j]==i){
				components[i]->pos[2*cursor]=g->pos[2*j];
				components[i]->pos[2*cursor+1]=g->pos[2*j+1];
				cursor++;
			}
		}
		generateEdges(*components[i], threshold);
	}
	return which;
}


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
		g.pos[2*i]=a;
		g.pos[2*i+1]=b;
	}
	fclose(f);
	generateEdges(g,threshold);
	return g;
}


void generateEdges (graph g, int threshold){

	for (int i = 0; i<g.size; i++) {
		for (int j= 0; j<i; j++) {
			if ((g.pos[2*i]-g.pos[2*j])*(g.pos[2*i]-g.pos[2*j])+(g.pos[2*i+1]-g.pos[2*j+1])*(g.pos[2*i+1]-g.pos[2*j+1])<=(threshold*threshold)){
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
