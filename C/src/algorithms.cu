#include "../include/algorithms.h"
#include "../include/treeprimitives.h"


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <errno.h>
#include <sys/stat.h>
#include <pthread.h>

int nnodes;
dectree** nodestocompute;
int tocompute;
pthread_mutex_t mut;
graph* gwork;

int fillThevoid (dectree* t, graph* g){
	gwork=g;
	t->computed=0;
	nodestocompute[tocompute]=t;
	tocompute++;
	if (t->left!=NULL)
		fillThevoid(t->left, g);
	if (t->right!=NULL)
		fillThevoid(t->right, g);
	return EXIT_SUCCESS;
}

int getnumberofnodes(dectree* t){
	int n=1;
	if (t->left!=NULL)
		n+=getnumberofnodes(t->left);
	if (t->right!=NULL)
		n+=getnumberofnodes(t->right);
	return n;
}

__global__
void computeTwins (int* mat, int* width, int* result){
	int twin=1;
	for (int i=0; i<*width; i++){
		if (mat[blockIdx.x*(*width) + i]!=mat[threadIdx.x*(*width)+i])
			twin = 0;
	}
	if (twin==1)
		result[threadIdx.x*blockDim.x+blockIdx.x]=1;
}

__global__
void computeRepresentatives(int* result, int* ptr){
	for (int i=0; i<blockDim.x; i++){
		if (result[threadIdx.x*blockDim.x+i]==1){
			ptr[threadIdx.x*2]=threadIdx.x;
			if (i<threadIdx.x)
				ptr[threadIdx.x*2+1]=i;
			else
				ptr[threadIdx.x*2+1]=threadIdx.x;
			break;
		}
	}
}
	
__global__
void computeNeighboorhoods(int* lc, int* rc, int* lrc, int* lrcard, int* mat, int* width, int* na, int* nacomp, int*a, int* acomp, int* rep, int* repc){
	for (int i = 0; i<*width; i++){
		int indexi=-1;
		for (int k=0; k<(*nacomp); k++){
			if (acomp[k]==repc[i]){
				indexi=k;
				break;	
			}
		}
		int isin = 0;
		for (int j = 0 ; j<blockDim.x; j++){
			int indexj=-1;
			for (int k=0; k<(*na); k++){
				if (a[k]==rep[j]){
					indexj=k;
					break;
				}
			}
			if ((lc[(blockIdx.x*blockDim.x+threadIdx.x)*blockDim.x+j]==1)&&(mat[indexj*(*nacomp)+indexi]==1)){
				isin=1;
				break;
			}
		}
		if (isin==1)
			rc[(blockIdx.x*blockDim.x+threadIdx.x)*(*width)+i]=1;
		else
			rc[(blockIdx.x*blockDim.x+threadIdx.x)*(*width)+i]=0;
	}
	/*
	int already=0;
	for (int i=0; i<*lrcard; i++){
		int id=1;
		for (int j = 0; j<*width;j++){
			if (rc[(blockIdx.x*blockDim.x+threadIdx.x)*(*width)+j]!=lrc[i*(*width)+j]){
				id=0;
				break;
			}
		}
		if (id == 1){
			already=1;
			break;
		}
	}
	if (already==1)
		rc[(blockIdx.x*blockDim.x+threadIdx.x)*(*width)]=-1;
	*/
}

__global__
void computeAlgorithm (int* tabg, int* lra, int* lrb, int* lrw, int* lracard, int* lrbcard, int* lrwcard, int* lnracard, int* lnrbcard, int* lnrwcard, int* mw, int* macomp, int* mbcomp, int* nrepa, int* nrepb, int* nrepw, int* repacomp, int* repbcomp, int* repw, int* taba, int* tabb, int* ptrac, int* ptrbc, int* ptrw, int* nacomp, int* nbcomp, int* nw){
	int indexa = threadIdx.x/(*lrbcard);
	int indexb = threadIdx.x%(*lrbcard);
	int indexbc = 0;
	int indexac = 0;
	int indexw = 0;

	for (int i = 0; i<(*nrepw); i++){
		if (lrw[blockIdx.x*(*nrepw)+i]==1){
			int rep = -1;
			for (int j=0; j< (*nbcomp); j++){
				if (ptrbc[2*j]==repw[i]){
					rep=ptrbc[2*j+1];
					break;
				}
			}
			for (int j = 0; j< (*nrepb); j++){
				if (repbcomp[j]==rep){
					rep=j;
					break;
				}
			}
			indexbc = mbcomp[indexbc*(*nrepb)+i];
		}
	}
	for (int i = 0; i<(*nrepa); i++){
		if (lra[indexa*(*nrepa)+i]==1){
			int rep = -1;
			for (int j=0; j< (*nbcomp); j++){
				if (ptrbc[2*j]==repacomp[i]){
					rep=ptrbc[2*j+1];
					break;
				}
			}
			for (int j = 0; j< (*nrepb); j++){
				if (repbcomp[j]==rep){
					rep=j;
					break;
				}
			}
			indexbc = mbcomp[indexbc*(*nrepb)+i];
		}
	}
	for (int i = 0; i<(*nrepw); i++){
		if (lrw[blockIdx.x*(*nrepw)+i]==1){
			int rep = -1;
			for (int j=0; j< (*nacomp); j++){
				if (ptrac[2*j]==repw[i]){
					rep=ptrac[2*j+1];
					break;
				}
			}
			for (int j = 0; j< (*nrepa); j++){
				if (repacomp[j]==rep){
					rep=j;
					break;
				}
			}
			indexac = macomp[indexac*(*nrepa)+i];
		}
	}
	for (int i = 0; i<(*nrepb); i++){
		if (lrb[indexb*(*nrepb)+i]==1){
			int rep = -1;
			for (int j=0; j< (*nacomp); j++){
				if (ptrac[2*j]==repbcomp[i]){
					rep=ptrac[2*j+1];
					break;
				}
			}
			for (int j = 0; j< (*nrepa); j++){
				if (repacomp[j]==rep){
					rep=j;
					break;
				}
			}
			indexac = macomp[indexac*(*nrepa)+i];
		}
	}
	for (int i = 0; i<(*nrepa); i++){
		if (lra[indexa*(*nrepa)+i]==1){
			int rep = -1;
			for (int j=0; j< (*nw); j++){
				if (ptrw[2*j]==repacomp[i]){
					rep=ptrw[2*j+1];
					break;
				}
			}
			for (int j = 0; j< (*nrepw); j++){
				if (repw[j]==rep){
					rep=j;
					break;
				}
			}
			indexw = mw[indexw*(*nrepw)+i];
		}
	}
	for (int i = 0; i<(*nrepb); i++){
		if (lrb[indexb*(*nrepb)+i]==1){
			int rep = -1;
			for (int j=0; j< (*nw); j++){
				if (ptrw[2*j]==repbcomp[i]){
					rep=ptrw[2*j+1];
					break;
				}
			}
			for (int j = 0; j< (*nrepw); j++){
				if (repw[j]==rep){
					rep=j;
					break;
				}
			}
			indexw = mw[indexw*(*nrepw)+i];
		}
	}

	tabg[5*(*lra)*(*lrb)*blockIdx.x+5*(*lrb)*indexa+5*indexb]=indexw;
	tabg[5*(*lra)*(*lrb)*blockIdx.x+5*(*lrb)*indexa+5*indexb+1]=indexac;
	tabg[5*(*lra)*(*lrb)*blockIdx.x+5*(*lrb)*indexa+5*indexb+2]=indexbc;
	tabg[5*(*lra)*(*lrb)*blockIdx.x+5*(*lrb)*indexa+5*indexb+3]=taba[indexa*(*lnracard)+indexac];
	tabg[5*(*lra)*(*lrb)*blockIdx.x+5*(*lrb)*indexa+5*indexb+3]=tabb[indexb*(*lnrbcard)+indexbc];
}

cutdata cutThatTree (graph* g, dectree* t){
	cutdata c;
	c.na=0;
	c.nacomp=0;
	c.a=NULL;
	c.acomp=NULL;

	c.na=getnumberofleaves (*(t));
	c.a=(int*)malloc(c.na*sizeof(int));
	getallleaves(*(t), c.a);
	
	c.nacomp=g->size-c.na;

	c.acomp=(int*)malloc(c.nacomp*sizeof(int));

	int i=0;
	int j=0;
	int k=0;
	for (i=0;i<g->size;i++){
		int ina=0;
		for (j=0;j<c.na;j++){
			if (c.a[j]==i){
				ina=1;
				break;
			}
		}
		if (ina==0){
			c.acomp[k]=i;
			k++;
			if (k==c.nacomp)
				break;
		}
	}

	c.matrixrevisited = (int*)malloc(c.na*c.nacomp*sizeof(int));


	for (int i=0;i<c.na;i++){
		for (int j=0;j<c.nacomp;j++)
			c.matrixrevisited[i*c.nacomp+j]=g->matrix[c.a[i]*g->size+c.acomp[j]];
	}
	/*	printf("M=\n");
		for (int i =0; i<c.na; i++){
			for (int j= 0; j<c.nacomp; j++){
				printf("%d ", c.matrixrevisited[i*c.nacomp+j]);
			}
			printf("\n");
		}*/
	return c;
}


__global__ 	
void computeMatrix(int* lr, int* ln, int* lrcard, int* tc, int* mg, int* mat, int* na, int* nacomp, int* a, int* acomp, int* rep, int* repc, int* width, int* r){
	int index= -1;
	for (int i=0; i<(*na);i++){
		if (a[i]==rep[threadIdx.x]){
			index=i;
			break;
		}
	}
	for (int i = 0; i <*width; i++){
		int indexi=-1;
		for (int j=0; j<(*nacomp);j++){
			if (acomp[j]==repc[i]){
				indexi=j;
				break;
			}
		}
		int isin=0;
		if (mat[index*(*nacomp)+indexi]==1)
			isin = 1;
		if (isin==0){
			for (int j = 0; j<blockDim.x; j++){
				int indexj=-1;
				for (int k=0; k<(*na);k++){
					if (a[k]==rep[j]){
						indexj=k;
						break;
					}
				}
				if ((lr[blockIdx.x*blockDim.x+j]==1)&&(mat[indexj*(*nacomp)+indexi]==1)){
					isin = 1;
					break;
				}
			}
		}		
		if (isin==1)
			r[i]=1;
		else
			r[i]=0;
	}
	int answer=-1;
	for (int i = 0; i <(*lrcard); i++){
		int id = 1;
		for (int j = 0; j<(*width); j++){
			if (r[j]!=lr[i*(*width)+j]){
				id = 0;
				break;
			}
		}
		if (id==1){
			answer=i;
			break;
		}
	}
	mg[blockIdx.x*blockDim.x+threadIdx.x]=answer;
}

int firstpreprocess(graph* g,  cutdata* c){

	int* res;

	res= (int*) malloc(c->na*c->na*sizeof(int));

	for (int i=0; i<c->na*c->na; i++)
		res[i]=0;

	int* mat;
	int* width;
	int * result;

	cudaMalloc((void**)&mat, c->na*c->nacomp*sizeof(int));
	cudaMalloc((void**)&width, sizeof(int));
	cudaMalloc((void**)&result, c->na*c->na*sizeof(int));
	
	cudaMemcpy(result, res, c->na*c->na*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(mat, c->matrixrevisited, c->na*c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(width, &(c->nacomp), sizeof(int), cudaMemcpyHostToDevice);

	computeTwins<<<c->na,c->na>>>(mat, width, result);

	cudaMemcpy(res, result, c->na*c->na*sizeof(int), cudaMemcpyDeviceToHost);

	c->pointtorep=(int*)malloc(2*c->na*sizeof(int));
	c->pointtorepincomp=(int*)malloc(2*c->nacomp*sizeof(int));
	
	for (int i=0; i<2*c->na;i++){
		c->pointtorep[i]=-1;
	}

	for (int i=0;i<2*c->nacomp;i++){
		c->pointtorepincomp[i]=-1;
	}

	int* ptr;
	
	cudaMalloc((void**)&ptr, 2*c->na*sizeof(int));

	cudaMemcpy(ptr, c->pointtorep, 2*c->na*sizeof(int), cudaMemcpyHostToDevice);

	computeRepresentatives<<<1,c->na>>>(result, ptr);	
	
	cudaMemcpy(c->pointtorep, ptr, 2*c->na*sizeof(int), cudaMemcpyDeviceToHost);


	int* res2= (int*) malloc(c->nacomp*c->nacomp*sizeof(int));
	int* result2;

	int *reversedMatrix = (int*)malloc(c->na*c->nacomp*sizeof(int));
	for (int i =0; i<c->na; i++){
		for (int j=0; j<c->nacomp; j++){
			reversedMatrix[j*c->na+i]=c->matrixrevisited[i*c->nacomp+j];
		}
	}

	for (int i = 0; i< c-> nacomp*c->nacomp;i++)
		res2[i]=0;

	cudaMalloc((void**)&result2, c->nacomp*c->nacomp*sizeof(int));

	cudaMemcpy(result2, res2, c->nacomp*c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(mat, reversedMatrix, c->na*c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(width, &(c->na), sizeof(int), cudaMemcpyHostToDevice);

	computeTwins<<<c->nacomp,c->nacomp>>>(mat, width, result2);

	int* ptr2;

	cudaMalloc((void**)&ptr2, 2*c->nacomp*sizeof(int));

	cudaMemcpy(ptr2, c->pointtorepincomp, 2*c->nacomp*sizeof(int), cudaMemcpyHostToDevice);

	computeRepresentatives<<<1,c->nacomp>>>(result2, ptr2);

	cudaMemcpy(c->pointtorepincomp, ptr2, 2*c->nacomp*sizeof(int), cudaMemcpyDeviceToHost);
		
	for (int i = 0; i< 2*c->na; i++)
		c->pointtorep[i]=c->a[c->pointtorep[i]];
	for (int i = 0; i< 2*c->nacomp; i++)
		c->pointtorepincomp[i]=c->acomp[c->pointtorepincomp[i]];



	c->tc=(int*)malloc(c->na*sizeof(int));
	for (int i=0;i<c->na;i++)
		c->tc[i]=-1;
	c->complementtc=(int*)malloc(c->nacomp*sizeof(int));
	int cursor = 0;
	for (int i=0;i<c->na;i++){
		int here=0;
		for (int j=0;j<c->na;j++){
			if (c->pointtorep[2*i+1]==c->tc[j]){
				here=1;
				break;
			}
		}
		if (here==0){
			c->tc[cursor]=c->pointtorep[2*i+1];			
			cursor++;
		}
	}

	c->nrep=cursor;

	for (int i=0;i<c->nacomp;i++)
		c->complementtc[i]=-1;

	cursor=0;
	for (int i=0;i<c->nacomp;i++){
		int here=0;
		for (int j=0;j<c->nacomp;j++){
			if (c->pointtorepincomp[2*i+1]==c->complementtc[j]){
				here=1;
				break;
			}
		}
		if (here==0){
			c->complementtc[cursor]=c->pointtorepincomp[2*i+1];
			cursor++;
		}
	}
	c->nrepincomp=cursor;

	return EXIT_SUCCESS;

}



int secondpreprocess (cutdata* c, graph* g){


	c->lra = (int*) malloc (c->nrep*sizeof(int));
	c->lnra = (int*) malloc (c->nrepincomp*sizeof(int));
	c->lracard=1;
	c->lnracard=1;

	for (int i = 0; i<c->nrep; i++)
	c->lra[i]=0;
	for (int i = 0; i< c->nrepincomp; i++)
	c->lnra[i]=0;

	int *nextLevel;
	int *lastLevel=(int*)malloc(sizeof(int));
	lastLevel[0]=0;

	int sizeoflast=1;
	int sizeofnext=0;

	while (sizeoflast!=0){
		int* l = (int*)malloc(sizeoflast*c->nrep*c->nrep*sizeof(int));
		int* r = (int*)malloc(sizeoflast*c->nrep*c->nrepincomp*sizeof(int));

		for (int i = 0; i< sizeoflast; i ++){
			for (int j = 0; j<c->nrep; j++) {
				for (int k=0; k<c->nrep; k++){
					if (k==j)
						l[(i*c->nrep+j)*c->nrep+k]=1;
					else
						l[(i*c->nrep+j)*c->nrep+k]=c->lra[lastLevel[i]*c->nrep+k];
				}
			}
		}
		int* lc;
		int* rc;
		int* lrc;
		int* lrcard;
		int* mat;
		int* width;
		int* na;
		int* nacomp;
		int* a;
		int* acomp;
		int* rep;
		int* repc;


		cudaMalloc((void**)&lc, sizeoflast*c->nrep*c->nrep*sizeof(int));
		cudaMalloc((void**)&rc, sizeoflast*c->nrep*c->nrepincomp*sizeof(int));
		cudaMalloc((void**)&lrc, c->lnracard*c->nrepincomp*sizeof(int));
		cudaMalloc((void**)&mat, c->na * c-> nacomp * sizeof(int));
		cudaMalloc((void**)&lrcard, sizeof(int));
		cudaMalloc((void**)&width, sizeof(int));
		cudaMalloc((void**)&na, sizeof(int));
		cudaMalloc((void**)&nacomp, sizeof(int));
		cudaMalloc((void**)&a, c->na*sizeof(int));
		cudaMalloc((void**)&acomp, c->nacomp*sizeof(int));
		cudaMalloc((void**)&rep, c->nrep*sizeof(int));
		cudaMalloc((void**)&repc, c->nrepincomp*sizeof(int));

		cudaMemcpy(lc, l, sizeoflast*c->nrep*c->nrep*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(rc, r, sizeoflast*c->nrep*c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(lrc, c->lnra, c->lnracard*c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(lrcard, &(c->lnracard), sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(mat, c->matrixrevisited, c->na*c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(width, &(c->nrepincomp), sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(na, &(c->na), sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(nacomp, &(c->nacomp), sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(a, c->a, c->na*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(acomp, c->acomp, c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(rep, c->tc, c->nrep*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(repc, c->complementtc, c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);

		computeNeighboorhoods<<<sizeoflast,c->nrep>>>(lc, rc, lrc, lrcard, mat, width, na, nacomp, a, acomp, rep, repc);

		cudaMemcpy(l, lc, sizeoflast*c->nrep*c->nrep*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(r, rc, sizeoflast*c->nrep*c->nrepincomp*sizeof(int), cudaMemcpyDeviceToHost);
		int *ltemp = (int*)malloc((c->lracard+sizeoflast*c->nrep)*c->nrep*sizeof(int));
		int *lrtemp = (int*)malloc((c->lnracard+sizeoflast*c->nrep)*c->nrepincomp*sizeof(int));
		nextLevel = (int*)malloc(sizeoflast*c->nrep*c->nrep*sizeof(int));
		printf("a= ");
		for (int i=0; i<c->na; i++)
			printf("%d, ", c->a[i]);
		printf("\n");
		printf("l=\n");
		for (int i =0; i< sizeoflast* c->nrep; i++){
			for (int j = 0; j<c->nrep; j++){
				printf("%d, ", l[i*c->nrep+j]);
			}
			printf("\n");
		}
		printf("r=\n");
		for (int i =0; i<sizeoflast*c->nrep; i++){
			for (int j=0; j<c->nrepincomp; j++)
				printf("%d, ", r[i*c->nrepincomp+j]);
			printf("\n");
		}


		for (int i= 0; i< c->lracard*c->nrep;i++)
			ltemp[i]=c->lra[i];
		for (int i = 0; i<c->lnracard*c->nrepincomp;i++)
			lrtemp[i]=c->lnra[i];
		sizeofnext = 0;
		for (int i = 0; i < sizeoflast*c->nrep; i++){
			if (r[i*c->nrepincomp]!=-1){
				int alreadyin = 0;
				for (int j=0; j<c->lnracard+sizeofnext; j++){
					int id=1;
					for (int k= 0; k<c->nrepincomp; k++){
						if (r[i*c->nrepincomp+k]!=lrtemp[j*c->nrepincomp+k]){
							id=0;
							break;
						}
					}
					if (id==1){
						alreadyin=1;
						break;
					}
				}
				if (alreadyin ==0){
					for (int j =0; j<c->nrep; j++){
						nextLevel[sizeofnext*c->nrep+j]=l[i*c->nrep+j];
						ltemp[(c->lracard+sizeofnext)*c->nrep+j]=l[i*c->nrep+j];
					}
					for (int j =0; j<c->nrepincomp; j++){
						lrtemp[(c->lnracard+sizeofnext)*c->nrepincomp+j]=r[i*c->nrepincomp+j];
					}
					sizeofnext++;
				}
			}
		}
		c->lracard=c->lracard+sizeofnext;
		c->lnracard=c->lnracard+sizeofnext;
		c->lra=ltemp;
		c->lnra=lrtemp;



		lastLevel=nextLevel;
		sizeoflast=sizeofnext;
	}
	printf("comp\n");
	c->lracomp = (int*) malloc (c->nrepincomp*sizeof(int));
	c->lnracomp = (int*) malloc (c->nrep*sizeof(int));
	c->lracompcard=1;
	c->lnracompcard=1;

	for (int i = 0; i<c->nrepincomp; i++)
	c->lracomp[i]=0;
	for (int i = 0; i< c->nrep; i++)
	c->lnracomp[i]=0;

	nextLevel=NULL;
	lastLevel=NULL;
	lastLevel=(int*)malloc(sizeof(int));
	lastLevel[0]=0;

	sizeoflast=1;
	sizeofnext=0;

	while (sizeoflast!=0){
		int* l = (int*)malloc(sizeoflast*c->nrepincomp*c->nrepincomp*sizeof(int));
		int *r = (int*)malloc(sizeoflast*c->nrep*c->nrepincomp*sizeof(int));

		for (int i = 0; i< sizeoflast; i ++){
			for (int j = 0; j<c->nrepincomp; j++) {
				for (int k=0; k<c->nrepincomp; k++){
					if (k==j)
						l[(i*c->nrepincomp+j)*c->nrepincomp+k]=1;
					else
						l[(i*c->nrepincomp+j)*c->nrepincomp+k]=c->lracomp[lastLevel[i]*c->nrepincomp+k];
				}
			}
		}
		int* lc;
		int* rc;
		int* lrc;
		int* lrcard;
		int* mat;
		int* width;
		int* na;
		int* nacomp;
		int* a;
		int* acomp;
		int* rep;
		int* repc;
		int* revMatrix = (int*)malloc(c->na*c->nacomp*sizeof(int));

		for (int i = 0; i<c->na; i++){
			for (int j =0 ;j <c->nacomp; j++){
				revMatrix[j*c->na+i]=c->matrixrevisited[i*c->nacomp+j];
			}
		}


		cudaMalloc((void**)&lc, sizeoflast*c->nrepincomp*c->nrepincomp*sizeof(int));
		cudaMalloc((void**)&rc, sizeoflast*c->nrepincomp*c->nrep*sizeof(int));
		cudaMalloc((void**)&lrc, c->lnracompcard*c->nrep*sizeof(int));
		cudaMalloc((void**)&mat, c->nacomp * c-> na * sizeof(int));
		cudaMalloc((void**)&lrcard, sizeof(int));
		cudaMalloc((void**)&width, sizeof(int));
		cudaMalloc((void**)&na, sizeof(int));
		cudaMalloc((void**)&nacomp, sizeof(int));
		cudaMalloc((void**)&a, c->nacomp*sizeof(int));
		cudaMalloc((void**)&acomp, c->na*sizeof(int));
		cudaMalloc((void**)&rep, c->nrepincomp*sizeof(int));
		cudaMalloc((void**)&repc, c->nrep*sizeof(int));

		cudaMemcpy(lc, l, sizeoflast*c->nrepincomp*c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(rc, r, sizeoflast*c->nrep*c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(lrc, c->lnracomp, c->lnracompcard*c->nrep*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(lrcard, &(c->lnracompcard), sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(mat, revMatrix, c->na*c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(width, &(c->nrep), sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(na, &(c->nacomp), sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(nacomp, &(c->na), sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(a, c->acomp, c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(acomp, c->a, c->na*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(rep, c->complementtc, c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(repc, c->tc, c->nrep*sizeof(int), cudaMemcpyHostToDevice);

		computeNeighboorhoods<<<sizeoflast,c->nrepincomp>>>(lc, rc, lrc, lrcard, mat, width, na, nacomp, a, acomp, rep, repc);

		cudaMemcpy(l, lc, sizeoflast*c->nrepincomp*c->nrepincomp*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(r, rc, sizeoflast*c->nrep*c->nrepincomp*sizeof(int), cudaMemcpyDeviceToHost);
		int *ltemp = (int*)malloc((c->lracompcard+sizeoflast*c->nrepincomp)*c->nrepincomp*sizeof(int));
		int *lrtemp = (int*)malloc((c->lnracompcard+sizeoflast*c->nrepincomp)*c->nrep*sizeof(int));
		nextLevel = (int*)malloc(sizeoflast*c->nrepincomp*c->nrepincomp*sizeof(int));

		for (int i= 0; i< c->lracompcard*c->nrepincomp;i++)
			ltemp[i]=c->lracomp[i];
		for (int i = 0; i<c->lnracompcard*c->nrep;i++)
			lrtemp[i]=c->lnracomp[i];
		printf("a= ");
		for (int i =0; i<c->na; i++)
			printf("%d, ",c->a[i]);
		printf("\n");
		printf("l = \n");
		for (int i = 0; i<c->nrepincomp*sizeoflast; i++){
			for (int j = 0; j<c->nrepincomp; j++){
				printf ("%d, ", l[i*c->nrepincomp+j]);
			}
			printf("\n");
		}
		printf("r = \n");
		for (int i = 0; i<c->nrepincomp*sizeoflast; i++){
			for (int j = 0; j<c->nrep; j++)
				printf("%d, ", r[i*c->nrep+j]);
			printf("\n");
		}
		sizeofnext = 0;
		for (int i = 0; i < sizeoflast*c->nrepincomp; i++){
			if (r[i*c->nrep]!=-1){
				int alreadyin = 0;
				for (int j=0; j<c->lnracompcard+sizeofnext; j++){
					int id=1;
					for (int k= 0; k<c->nrep; k++){
						if (r[i*c->nrep+k]!=lrtemp[j*c->nrep+k]){
							id=0;
							break;
						}
					}
					if (id==1){
						alreadyin=1;
						break;
					}
				}
				if (alreadyin ==0){
					for (int j =0; j<c->nrepincomp; j++){
						nextLevel[sizeofnext*c->nrepincomp+j]=l[i*c->nrepincomp+j];
						ltemp[(c->lracompcard+sizeofnext)*c->nrepincomp+j]=l[i*c->nrepincomp+j];
					}
					for (int j =0; j<c->nrep; j++){
						lrtemp[(c->lnracompcard+sizeofnext)*c->nrep+j]=r[i*c->nrep+j];
					}
					sizeofnext++;
				}
			}
		}
		c->lracompcard=c->lracompcard+sizeofnext;
		c->lnracompcard=c->lnracompcard+sizeofnext;
		c->lracomp=ltemp;
		c->lnracomp=lrtemp;
		lastLevel=nextLevel;
		sizeoflast=sizeofnext;
		//printf("Yo\n");
	}



	return EXIT_SUCCESS;
}


int thirdpreprocess (cutdata* c, graph* g){

	c->m=(int*)malloc(c->lracard*c->nrep*sizeof(int));
	for (int i=0; i<c->lracard*c->nrep; i++)
		c->m[i]=-1;
	int* lr;
	int* ln;
	int* lrcard;
	int* tc;
	int* mg;
	int* mat;
	int* na;
	int* nacomp;
	int* a;
	int* acomp;
	int* rep;
	int* repc;
	int* width;
	int* r;

	cudaMalloc((void**)&lr, c->lracard*c->nrep*sizeof(int));
	cudaMalloc((void**)&ln, c->lnracard*c->nrepincomp*sizeof(int));
	cudaMalloc((void**)&lrcard, sizeof(int));
	cudaMalloc((void**)&tc, c->nrep*sizeof(int));
	cudaMalloc((void**)&mg, c->lracard*c->nrep*sizeof(int));
	cudaMalloc((void**)&mat, c->na*c->nacomp*sizeof(int));
	cudaMalloc((void**)&na, sizeof(int));
	cudaMalloc((void**)&nacomp, sizeof(int));
	cudaMalloc((void**)&a, c->na*sizeof(int));
	cudaMalloc((void**)&acomp, c->nacomp*sizeof(int));
	cudaMalloc((void**)&rep, c->nrep*sizeof(int));
	cudaMalloc((void**)&repc, c->nrepincomp*sizeof(int));
	cudaMalloc((void**)&width, sizeof(int));	
	cudaMalloc((void**)&r, c->nrepincomp*sizeof(int));

	cudaMemcpy(lr, c->lra, c->lracard*c->nrep*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ln, c->lnra, c->lnracard*c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(lrcard, &(c->lracard), sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(tc, c->tc, c->nrep*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(mg, c->m, c->lracard*c->nrep*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(mat, c->matrixrevisited, c->na*c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(na, &(c->na), sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(nacomp, &(c->nacomp), sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(a, c->a, c->na*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(acomp, c->acomp, c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(rep, c->tc, c->nrep*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(repc, c->complementtc, c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(width, &(c->nrepincomp), sizeof(int), cudaMemcpyHostToDevice);
	
	computeMatrix<<<c->lracard, c->nrep>>>(lr, ln, lrcard, tc, mg, mat, na, nacomp, a, acomp, rep, repc, width, r);

	cudaMemcpy(c->m, mg, c->lracard*c->nrep*sizeof(int), cudaMemcpyDeviceToHost);

	c->mcomp=(int*)malloc(c->lracompcard*c->nrepincomp*sizeof(int));
	for (int i=0; i<c->lracompcard*c->nrepincomp; i++)
		c->mcomp[i]=-1;

	cudaMalloc((void**)&lr, c->lracompcard*c->nrepincomp*sizeof(int));
	cudaMalloc((void**)&ln, c->lnracompcard*c->nrep*sizeof(int));
	cudaMalloc((void**)&lrcard, sizeof(int));
	cudaMalloc((void**)&tc, c->nrepincomp*sizeof(int));
	cudaMalloc((void**)&mg, c->lracompcard*c->nrepincomp*sizeof(int));
	cudaMalloc((void**)&mat, c->nacomp*c->na*sizeof(int));
	cudaMalloc((void**)&na, sizeof(int));
	cudaMalloc((void**)&nacomp, sizeof(int));
	cudaMalloc((void**)&a, c->nacomp*sizeof(int));
	cudaMalloc((void**)&acomp, c->na*sizeof(int));
	cudaMalloc((void**)&rep, c->nrepincomp*sizeof(int));
	cudaMalloc((void**)&repc, c->nrep*sizeof(int));
	cudaMalloc((void**)&width, sizeof(int));
	cudaMalloc((void**)&r, c->nrep*sizeof(int));

	int* reversedMatrix = (int*) malloc(c->na*c->nacomp*sizeof(int));
	for (int i = 0; i<c->na; i++){
		for (int j=0; j<c->nacomp; j++){
			reversedMatrix[j*c->na+i]=c->matrixrevisited[i*c->nacomp+j];
		}
	}
	

	cudaMemcpy(lr, c->lracomp, c->lracompcard*c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ln, c->lnracomp, c->lnracompcard*c->nrep*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(lrcard, &(c->lracompcard), sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(tc, c->complementtc, c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(mg, c->mcomp, c->lracompcard*c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(mat, reversedMatrix, c->na*c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(na, &(c->nacomp), sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(nacomp, &(c->na), sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(a, c->a, c->na*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(acomp, c->acomp, c->nacomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(rep, c->complementtc, c->nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(repc, c->tc, c->nrep*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(width, &(c->nrep), sizeof(int), cudaMemcpyHostToDevice);

	computeMatrix<<<c->lracompcard, c->nrepincomp>>>(lr, ln, lrcard, tc, mg, mat, na,nacomp, a, acomp, rep, repc, width, r);

	cudaMemcpy(c->mcomp, mg, c->lracompcard*c->nrepincomp*sizeof(int), cudaMemcpyDeviceToHost);

	return EXIT_SUCCESS;
}

void *threadAlgorithm ( void *arg){
	int i = *((int*)arg);
	int j = i;
	int treating=1;
	while (treating!=0){
		pthread_mutex_lock(&mut);
		if (tocompute==0){
			pthread_mutex_unlock(&mut);
			treating=0;
		}
		else {
			if (nodestocompute[i]->computed==0){
			//	printf("Thread %d computing %d\n", j, i);
				nodestocompute[i]->computed=2;
				pthread_mutex_unlock(&mut);
				stepalgorithm(nodestocompute[i],gwork);
			}
			else{
				pthread_mutex_unlock(&mut);
			//	printf("Thread %d: %d has already been computed\n", j, i);
			}
			i=(i+1)%nnodes;
			
			/*pthread_mutex_lock(&mut);
			printf("State of the art with %d left to compute\n", tocompute);
			for (int k = 0;k<nnodes; k++)
				printf("Node %d, computed = %d\n", k, nodestocompute[k]->computed);
			pthread_mutex_unlock(&mut);
			*/
		}

	}
	printf("Thread %d closing\n", j);
	pthread_exit(NULL);
}

int* toplevelalgorithm (dectree* t, graph* g, int n){

	if ((t->right==NULL)||(t->left==NULL)){
		return NULL;
	}

	nnodes=getnumberofnodes(t)-1;
	tocompute=0;
	nodestocompute=(dectree**)malloc(nnodes*sizeof(dectree*));
	if (t->right!=NULL)
		fillThevoid(t->right,g);
	else
		return toplevelalgorithm(t->left, g, n);

	if (t->left!=NULL)
		fillThevoid(t->left, g);
	
	pthread_t* threads = (pthread_t*)malloc(n*sizeof(pthread_t));

	for (int i=0; i<n; i++){
		if (pthread_create(&threads[i], NULL, threadAlgorithm, &i)){
			perror("pthread_create");
			return NULL;
		}
	}
	printf("Initiating threads\n");

	int finished=0;
	while (finished==0){
		pthread_mutex_lock(&mut);
		printf("----------------------------------Main Thread sees tocompute at %d\n",tocompute);
		if (tocompute<=0)
			finished=1;
		pthread_mutex_unlock(&mut);
	}
	printf("-----------------------------------------------Preparing to close\n");
	//for (int i=0; i<n; i++)
	//	pthread_join(threads[i],NULL);
	printf("--------------------------------------------------------closing\n");
	sleep(4);
	for (int i =0; i<nnodes; i++){
		printf("Data of node %d-----------------------\n",i);
		printf("a = ");
		for (int j=0; j<nodestocompute[i]->c.na; j++)
			printf("%d, ",nodestocompute[i]->c.a[j]);
		printf("\n");
		printf("acomp = ");
		for (int j=0; j<nodestocompute[i]->c.nacomp; j++)
			printf("%d, ", nodestocompute[i]->c.acomp[j]);
		printf("\n");
		printf("tc = ");
		for (int j=0; j<nodestocompute[i]->c.nrep; j++)
			printf("%d, ", nodestocompute[i]->c.tc[j]);
		printf("\n");
		printf("complementtc = ");
		for (int j=0; j<nodestocompute[i]->c.nrepincomp; j++)
			printf("%d, ", nodestocompute[i]->c.complementtc[j]);
		printf("\n");
		printf("pointtorep = \n");
		for (int j=0; j<nodestocompute[i]->c.na; j++)
			printf("%d->%d\n", nodestocompute[i]->c.pointtorep[2*j], nodestocompute[i]->c.pointtorep[2*j+1]);
		printf("pointtorepincomp = \n");
		for (int j=0; j<nodestocompute[i]->c.nacomp; j++)
			printf("%d->%d\n", nodestocompute[i]->c.pointtorepincomp[2*j], nodestocompute[i]->c.pointtorepincomp[2*j+1]);
		printf("matrixrevisited = \n");
		for (int j=0; j<nodestocompute[i]->c.na; j++){
			for (int k=0; k<nodestocompute[i]->c.nacomp; k++)
				printf("%d ", nodestocompute[i]->c.matrixrevisited[j*nodestocompute[i]->c.nacomp+k]);
			printf("\n");
		}
		printf("lra = \n");
		for (int j=0; j<nodestocompute[i]->c.lracard; j++){
			for (int k=0; k<nodestocompute[i]->c.nrep; k++)
				printf("%d, ", nodestocompute[i]->c.lra[j*nodestocompute[i]->c.nrep+k]);
			printf("\n");
		}
		printf("lnra = \n");
		for (int j=0; j<nodestocompute[i]->c.lnracard; j++){
			for (int k=0; k<nodestocompute[i]->c.nrepincomp; k++)
				printf("%d, ", nodestocompute[i]->c.lnra[j*nodestocompute[i]->c.nrepincomp+k]);
			printf("\n");
		}
		printf("lracomp = \n");
		for (int j=0; j<nodestocompute[i]->c.lracompcard; j++){
			for (int k=0; k<nodestocompute[i]->c.nrepincomp;k++)
				printf("%d, " , nodestocompute[i]->c.lracomp[j*nodestocompute[i]->c.nrepincomp+k]);
			printf("\n");
		}
		printf("lnracomp = \n");
		for (int j=0; j<nodestocompute[i]->c.lnracompcard; j++){
			for (int k=0; k<nodestocompute[i]->c.nrep; k++)
				printf("%d, ", nodestocompute[i]->c.lnracomp[j*nodestocompute[i]->c.nrep+k]);
			printf("\n");
		}
		printf("m =\n");
		for (int j=0; j<nodestocompute[i]->c.lracard; j++){
			for (int k=0; k<nodestocompute[i]->c.nrep; k++)
				printf("%d ", nodestocompute[i]->c.m[j*nodestocompute[i]->c.nrep+k]);
			printf("\n");
		}
		printf("mcomp =\n");
		for (int j=0; j<nodestocompute[i]->c.lracompcard; j++){
			for (int k=0; k<nodestocompute[i]->c.nrepincomp; k++)
				printf("%d ", nodestocompute[i]->c.mcomp[j*nodestocompute[i]->c.nrepincomp+k]);
			printf("\n");
		}
	}

	cutdata *c = (cutdata*)malloc(2*sizeof(cutdata));
	cutdata c1 = t->left->c;
	cutdata c2= t->right->c;
	int size=-1;	
	int amax=-1;
	int acmax=-1;
	int bmax=-1;
	int bcmax=-1;

	for (int i=0;i<c1.lracard;i++){
		for (int j=0;j<c2.lracard;j++){
			pointset p;
			p.size=0;
			p.members=(int*)malloc(c1.nrep*sizeof(int));
			for (int k=0; k<c1.nrep; k++){
				if (c1.lra[i*c1.nrep+k]==1){
					p.size++;
					p.members[p.size-1]=c1.tc[k];
				}
			} 
			pointset q;
			q.size=0;
			q.members=(int*)malloc(c2.nrep*sizeof(int));
			for (int k=0; k<c2.nrep; k++){
				if (c2.lra[i*c2.nrep+k]==1){
					q.size++;
					q.members[q.size-1]=c2.tc[k];
				}
			}
			int indexac=0;
			int indexbc=0;
			for (int k = 0; k<c2.nrep; k++){
				if (c2.lra[c2.nrep*j+k]==1){
					int rep = -1;
					for (int l=0; l< c1.nacomp; l++){
						if (c1.pointtorepincomp[2*l]==c2.tc[k]){
							rep=c1.pointtorepincomp[2*l+1];
							break;
						}
					}
					for (int l = 0; l< c1.nrepincomp; l++){
						if (c1.complementtc[l]==rep){
							rep=l;
							break;
						}
					}
				indexac = c1.mcomp[indexac*c1.nrepincomp+k];
				}
			}
			for (int k = 0; k<c1.nrep; k++){
				if (c1.lra[c1.nrep*j+k]==1){
					int rep = -1;
					for (int l=0; l< c2.nacomp; l++){
						if (c2.pointtorepincomp[2*l]==c1.tc[k]){
							rep=c2.pointtorepincomp[2*l+1];
							break;
						}
					}
					for (int l = 0; l< c2.nrepincomp; l++){
						if (c2.complementtc[l]==rep){
							rep=l;
							break;
						}
					}
				indexbc = c2.mcomp[indexbc*c2.nrepincomp+k];
				}
			}

			if ((size==-1)&&(c1.tab[i*c1.lracompcard+indexac]!=-1)&&(c2.tab[j*c2.lracompcard+indexbc]!=-1)){
				size=c1.tab[i*c1.lracompcard+indexac]+c2.tab[j*c2.lracompcard+indexbc];
				amax=i;
				acmax=indexac;
				bmax=j;
				bcmax=indexbc;
			}
			else {
				if ((c1.tab[i*c1.lracompcard+indexac]!=-1)&&(c2.tab[j*c2.lracompcard+indexbc]!=-1)){
					if (c1.tab[i*c1.lracompcard+indexac]+c2.tab[j*c2.lracompcard+indexbc]<size){
						size=c1.tab[i*c1.lracompcard+indexac]+c2.tab[j*c2.lracompcard+indexbc];
						amax=i;
						acmax=indexac;
						bmax=j;
						bcmax=indexbc;
					}
				}
			}
	
		}
		
	}
	printf("bmax= %d, bcmax= %d, amax= %d, acmax= %d\n", bmax, bcmax, amax, acmax);
	printf("c1.tab=\n");
	for (int i= 0; i<c1.lracard; i++){
		for (int j=0; j<c1.lracompcard; j++){
			printf("%d ", c1.tab[i*c1.lracompcard+j]);
		}
		printf("\n");
	}
	printf("c2.tab=\n");
	for (int i= 0; i<c2.lracard; i++){
		for (int j=0; j<c2.lracompcard; j++){
			printf("%d ", c2.tab[i*c2.lracompcard+j]);
		}
		printf("\n");
	}
	int * sol = (int*)malloc((c2.tab[bmax*c2.lracompcard+bcmax]+c1.tab[amax*c1.lracompcard+acmax])*sizeof(int));
	int* left = (int*)malloc(c1.tab[amax*c1.lracompcard+acmax]*sizeof(int));
	int* right= (int*)malloc(c2.tab[bmax*c2.lracompcard+bcmax]*sizeof(int));
	left = computeDS (t->left, c1.tab[amax*c1.lracompcard+acmax], amax, acmax);
	right= computeDS (t->right, c2.tab[bmax*c2.lracompcard+bcmax], bmax, bcmax);

	for (int i = 0; i< c1.tab[amax*c1.lracompcard+acmax]; i++)
		sol[i]=left[i];
	
	for (int i = 0; i< c2.tab[bmax*c2.lracompcard+bcmax]; i++)
		sol[i+c1.tab[amax*c1.lracompcard+acmax]]=right[i];
	

	for (int i = 0; i<c2.tab[bmax*c2.lracompcard+bcmax]+c1.tab[amax*c1.lracompcard+acmax];i++)
		printf("%d, ", sol[i]);
	
	printf("\n");
	return sol;
}


int stepalgorithm (dectree* t, graph* g){
	
	if ((t->right==NULL)||(t->left==NULL)){
		t->c = cutThatTree (g, t);
		firstpreprocess (g,&(t->c));
		secondpreprocess (&(t->c), g);
		thirdpreprocess (&(t->c), g);

		t->c.tab = (int*)malloc(t->c.lracard*t->c.lracompcard*sizeof(int));
		t->c.tab[0]=-1;
		t->c.tab[1]=0;
		t->c.tab[2]=1;
		t->c.tab[3]=1;
		t->computed=1;
		pthread_mutex_lock(&mut);
		tocompute--;
		pthread_mutex_unlock(&mut);
	}

	else {
		if ((t->right->computed==1)&&(t->left->computed==1)){
			t->c = cutThatTree (g, t);
			firstpreprocess (g,&(t->c));
			secondpreprocess (&(t->c), g);
			thirdpreprocess (&(t->c), g);
			t->c.tab = (int*)malloc(t->c.lracard*t->c.lracompcard*sizeof(int));
			t->c.box = (int*)malloc(t->c.lracard*t->c.lracompcard*6*sizeof(int));
			for (int i = 0; i<t->c.lracard*t->c.lracompcard;i++){
				t->c.tab[i]=-1;
				t->c.box[6*i]=-1;
				t->c.box[6*i+1]=-1;
				t->c.box[6*i+2]=-1;
				t->c.box[6*i+3]=-1;
				t->c.box[6*i+4]=-1;
				t->c.box[6*i+5]=-1;
			}

			int *tmptab = (int*)malloc(5*t->right->c.lracard*t->left->c.lracard*t->c.lracompcard*sizeof(int));
			for (int i = 0; i< 5*t->right->c.lracard*t->left->c.lracard*t->c.lracompcard; i++)
				tmptab[i]=-1;

			int* tabg;
			int* lra;
			int* lrb;
			int* lrw;
			int* lracard;
			int* lrbcard;
			int* lrwcard;
			int* lnracard;
			int* lnrbcard;
			int* lnrwcard;
			int* mw;
			int* macomp;
			int* mbcomp;
			int* nrepa;
			int* nrepb;
			int* nrepw;
			int* repacomp;
			int* repbcomp;
			int* repw;
			int* taba;
			int* tabb;
			int* ptrac;
			int* ptrbc;
			int* ptrw;
			int* nacomp;
			int* nbcomp;
			int* nw;

			cudaMalloc((void**)&tabg, 5*t->right->c.lracard*t->left->c.lracard*t->c.lracompcard*sizeof(int)); 
			cudaMalloc((void**)&lra, t->left->c.lracard*t->left->c.nrep*sizeof(int)); 
			cudaMalloc((void**)&lrb, t->right->c.lracard*t->right->c.nrep*sizeof(int));
			cudaMalloc((void**)&lrw, t->c.lracompcard*t->c.nrepincomp*sizeof(int));
			cudaMalloc((void**)&lracard, sizeof(int));
			cudaMalloc((void**)&lrbcard, sizeof(int));
			cudaMalloc((void**)&lrwcard, sizeof(int));
			cudaMalloc((void**)&lnracard, sizeof(int));
			cudaMalloc((void**)&lnrbcard, sizeof(int));
			cudaMalloc((void**)&lnrwcard, sizeof(int));
			cudaMalloc((void**)&mw, t->c.lracard*t->c.nrep*sizeof(int));
			cudaMalloc((void**)&macomp, t->left->c.lracompcard*t->left->c.nrepincomp*sizeof(int));
			cudaMalloc((void**)&mbcomp, t->right->c.lracompcard*t->right->c.nrepincomp*sizeof(int));
			cudaMalloc((void**)&nrepa, sizeof(int));	
			cudaMalloc((void**)&nrepb, sizeof(int));
			cudaMalloc((void**)&nrepw, sizeof(int));
			cudaMalloc((void**)&repacomp, t->left->c.nrepincomp*sizeof(int));
			cudaMalloc((void**)&repbcomp, t->right->c.nrepincomp*sizeof(int));
			cudaMalloc((void**)&repw, t->c.nrep*sizeof(int));
			cudaMalloc((void**)&taba, t->left->c.lracard*t->left->c.lracompcard*sizeof(int));
			cudaMalloc((void**)&tabb, t->right->c.lracard*t->right->c.lracompcard*sizeof(int));
			cudaMalloc((void**)&ptrac, 2*t->left->c.nacomp*sizeof(int));
			cudaMalloc((void**)&ptrbc, 2*t->right->c.nacomp*sizeof(int));
			cudaMalloc((void**)&ptrw, 2*t->c.na*sizeof(int));
			cudaMalloc((void**)&nacomp, sizeof(int));
			cudaMalloc((void**)&nbcomp, sizeof(int));
			cudaMalloc((void**)&nw, sizeof(int));

			cudaMemcpy(tabg, tmptab, 5*t->right->c.lracard*t->left->c.lracard*t->c.lracompcard*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(lra, t->left->c.lra, t->left->c.lracard*t->left->c.nrep*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(lrb, t->right->c.lra, t->right->c.lracard*t->right->c.nrep*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(lrw, t->c.lracomp, t->c.lracompcard*t->c.nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(lracard, &(t->left->c.lracard), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(lrbcard, &(t->right->c.lracard), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(lrwcard, &(t->c.lracompcard), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(lnracard, &(t->left->c.lracompcard), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(lnrbcard, &(t->right->c.lracompcard), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(lnrwcard, &(t->c.lnracompcard), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(mw, t->c.m, t->c.lracard*t->c.nrep*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(macomp, t->left->c.mcomp, t->left->c.lracompcard*t->left->c.nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(mbcomp, t->right->c.mcomp, t->right->c.lracompcard*t->right->c.nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(nrepa, &(t->left->c.nrepincomp), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(nrepb, &(t->right->c.nrepincomp), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(nrepw, &(t->c.nrep), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(repacomp, t->left->c.complementtc, t->left->c.nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(repbcomp, t->right->c.complementtc, t->right->c.nrepincomp*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(repw, t->c.tc, t->c.nrep*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(taba, t->left->c.tab, t->left->c.lracard*t->left->c.lracompcard*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(tabb, t->right->c.tab, t->right->c.lracard*t->right->c.lracompcard*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(ptrac, t->left->c.pointtorepincomp, 2*t->left->c.nacomp*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(ptrbc, t->right->c.pointtorepincomp, 2*t->right->c.nacomp*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(ptrw, t->c.pointtorep, 2*t->right->c.na*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(nacomp, &(t->left->c.nacomp), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(nbcomp, &(t->right->c.nacomp), sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(nw, &(t->c.na), sizeof(int), cudaMemcpyHostToDevice);

			computeAlgorithm <<<t->c.lracompcard, t->left->c.lracard*t->right->c.lracard>>> (tabg, lra, lrb, lrw, lracard, lrbcard, lrwcard, lnracard, lnrbcard, lnrwcard, mw, macomp, mbcomp, nrepa, nrepb, nrepw, repacomp, repbcomp, nrepw, taba, tabb, ptrac, ptrbc, ptrw, nacomp, nbcomp, nw);		

			cudaMemcpy(tmptab, tabg, 5*t->right->c.lracard*t->left->c.lracard*t->c.lracompcard*sizeof(int), cudaMemcpyDeviceToHost);			

			for (int i =0; i<t->c.lracompcard; i++){
				for (int j = 0; j<t->left->c.lracard; j++){
					for (int k = 0; k<t->right->c.lracard; k++){
					//	printf(" SALUT %d/%d, %d/%d, %d, %d\n", i, t->c.lracompcard, j, t->left->c.lracard, k, t->right->c.lracard);
						if (t->c.tab[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*t->c.lracompcard+i]==-1){
							t->c.tab[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*t->c.lracompcard+i]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+3]+tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+4];
							t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*6*t->c.lracompcard+6*i]=j;
							t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*6*t->c.lracompcard+6*i+1]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+3]+tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+1];
							t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*2*t->c.lracompcard+6*i+2]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+3];
							t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*2*t->c.lracompcard+2*i+3]=k;
							t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*2*t->c.lracompcard+2*i+4]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+2];
							t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*2*t->c.lracompcard+2*i+5]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+4];
						}
						else {
							if ((tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+3]!=-1)&&(tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+4]!=-1)&&(tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+3]+tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+4]<t->c.tab[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*t->c.lracompcard+i])){
								t->c.tab[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*t->c.lracompcard+i]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+3]+tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+4];
								t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*6*t->c.lracompcard+6*i]=j;
								t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*6*t->c.lracompcard+6*i+1]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+3]+tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+1];
								t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*2*t->c.lracompcard+6*i+2]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+3];
								t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*2*t->c.lracompcard+2*i+3]=k;
								t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*2*t->c.lracompcard+2*i+4]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+2];
								t->c.box[tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5]*2*t->c.lracompcard+2*i+5]=tmptab[5*t->left->c.lracard*t->right->c.lracard*i+j*5*t->right->c.lracard+k*5+4];
							}
						}
					}
				}
				t->computed=1;
				pthread_mutex_lock(&mut);
				tocompute--;
				pthread_mutex_unlock(&mut);
			}
		}	
		else {
			t->computed=0;
		//	printf("Can't compute for the moment\n");
		}	
	}

	return EXIT_SUCCESS;
}

int* computeDS (dectree* t, int much, int aleft, int acleft){
	int* sol = (int*)malloc(much*sizeof(int));

	if ((t->left==NULL)&&(t->right==NULL)&&(much==1)){
		sol[0]=t->c.tc[0];
	}

	else if (much!=0){
		int *left=(int*)malloc(t->c.box[6*t->c.lracompcard*aleft+6*acleft+2]);
		left = computeDS (t->left, t->c.box[6*t->c.lracompcard*aleft+6*acleft+2], t->c.box[6*t->c.lracompcard*aleft+6*acleft+1], t->c.box[6*t->c.lracompcard*aleft+6*acleft+1]);
	

		int *right=(int*)malloc(t->c.box[6*t->c.lracompcard*aleft+6*acleft+5]);
		right = computeDS (t->right, t->c.box[6*t->c.lracompcard*aleft+6*acleft+5], t->c.box[6*t->c.lracompcard*aleft+6*acleft+3], t->c.box[6*t->c.lracompcard*aleft+6*acleft+4]);
		for (int i = 0; i<t->c.box[6*t->c.lracompcard*aleft+6*acleft+2]; i++){
			sol[i]=left[i];
		}
		for (int i = 0; i<t->c.box[6*t->c.lracompcard*aleft+6*acleft+5]; i++){
			sol[t->c.box[6*t->c.lracompcard*aleft+6*acleft+2]+i]=right[i];
		}
	
	}

	
	return sol;
}

int getBW (dectree t, graph g){
	int bwmax=-1;
	if ((t.right==NULL)||(t.left==NULL))
		bwmax=-1;
	else {
		cutdata c1 = cutThatTree (&g, &t);

		firstpreprocess (&g,&c1);
		secondpreprocess (&c1, &g);

		cutdata c2 = cutThatTree (&g, &t);

		firstpreprocess (&g,&c2);
		secondpreprocess (&c2, &g);
		
		int p  = getBW(*(t.left), g);
		int q  = getBW(*(t.right), g);
		

		if (p>bwmax)
			bwmax=p;
		if (q>bwmax)
			bwmax=q;
		if (c1.lracard>bwmax)
			bwmax=c1.lracard;
		if (c2.lracard>bwmax)
			bwmax=c2.lracard;
	}
	return bwmax;
}
