#include "../include/treeprimitives.h"
#include "../include/algorithms.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <errno.h>
#include <sys/stat.h>


int getallleaves(dectree t, int *list){
	int n=0;
	int *lleft;
	int *lright;
	int nleft=0;
	int nright=0;
	if (t.left!=NULL){
		nleft= getnumberofleaves (*(t.left));
		lleft=(int*)malloc(nleft*sizeof(int));
		getallleaves(*(t.left), lleft);
	}

	if (t.right!=NULL){
		nright=getnumberofleaves(*(t.right));
		lright=(int*)malloc(nright*sizeof(int));
		getallleaves(*(t.right),lright);
	}

	if ((t.left==NULL)&&(t.right==NULL)){
		list[0]=t.label;
		n+=1;
	}
	else {
		n=nleft+nright;
		for (int i=0; i<nleft;i++)
			list[i]=lleft[i];

		for (int i=0;i<nright;i++)
			list[nleft+i]=lright[i];
	}
	return n;
}


int getnumberofleaves (dectree t){
	int n=0;
	int nleft=0;
	int nright=0;

	if (t.left!=NULL)
		nleft=getnumberofleaves(*(t.left));
	if (t.right!=NULL)
		nright=getnumberofleaves(*(t.right));
	if ((t.left==NULL)&&(t.right==NULL))
		n+=1;
	else 
		n=nleft+nright;
	return n;
}

dectree *generateTree (pointset p, graph g, int verthor){
	dectree *t;
	t = (dectree*)malloc(sizeof(dectree));
	pointset p1;
	pointset p2;
	p1.size=0;
	p2.size=0;
	p1.members=(int*)malloc(p.size*sizeof(int));
	p2.members=(int*)malloc(p.size*sizeof(int));
	int vnext=0;
	if (p.size==0)
		return NULL;
	else if (p.size==1){
		(*t).label=p.members[0];
		(*t).left=NULL;
		(*t).right=NULL;
	}
	else if (p.size==2){
		(*t).label=-1;
		(*t).left=(dectree*)malloc(sizeof(dectree));
		(*t).right=(dectree*)malloc(sizeof(dectree));
		dectree *t1;
		dectree *t2;
		t1=(dectree*)malloc(sizeof(dectree));
		t2=(dectree*)malloc(sizeof(dectree));
		(*t).left=t1;
		(*t).right=t2;
		(*t1).label=p.members[0];
		(*t1).left=NULL;
		(*t1).right=NULL;
		(*t2).label=p.members[1];
		(*t2).left=NULL;
		(*t2).right=NULL;
	}
	else {
		if (verthor==0){
			int xmoy=0;
			for (int i=0;i<p.size;i++)
				xmoy+=g.pos[2*p.members[i]];
			xmoy=xmoy/p.size;
			
			for (int i=0;i<p.size;i++){
				if (g.pos[2*p.members[i]]<xmoy){
					p1.size++;
					p1.members[p1.size-1]=p.members[i];
				}
				else {
					p2.size++;
					p2.members[p2.size-1]=p.members[i];
				}
			}
			if (p1.size==0){
				p1.size++;
				p1.members[0]=p2.members[p2.size-1];
				p2.size--;
			}
			if (p2.size==0){
				p2.size++;
				p2.members[0]=p1.members[p1.size-1];
				p1.size--;
			}
			vnext=1;
		}

		if (verthor==1){
			int ymoy=0;
			for (int i=0;i<p.size;i++)
				ymoy+=g.pos[2*p.members[i]+1];
			ymoy=ymoy/p.size;
			
			for (int i=0;i<p.size;i++){
				if (g.pos[2*p.members[i]+1]<ymoy){
					p1.size++;
					p1.members[p1.size-1]=p.members[i];
				}
				else {
					p2.size++;
					p2.members[p2.size-1]=p.members[i];
				}
			}
			if (p1.size==0){
				p1.size++;
				p1.members[0]=p2.members[p2.size-1];
				p2.size--;
			}
			if (p2.size==0){
				p2.size++;
				p2.members[0]=p1.members[p1.size-1];
				p1.size--;
			}
			vnext=0;
		}

		(*t).label=-1;
		(*t).left=(dectree*)malloc(sizeof(dectree));
		(*t).left=generateTree(p1,g,vnext);
		(*t).right=(dectree*)malloc(sizeof(dectree));
		(*t).right=generateTree(p2,g,vnext);
	}
	return t;
}


pointset getCandidates (pointset left, pointset right, graph g){
	pointset cand;
	cand.size=0;
	cand.members=(int*)malloc(right.size*sizeof(int));
	pointset nleft;
	nleft.size=0;
	nleft.members=(int*)malloc((right.size+left.size)*sizeof(int));
	for (int j=0;j<right.size;j++){
		for (int i=0;i<left.size;i++){	
			if (g.matrix[left.members[i]*g.size+right.members[j]]==1){
				nleft.size++;
				nleft.members[nleft.size-1]=right.members[j];
				break;
			}
		}
	}

	for (int i=0;i<left.size;i++){
		nleft.size++;
		nleft.members[nleft.size-1]=left.members[i];
	}

	pointset nnleft;
	nnleft.size=0;
	nnleft.members=(int*)malloc((right.size+left.size)*sizeof(int));

	for (int i=0;i<right.size;i++){
		for (int j=0;j<nleft.size;j++){
			if (g.matrix[right.members[i]*g.size+nleft.members[j]]==1){
				nnleft.size++;
				nnleft.members[nnleft.size-1]=right.members[i];
				break;
			}
		}
	}

	return nnleft;
		
}


setwithinsets incrementun(graph g, pointset x, setwithinsets unx, int v){
	setwithinsets unv;
	unv.size=0;
	unv.set=(pointset*)malloc((g.size*g.size*g.size*g.size+g.size*g.size*g.size)*sizeof(pointset));

	for (int i=0; i<unx.size; i++){
		pointset s;
		s.size=0;
		s.members=(int*)malloc(g.size*sizeof(int));
	
		for (int j=0; j<unx.set[i].size; j++){
			if (unx.set[i].members[j]!=v){
				s.size++;
				s.members[s.size-1]=unx.set[i].members[j];
			}			
		}

		int alreadyin = 0;

		for (int j=0; j<unv.size; j++){
			if (s.size==unv.set[j].size){
				int common=0;
				for (int k=0;k<s.size;k++){
					for (int l=0;l<unv.set[j].size;l++){
						if (s.members[k]==unv.set[j].members[l]){
							common++;
							break;
						}
					}
				}
				if (common==s.size){
					alreadyin = 1;
					break;
				}
			}
		}
		if (alreadyin==0){
			unv.size++;
			unv.set[unv.size-1]=s;
		}

		pointset t;
		t.size=0;
		t.members=(int*)malloc(g.size*sizeof(int));

		for (int j=0; j<g.size; j++){
			if (j!=v){
				int inX=0;
				for (int k=0; k<x.size;k++){
					if(x.members[k]==j){
						inX=1;
						break;
					}
				}
				if (inX==0){
					if (g.matrix[j*g.size+v]==1){
						t.size++;
						t.members[t.size-1]=j;
					}
				}
			}
		}

		alreadyin = 0;

		for (int j=0; j<unv.size; j++){
			if (t.size==unv.set[j].size){
				int common=0;
				for (int k=0;k<t.size;k++){
					for (int l=0;l<unv.set[j].size;l++){
						if (t.members[k]==unv.set[j].members[l]){
							common++;
							break;
						}
					}
				}
				if (common==s.size){
					alreadyin = 1;
					break;
				}
			}
		}
		if (alreadyin==0){
			unv.size++;
			unv.set[unv.size-1]=t;
		}
		
	}

	return unv;
}


pointset incrementalUNheuristic (graph g, int init){

	pointset dec;
	dec.size=1;
	dec.members=(int*)malloc(g.size*sizeof(int));
	dec.members[0]=init;
	pointset left;
	left.size=1;
	left.members=(int*)malloc(g.size*sizeof(int));
	pointset right;
	right.size=0;
	right.members=(int*)malloc(g.size*sizeof(int));
	left.members[0]=init;

	for (int i=0;i<g.size;i++){
		if (i!=init){
			right.size++;
			right.members[right.size-1]=i;
		}
	}

	setwithinsets unleft;
	unleft.size=2;
	unleft.set=(pointset*)malloc((g.size*g.size*g.size*g.size+g.size*g.size*g.size)*sizeof(pointset));
	unleft.set[0].size=0;
	unleft.set[1].size=0;
	unleft.set[1].members=(int*)malloc(g.size*sizeof(int));
	
	for (int i=0; i<g.size;i++){
		if ((g.matrix[init*g.size+i]==1)||(i!=init)){
			unleft.set[1].size++;
			unleft.set[1].members[unleft.set[1].size-1]=i;
		}
	}

	while (right.size!=0){
		if (right.size==1){
			dec.size++;
			dec.members[dec.size-1]=right.members[0];
			right.size--;
		}
		else {
			pointset candidates = getCandidates (left, right, g);
			int chosen=-1;
			setwithinsets unchosen;
			unchosen.size=0;
			unchosen.set=(pointset*)malloc(200*(g.size*g.size*g.size*g.size+g.size*g.size*g.size)*sizeof(pointset));

			if (candidates.size==0){
				candidates=right;
			}
			for (int i=0; i<candidates.size; i++){
				setwithinsets unv;
				unv.size=0;
				unv.set=(pointset*)malloc((g.size*g.size*g.size*g.size+g.size*g.size*g.size)*sizeof(pointset));
				unv = incrementun(g, left, unleft, candidates.members[i]);
				if ((chosen==-1)||(unv.size<unchosen.size)){
					chosen=candidates.members[i];
					unchosen=unv;
				}
			}

			dec.size++;
			dec.members[dec.size-1]=chosen;
			left.size++;
			left.members[left.size-1]=chosen;
			pointset newright;
			newright.size=0;
			newright.members=(int*)malloc(g.size*sizeof(int));
			for (int i=0;i<right.size;i++){
				if (right.members[i]!=chosen){
					newright.size++;
					newright.members[newright.size-1]=right.members[i];
				}
			}
			right = newright;
			unleft=unchosen;
		}
	}

	return dec;
}


dectree *generateTreeBWstep (graph g, pointset dec, int i){
	dectree *t;
	t=(dectree*)malloc(sizeof(dectree));

	if (dec.size-i==0)
		return NULL;

	if (dec.size-1==i){
		(*t).label=dec.members[dec.size-i-1];
		(*t).right=NULL;
		(*t).left=NULL;
	}
	else {
		dectree *tleft;
		tleft=(dectree*)malloc(sizeof(dectree));
		dectree *tright;
		tright=(dectree*)malloc(sizeof(dectree));
		(*tright).label=dec.members[dec.size-i-1];
		(*tright).left=NULL;
		(*tright).right=NULL;
		tleft=generateTreeBWstep (g, dec, i+1);
		(*t).left=tleft;
		(*t).right=tright;
		(*t).label=-1; 
	}
	return t;
}


dectree *generateTreeBW (graph g){
	dectree **t = (dectree**)malloc(g.size*sizeof(dectree*));
	int *bw = (int*)malloc(g.size*sizeof(int));
	for (int i=0; i<g.size; i++){
		pointset dec=incrementalUNheuristic (g, i);
		t[i]=(dectree*)malloc(sizeof(dectree));
		dectree *tleft;
		tleft=(dectree*)malloc(sizeof(dectree));
		dectree *tright;
		tright=(dectree*)malloc(sizeof(dectree));
		(*tright).label=dec.members[dec.size-1];
		(*tright).left=NULL;
		(*tright).right=NULL;
		tleft=generateTreeBWstep (g, dec, 1);
		t[i]->left=tleft;
		t[i]->right=tright;
		t[i]->label=-1;
		bw[i]=getBW(t[i]->left, &g);
		bw[i]+=getBW(t[i]->right, &g);
	} 

	int size=bw[0];
	int min=0;

	for (int i=1; i<g.size; i++){
		if (bw[i]<size){
			size=bw[i];
			min=i;
		}
	}
	return t[min];
}


int printTree (dectree t){
	printf("Un arbre");
	if (t.label!=-1)
		printf(" de label %d.", t.label);
	else {
		if (t.left!=NULL){
			printf(" dont le fils gauche est {");
			printTree(*(t.left));
			printf("}");
		}
		if (t.right != NULL){
			printf(" et dont le fils droit est {");
			printTree(*(t.right));
			printf("}");
		}
	}
	printf("\n");
	return EXIT_SUCCESS;
}
