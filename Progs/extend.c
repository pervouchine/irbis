//	Copyright 2011,2012 Dmitri Pervouchine (dp@crg.eu)
//	This file is a part of the IRBIS package.
//	IRBIS package is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//	
//	IRBIS package is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with IRBIS package.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#define MAX2(A,B) ((A<B) ? (B) : (A))
#define MIN2(A,B) ((A<B) ? (A) : (B))

#define MAX3(A,B,C) MAX2((A),MAX2(B,C))
#define MAXBUFFLENGTH 1000000
#define GAP -10
#define INFTY 65535

#define imatrix(A,X,Y) (((X)>=0 && (Y)>=0 && (Y)<=(X)) ? A[X][Y] : -INFTY)
#define zmatrix(A,X,Y) (((X)>=0 && (Y)>=0) ? A[X][Y] : 0)

char buffr[100];

char* sc(int i) {
    if(i<0) return((char*)"*");
    sprintf(buffr,"%i",i);
    return(buffr);
} 

int bp(char, char);
int gap(int);
int debug = 0;
int min_helix = 3;
int max_loop = 4;

void extend(char* left, int n, char *right, int m) {
    int i,j,k,d;
    int l,g;
    int **V;
    int **W;
    int *M;
    int s;
    char c;

    int maxv, maxi, maxd;

    l = 2*MIN2(n,m);

    V = (int**)malloc(sizeof(int*)*(l+10));
    W = (int**)malloc(sizeof(int*)*(l+10));
    M =  (int*)malloc(sizeof(int)*(l+10));

    maxv = -INFTY;
    maxi = maxd = -1;

    if(debug) {
	c = left[n];left[n]=0;
	fprintf(stderr,"Seq1:[%s]\n",left);
	left[n]=c;
	c = right[m];right[m]=0;
	fprintf(stderr,"Seq2:[%s]\n",right);
	right[m]=c;
    }

    for(d=0;d<l;d++) {
	V[d] = (int*)malloc(sizeof(int)*(d+10));
	W[d] = (int*)malloc(sizeof(int)*(d+10));
	M[d] = -INFTY;
	if(debug) fprintf(stderr, "%i\t",d);
	for(i=0;i<=d; i++) {
	    if(i<n && (j=d-i)<m) {
		s=(i==min_helix-1 && j==min_helix-1) ? 0 : imatrix(W,d-2*min_helix,i-min_helix);
	    	for(k=0;k<min_helix;k++) {
		    s+=(i-k>=0 && j-k>=0) ? bp(left[i-k],right[m-1-(j-k)]) : -INFTY;
	    	}
		if(debug) fprintf(stderr, "[%i]",s);
	    	V[d][i] = MAX2(s, imatrix(V,d-2,i-1)) + bp(left[i],right[m-1-j]);
	    	s = -INFTY;
	    	for(g=1;g<=max_loop;g++) {
		    for(k=0;k<=g;k++) {
	   	    	if(2*k<=max_loop && 2*(g-k)<=max_loop) s = MAX2(imatrix(V,d-g,i-k) + gap(g),s);
		    }
	    	}
	    	W[d][i]=s; 
		if(V[d][i]>maxv) {
		    maxv = V[d][i];
		    maxd = d;
		    maxi = i;
		}
		M[d] = MAX3(M[d], V[d][i], W[d][i]);
		if(debug) fprintf(stderr, "(%i %i %c %c %s %s)",i,j,left[i],right[m-1-j],sc(V[d][i]),sc(W[d][i]));
	    }
	}
	s = (d-2*min_helix>=0) ? M[d-2*min_helix] : INFTY;
	for(g=0;g<=max_loop;g++) {
	    s = MAX2(((d-g>=0) ? M[d-g] : INFTY), s);
	}
	if(debug) fprintf(stderr, "\n");
	if(s<0) break;
    } 
    l = d;

    d = maxd;
    i = maxi;

    if(debug) fprintf(stderr, "maxv=%i\n", maxv);

    if(maxv<0) return;

    while(i>=0 && (j=d-i)>=0) {
	if(V[d][i] == imatrix(V,d-2,i-1) + bp(left[i],right[m-1-j])) {
	    if(debug) fprintf(stderr,"traceback:\ti=%i\tj=%i\td=%i\t%c~%c\t%i\n",i,j,d,left[i],right[m-1-j],V[d][i]);
	    left[i] = toupper(left[i]);
	    right[m-1-j] = toupper(right[m-1-j]);
	    i--;j--;
	    d = i + j;
	}
	else {
	    for(k=0;k<min_helix;k++) {
		if(debug) fprintf(stderr,"traceback:\ti=%i\tj=%i\td=%i\t%c:%c\thelix\n",i,j,d,left[i],right[m-1-j]);
		left[i] = toupper(left[i]);
            	right[m-1-j] = toupper(right[m-1-j]);
            	i--;j--;
	    }
	    d = i + j;
	    if(i>=0 && j>=0) {	    
            	for(g=1;g<=max_loop;g++) {
                    for(k=0;k<=g;k++) {
                    	if(2*k<=max_loop && 2*(g-k)<=max_loop) if(W[d][i] == imatrix(V,d-g,i-k) + gap(g)) {
			    if(debug) fprintf(stderr,"traceback:\ti=%i\tj=%i\td=%i\tgap%ix%i\t%i\n",i,j,d,k,g-k,W[d][i]);
			    d = d-g;
			    i = i-k;
			    g = k = max_loop;
		    	}
                    }
            	}
		if(g==max_loop && k==max_loop) {
		    fprintf(stderr,"Something wring; can't ding the loop, exiting\n");
		    exit(1);
		}
	    }
	}
    }
    if(debug) fprintf(stderr, "traceback:\ti=%i\tj=%i\td=%i\texit\n",i,j,d);

    for(d=0;d<l;d++) {
	free(V[d]);
	free(W[d]);
    }
    free(V);
    free(W);
    free(M);
}

void sp(char *seq, int i) {
    char c;
    c = seq[i]; 
    seq[i]=0;
    printf("%s ",seq);
    seq[i]=c;
}

void fsp(char *seq, int i) {
    char c;
    c = seq[i];
    seq[i]=0;
    fprintf(stderr,"%s ",seq);
    seq[i]=c;
}


void report(char *seq, int n) {
    int i,j;
    char c;
    seq[0]   = tolower(seq[0]);
    seq[n-1] = tolower(seq[n-1]);
    for(i=0;i<n;i++) {
	if(isupper(seq[i])) {
	    sp(seq,i);
	    break;
	}
    }
    for(j=n-1;j>=0;j--) {
	if(isupper(seq[j])) {
	    break;
	}
    }
    j++;
    sp(seq+i,j-i);
    printf("%s\n",seq+j);
}

int main(int argc, char* argv[]) {
// the order in command line: sequence1 start1 end1 sequence2 start2 end2
    int i,j;
    int i1,i2,j1,j2;
    char seq1[MAXBUFFLENGTH]="";
    char seq2[MAXBUFFLENGTH]="";
    int n,m;
    char *pc;

    if(argc<=1) {
	fprintf(stderr,"Script for extention of complementary basepairs: written by dp Mar21 2012\n");
	fprintf(stderr,"Usage: %s -i seq1 i1 i2 seq2 j1 j2 -debug -helix min_helix (def=3) -loop max_loop (def=4)\n",argv[0]);
	exit(0);
    }

    for(i=1;i<argc;i++) {
        pc = argv[i];
        if(*pc == '-') {
            if(strcmp(pc+1,"i") == 0) {
		sscanf(argv[++i],"%s",&seq1[0]);
		sscanf(argv[++i],"%i",&i1);
		sscanf(argv[++i],"%i",&i2);
		sscanf(argv[++i],"%s",&seq2[0]);
		sscanf(argv[++i],"%i",&j1);
		sscanf(argv[++i],"%i",&j2);
	    }
	    if(strcmp(pc+1,"debug")==0) debug = 1;
	    if(strcmp(pc+1,"helix")==0) sscanf(argv[++i],"%i", &min_helix);
	    if(strcmp(pc+1,"loop")==0)  sscanf(argv[++i],"%i", &max_loop);
	}
    }

    n = strlen(seq1);
    m = strlen(seq2);

    if(debug) {
	fprintf(stderr,"i1=%i i2=%i j1=%i j2=%i n=%i m=%i\n",i1,i2,j1,j2,n,m);
    	fsp(seq1+i1,i2-i1);
    	fsp(seq2+j1,j2-j1);
    }

    if(0<=i1 && i1<i2 && i2<n && 0<=j1 && j1<j2 && j2<m) {
	for(i=0;i<n;i++) seq1[i] = tolower(seq1[i]);
	for(j=0;j<m;j++) seq2[j] = tolower(seq2[j]);
	i2++;j2++;
        extend(seq1 + i1, n-i1, seq2 + 0, j2);
	extend(seq2 + j1, m-j1, seq1 + 0, i2);
	report(seq1, n);
	report(seq2, m);
	exit(0);
    }
    else {
	fprintf(stderr,"Wrong input [%s %i %i %s %i %i], exiting\n",seq1, i1, i2, seq2, j1, j2);
	exit(1);
    }

}
/******************************************************************************************************************************/
int bp(char a, char b) {
    int s = (isupper(a) && isupper(b)) ? 100 : 0; 
    a = tolower(a);
    b = tolower(b);
    if(a=='a' && b == 't' || a=='t' && b == 'a') return(s + 30);
    if(a=='c' && b == 'g' || a=='g' && b == 'c') return(s + 50);
    if(a=='g' && b == 't' || a=='t' && b == 'g') return(s + 20);
    return(-INFTY);
}

int gap(int l) {
    if(l==0) return(0);
    if(l==1) return(-50);
    if(l==2) return(-65);
    if(l>=3) return(-27*l);
}
