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
#include <math.h>
#include <string.h>

#include "genutils.h"
#include "progressbar.h"

#define MAXALN 2000000

FILE *aln_file;
FILE *out_file;

int main(int argc, char* argv[]) {
    char aln_file_name[MAXBUFFLENGTH];
    char out_file_name[MAXBUFFLENGTH]="";
 
    char buff[MAXBUFFLENGTH];
    char chr1[MAXBUFFLENGTH];
    char chr2[MAXBUFFLENGTH];

    double dthreshold = 1.50;
    int dlimit = 100000;
    int max_depth = 4;

    int** chr_t;
    int** str_t;
    int** pos_t;

    int* chr_q;
    int* str_q;
    int* pos_q;

    int **score;
    int **jbest;
    int **lbest;

    int a,b,d,dmin,lmin,s;

    int *count;
    int *ptr;

    int max_rec=0;

    int pos1, pos2, str1, str2;
    int i, j, k, l, n;

    int s_max, k_max;

    if(argc==1) {
        fprintf(stderr,"This utility does ad-hoc filtering of the projected coordinates (aln) by maximum synteny\n");
        fprintf(stderr,"Last update by Dmitri Pervouchine (dp@crg.eu) on Mar 26, 2013\n");
        fprintf(stderr,"Keys:\n -in <aln_file> (remember to sort by position in ascending order)\n -out <output_file> [default=stdout]\n -maxdepth, -threshold, -lendiff params (see code)\n");
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {
        if(strcmp(argv[i],"-in")==0) {
            sscanf(argv[++i], "%s", &aln_file_name[0]);
        }
        if(strcmp(argv[i],"-out")==0) {
            sscanf(argv[++i], "%s", &out_file_name[0]);
        }
        if(strcmp(argv[i],"-lendiff")==0) {
            sscanf(argv[++i], "%i", &dlimit);
        }
        if(strcmp(argv[i],"-threshold")==0) {
            sscanf(argv[++i], "%lf", &dthreshold);
        }
        if(strcmp(argv[i],"-maxdepth")==0) {
            sscanf(argv[++i], "%i", &max_depth);
        }
        if(strcmp(argv[i],"-quiet")==0) {
            verbose=0;
        }
    }

    if(out_file_name[0]==0) {
        fprintf(stderr,"[WARNING: output file not specified, redirect to stdout]\n");
        out_file = stdout;
    }
    else {
        out_file = fopen(out_file_name,"w");
        if(out_file == NULL) {
            fprintf(stderr,"[ERROR: output file (%s) cannot be opened, exiting]\n", out_file_name);
            exit(1);
        }
	if(verbose) fprintf(stderr,"[>%s]\n",out_file_name);
    }

/*******************************************************************************************************/
    aln_file= fopen(aln_file_name,"r");
    if(aln_file==NULL) {
        fprintf(stderr,"[ERROR: cannot access %s, exiting]\n", aln_file_name);
	exit(1);
    }

    if(verbose) fprintf(stderr,"[<%s, pass 1", aln_file_name);
    while(fgets(buff, MAXBUFFLENGTH, aln_file)) {
	if(strlen(buff)<2) break;
	max_rec++;
    }
    if(verbose) fprintf(stderr,"]\n");

    chr_q = (int*)malloc(sizeof(int)*(max_rec + ARRAY_MARGIN));
    str_q = (int*)malloc(sizeof(int)*(max_rec + ARRAY_MARGIN));
    pos_q = (int*)malloc(sizeof(int)*(max_rec + ARRAY_MARGIN));

    chr_t = (int**)malloc(sizeof(int*)*(max_rec + ARRAY_MARGIN));
    str_t = (int**)malloc(sizeof(int*)*(max_rec + ARRAY_MARGIN));
    pos_t = (int**)malloc(sizeof(int*)*(max_rec + ARRAY_MARGIN));

    score = (int**)malloc(sizeof(int*)*(max_rec + ARRAY_MARGIN));
    lbest = (int**)malloc(sizeof(int*)*(max_rec + ARRAY_MARGIN));
    jbest = (int**)malloc(sizeof(int*)*(max_rec + ARRAY_MARGIN));

    count = (int*)malloc(sizeof(int)*(max_rec + ARRAY_MARGIN));
    ptr   = (int*)malloc(sizeof(int)*(max_rec + ARRAY_MARGIN));

    if(chr_q == NULL || str_q == NULL || pos_q == NULL || count == NULL || ptr == NULL || chr_t == NULL || str_t == NULL || pos_t == NULL) {
	fprintf(stderr,"[ERROR: not enough memory, exiting]\n");
	exit(1);
    }

    for(i=0;i<max_rec;i++) count[i] = ptr[i] = 0;

    if(verbose) fprintf(stderr,"[<%s, pass 2", aln_file_name);
    fseek (aln_file, 0, SEEK_SET);
    n=0;
    while(fgets(buff, MAXBUFFLENGTH, aln_file)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s %i %i" , &chr1[0], &pos1, &str1);
	if(assign_code(chr1) != chr_q[n] || str1 != str_q[n] || pos1 != pos_q[n]) {
	    n++;
	    chr_q[n] = assign_code(chr1);
            str_q[n] = str1;
            pos_q[n] = pos1;
	}
	count[n]++;
    }
    if(verbose) fprintf(stderr,"]\n");

    for(i=1; i<=n; i++) {
        chr_t[i] = (int*)malloc(sizeof(int)*(count[i] + ARRAY_MARGIN));
	pos_t[i] = (int*)malloc(sizeof(int)*(count[i] + ARRAY_MARGIN));
        str_t[i] = (int*)malloc(sizeof(int)*(count[i] + ARRAY_MARGIN));

	score[i] = (int*)malloc(sizeof(int)*(count[i] + ARRAY_MARGIN));
	lbest[i] = (int*)malloc(sizeof(int)*(count[i] + ARRAY_MARGIN));
        jbest[i] = (int*)malloc(sizeof(int)*(count[i] + ARRAY_MARGIN));

	if(chr_t[i] == NULL || str_t[i] == NULL || pos_t[i] == NULL) {
	    fprintf(stderr,"[ERROR: not enough memory, exiting]\n");
	    exit(1);
	}
    }

    if(verbose) fprintf(stderr,"[<%s, pass 3", aln_file_name);
    fseek (aln_file, 0, SEEK_SET);
    n=0;
    while(fgets(buff, MAXBUFFLENGTH, aln_file)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s %i %i %s %i %i" , &chr1[0], &pos1, &str1, &chr2[0], &pos2, &str2);
	if(assign_code(chr1) != chr_q[n] || str1 != str_q[n] || pos1 != pos_q[n]) n++;
	chr_t[n][ptr[n]] = assign_code(chr2);
	pos_t[n][ptr[n]] = pos2;
	str_t[n][ptr[n]] = str2;
	ptr[n]++;
    }
    if(verbose) fprintf(stderr,"]\n");

    for(i=1; i<=n; i++) {
	progressbar(i, n, (char*)"Filtering ", verbose);
	for(k=0; k<count[i]; k++) {
	    score[i][k] = 0;
	    lbest[i][k] = -1;
	    jbest[i][k] = -1;
	    for(j=i-1; j>0 && i-j<=max_depth; j--) {
		a = abs(pos_q[i] - pos_q[j]);
                dmin = INFTY;
                lmin = -1;
		for(l=0; l<count[j]; l++) {
		    if(chr_t[i][k] == chr_t[j][l] && str_t[i][k]*str_t[j][l] == str_q[i]*str_q[j]) {
			b = abs(pos_t[i][k] - pos_t[j][l]);
                        d = abs(b-a);
                        if(d<dmin) {
                            dmin = d;
                            lmin = l;
                        }
                    }
                }
		s = (lmin>=0 && (((double)dmin/a) < dthreshold || dmin < dlimit)) ? score[j][lmin] + a : 0;
		if(s > score[i][k]) {
		    score[i][k] = s;
		    jbest[i][k] = j;
		    lbest[i][k] = lmin;
		} 
	    }
	}
    }
    
    for(i=n;i>0;i--) {
	s_max = 0;
	k_max = -1;
	for(k=0; k<count[i]; k++) {
	    if(score[i][k]>s_max) {
		s_max = score[i][k];
		k_max = k;
	    }
	}
	if(k_max>=0) {
	    k = k_max; 
	    while(jbest[i][k]>=0 && lbest[i][k]>=0) {
		fprintf(out_file,"%s\t%i\t%i\t", get_chr_name(chr_q[i]), pos_q[i], str_q[i]);
		fprintf(out_file,"%s\t%i\t%i\n", get_chr_name(chr_t[i][k]), pos_t[i][k], str_t[i][k]);
		j = jbest[i][k];
		l = lbest[i][k];
		i = j;
		k = l;
	    }
	}

    }
    timestamp_report();
    exit(0);
}
