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

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/dir.h>
#include "genutils.h"

int main(int argc, char** argv) {
    char out_file_name[MAXBUFFLENGTH]="";
    char cps_file_name[MAXBUFFLENGTH]="";
    char chr[MAXBUFFLENGTH]; 
    char* file_name[MAXSPECIES];

    char buff[MAXLONGBUFFLENGTH + MAXBUFFLENGTH];
    char seq[MAXLONGBUFFLENGTH];

    int *key;

    FILE *inp_file;
    FILE *out_file;
    FILE *cps_file;

    int	x, max_key;
    int	i, j, nk, n_species;

    long pos, len, total;
    int strand;
    int *ind;
    char **muf;

    if(argc==1) {
        fprintf(stderr,"This utility creates multiple unaligned file (muf) from sequence files sus/suw based on a cps file\n");
        fprintf(stderr,"Last update by Dmitri Pervouchine (dp@crg.eu) on Mar 22, 2013\n");
	fprintf(stderr, "Usage: %s -in <cps_file> -files <input_file_1> <input_file_2> ... <input_file_n> -o <output_file>\n",argv[0]);
	exit(1);
    }

    timestamp_set();
    n_species=0;
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i],"-in")==0) {
            sscanf(argv[++i], "%s", &cps_file_name[0]);
        }
        if(strcmp(argv[i],"-out")==0) {
            sscanf(argv[++i], "%s", &out_file_name[0]);
        }
	if(strcmp(argv[i],"-files")==0) {
	    for(;i+1<argc;i++) {
		if(argv[i+1][0] == '-') break;
		file_name[n_species]    = (char*)malloc(sizeof(char)*MAXBUFFLENGTH);
		strcpy(file_name[n_species++], argv[i+1]);
	    }
	}
	if(strcmp(argv[i],"-quiet")==0) {
	    verbose = 0;
	}
    }

    if(out_file_name[0]==0) {
        fprintf(stderr,"[WARNING: output file not specified, redirect to stdout]\n");
        out_file = stdout;
    }
    else {
        out_file = fopen(out_file_name,"w");
        if(out_file == NULL) {
            fprintf(stderr,"[ERROR: output file %s cannot be opened for writing, exiting]\n", out_file_name);
            exit(1);
        }
        if(verbose) fprintf(stderr,"[>%s]\n",out_file_name);
    }

    cps_file= fopen(cps_file_name,"r");
    if(cps_file==NULL) {
        fprintf(stderr,"Can't access CPS file. Exiting\n");
        exit(1);
    }

    if(verbose) fprintf(stderr,"[<%s",cps_file_name);
    max_key=0;
    nk=0;
    while(fgets(buff,MAXBUFFLENGTH,cps_file)) {
        sscanf(buff,"%*s %*i %*i %*i %i",&x);
        if(x>max_key) max_key=x;
	nk++;
    }

    fseek (cps_file, 0, SEEK_SET);

    key = (int*)malloc(sizeof(int)*(nk + ARRAY_MARGIN));

    if(key ==NULL) {
	fprintf(stderr,"[ERROR: too many keys, exiting]\n");
        exit(1);
    }

    fseek (cps_file, 0, SEEK_SET);
    nk = 0;
    while(fgets(buff,MAXBUFFLENGTH,cps_file)) {
        sscanf(buff,"%*s %*i %*i %*i %i",&x);
        key[nk++]=x;
    }
    fclose(cps_file);
    if(verbose) fprintf(stderr,"]\n");

    ind = (int*)malloc(sizeof(int)*(max_key + ARRAY_MARGIN));
    muf = (char**)malloc(sizeof(char*)*(max_key + ARRAY_MARGIN));

    if(ind==NULL || muf==NULL) {
	fprintf(stderr,"[ERROR: failed to create index, exiting]\n");
        exit(1);
    }

    for(i=0;i<=max_key;i++) ind[i]=0;
    for(j=0;j<nk;j++) {
	ind[key[j]]=1;
	muf[key[j]]=(char*)malloc(sizeof(char)*MAXLONGBUFFLENGTH*MAXSPECIES);
	if(muf[key[j]]==NULL) {
            fprintf(stderr,"[ERROR: failed to allocate memory, exiting]\n");
	    exit(1);
	}
	muf[key[j]][0]=0;
	if(verbose) fprintf(stderr,"[%i]",key[j]);
    }
    if(verbose) fprintf(stderr,"\n");

    for(i=0;i<n_species;i++) {
	inp_file = fopen(file_name[i], "r");
	if(verbose) fprintf(stderr,"[%s, ", file_name[i]);
	if(inp_file==NULL) {
       	    fprintf(stderr,"[ERROR: cannot access %s, exiting]\n", file_name[i]);
            exit(1);
        }
	for(j=strlen(file_name[i])-1;file_name[i][j]!='.' && j>=0;j--);
	file_name[i][j]=0;
	for(j--;file_name[i][j]!='/' && j>=0;j--);
        if(verbose) fprintf(stderr,"tag=%s", file_name[i] + j + 1);
	while(fgets(buff, MAXLONGBUFFLENGTH + MAXBUFFLENGTH, inp_file)) {
	    x=0;
	    sscanf(buff,"%i %s %li %li %i %li %s", &x, &chr[0], &pos, &len, &strand, &total, &seq[0]);
	    if(x==0) {fprintf(stderr,"!");continue; }
	    if(x<=max_key) {
		if(ind[x]) {
                    sprintf(buff,"s %s.%s %li %li %i %li %s\n", file_name[i] + j + 1, chr, pos, len, strand, total, seq);	
		    strcat(muf[x],buff);
		}
 	    }
	}
	if(verbose) fprintf(stderr,"]\n");
    }

    for(i=0;i<=max_key;i++) {
	if(ind[i]) if(strlen(muf[i])>1) fprintf(out_file,"a id=%i\n%s\n", i, muf[i]);	
    }

    fclose(out_file);

    timestamp_report();
    exit(0);
}
