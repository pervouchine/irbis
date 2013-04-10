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

FILE *aln_file;
FILE *net_file;
FILE *out_file;

int main(int argc, char* argv[]) {
    char aln_file_name[MAXBUFFLENGTH];
    char net_file_name[MAXBUFFLENGTH];
    char out_file_name[MAXBUFFLENGTH]="";
 
    char buff[MAXBUFFLENGTH];
    char aux[MAXBUFFLENGTH];

    int qbest, kbest, jbest;
    int kprev;

    int max_genes;
    int max_sites;

    int i,j,m,s,x,k;

    int *position;
    int *chain;
    char **mapping;

    char chr[MAXBUFFLENGTH];
    int pos, str, len, chain_id;

    int flag;

    int chridx[MAXCHR+1];
    int chroff[MAXCHR+1];

    if(argc==1) {
        fprintf(stderr,"This utility does filtering of the projected coordinates (aln) by net alignment\n");
        fprintf(stderr,"Last update by Dmitri Pervouchine (dp@crg.eu) on Mar 22, 2013\n");
        fprintf(stderr,"Keys:\n -in <aln_file> (remember to sort by position in ascending order)\n -net <net_alignment_file>\n -out <output_file> [default=stdout]\n");
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {
        if(strcmp(argv[i],"-in")==0) {
            sscanf(argv[++i], "%s", &aln_file_name[0]);
        }
        if(strcmp(argv[i],"-net")==0) {
            sscanf(argv[++i], "%s", &net_file_name[0]);
        }
        if(strcmp(argv[i],"-out")==0) {
            sscanf(argv[++i], "%s", &out_file_name[0]);
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
            fprintf(stderr,"[ERROR: output file %s cannot be opened for writing, exiting]\n", out_file_name);
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

    for(i=0;i<MAXCHR;i++) chridx[i] = chroff[i] = 0;

    if(verbose) fprintf(stderr,"[<%s, pass 1",aln_file_name);
    while(fgets(buff,MAXBUFFLENGTH,aln_file)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s" , aux);
        chridx[assign_code(aux)]++;
    }
    if(verbose) fprintf(stderr,"]\n");

    for(s=i=0;i<MAXCHR;i++) {
        x = chridx[i];
        chridx[i] =s;
        s+=x;
    }
    chridx[i] = s;

    position = (int*) malloc(sizeof(int)*(s + ARRAY_MARGIN));
    chain    = (int*) malloc(sizeof(int)*(s + ARRAY_MARGIN));
    mapping  = (char**) malloc(sizeof(char*)*(s + ARRAY_MARGIN));

    if(position==NULL || chain==NULL || mapping==NULL) {
        fprintf(stderr,"[ERROR: failed to create index tables, exiting]\n");
        exit(1);
    }

    fseek (aln_file, 0, SEEK_SET);
    if(verbose) fprintf(stderr,"[<%s, pass 2",aln_file_name);
    while(fgets(buff,MAXBUFFLENGTH,aln_file)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s" , chr);
        i = assign_code(chr);
        j = chridx[i]+chroff[i];
        sscanf(buff,"%*s %i %*i %*s %*i %*i %i" , &position[j], &chain[j]);
	mapping[j] = (char*)malloc(sizeof(char)*(strlen(buff) + ARRAY_MARGIN));
	if(mapping[j]==NULL) {
	    fprintf(stderr,"[ERROR: not enough memory, exiting]\n");
	    exit(1);
	}
	strcpy(mapping[j], buff);
        chroff[i]++;
    }
    fclose(aln_file);
    if(verbose) fprintf(stderr,"]\n");

    if(verbose) fprintf(stderr,"[Sort by position (if not done before)");
    for(i=0;i<MAXCHR;i++) {
        flag=1;
        while(flag) {
            flag=0;
            for(j=chridx[i];j<chridx[i+1]-1;j++) {
                if(position[j]>position[j+1]) {
                    flag=1;
                    swapi(position+j,position+j+1);
		    swapi(chain+j,chain+j+1);
		    swaps(mapping+j,mapping+j+1);
                }
            }
        }
    }
    if(verbose) fprintf(stderr,"]\n");

    net_file= fopen(net_file_name,"r");
    if(net_file==NULL) {
	fprintf(stderr,"[ERROR: cannot access %s, exiting]\n", net_file_name);
        exit(1);
    }


    if(verbose) fprintf(stderr,"[<%s",net_file_name);
    i =-1;
    while(fgets(buff,MAXBUFFLENGTH,net_file)) {
        if(strlen(buff)<2) break;
	sscanf(buff,"%s" , &aux[0]);
	if(strcmp(aux,(char*)"net")==0) {
	    sscanf(buff,"%*s %s" , &chr[0]);
	    i = assign_code(chr);
	}
	if(strcmp(aux,(char*)"fill")==0) {
	    sscanf(buff,"%*s %i %i %*s %*s %*s %*s %*s %i\n", &pos, &len, &chain_id);
	    if(i<0) continue;
	    for(j=chridx[i];j<chridx[i+1];j++) {
		if(pos<=position[j] && position[j]<pos + len -1 && chain[j]==chain_id) {
		    fprintf(out_file, "%s", mapping[j]); 
		    for(k=j;position[k]==position[j];k--) chain[k]=-1;
		    for(k=j;position[k]==position[j];k++) chain[k]=-1;
		}
	    }
	}
    }
    fclose(net_file);
    if(verbose) fprintf(stderr,"]\n");

    for(i=0;i<MAXCHR;i++) {
	for(j=chridx[i];j<chridx[i+1];j++) {
	    free(mapping[j]);
	}
    }

    free(position);
    free(chain);
    free(mapping);
    timestamp_report();

}
