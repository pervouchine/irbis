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

FILE *cps_file;
FILE *out_file;
FILE *chain_file;

int main(int argc, char* argv[]) {
    char cps_file_name[MAXBUFFLENGTH];
    char chain_file_name[MAXBUFFLENGTH];
    char out_file_name[MAXBUFFLENGTH]="";
 
    int marginlength = 0;

    char buff[MAXBUFFLENGTH+1];
    char aux[MAXBUFFLENGTH+1];

    long score;
    int start1,end1,len1,start2,end2,len2;
    char strand1, strand2, chr1[MAXBUFFLENGTH], chr2[MAXBUFFLENGTH];

    int *size, *dq, *dt;
    int a,b,k,i,j,s,m,x;

    int *position;
    int *strand;

    int chridx[MAXCHR+1];
    int chroff[MAXCHR+1];

    char resstr;
    int  rescrd;
    long chain_id;


    if(argc==1) {
	fprintf(stderr,"This utility does liftOver of coordinates (cps) by  using chain alignment\n");
        fprintf(stderr,"Last update by Dmitri Pervouchine (dp@crg.eu) on Mar 22, 2013\n");
	fprintf(stderr,"Keys:\n -in <cps_file> (remember to sort by position in ascending order)\n -chain <chain_alignment_file>\n -out <output_file> [default=stdout]\n");
 	fprintf(stderr," -margin margin length [default=0]\n -quiet suppress verbose output [default=NO]\n");
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {
	if(strcmp(argv[i],"-in")==0) {
	   sscanf(argv[++i], "%s", &cps_file_name[0]);
	}
	if(strcmp(argv[i],"-chain")==0) {
	   sscanf(argv[++i], "%s", &chain_file_name[0]);
	}
	if(strcmp(argv[i],"-out")==0) {
           sscanf(argv[++i], "%s", &out_file_name[0]);
        }
	if(strcmp(argv[i],"-margin")==0) {
           sscanf(argv[++i], "%i", &marginlength);
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

    cps_file= fopen(cps_file_name,"r");
    if(cps_file==NULL) {
	fprintf(stderr,"[ERROR: cannot access %s, exiting]\n", cps_file_name);
	exit(1);
    }

    for(i=0;i<MAXCHR;i++) chridx[i] = chroff[i] = 0;

    if(verbose) fprintf(stderr,"[<%s, pass 1",cps_file_name);
    while(fgets(buff,MAXBUFFLENGTH,cps_file)) {
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

    position   = (int*) malloc(sizeof(int)*(s + ARRAY_MARGIN));
    strand     = (int*) malloc(sizeof(int)*(s + ARRAY_MARGIN));

    if(position==NULL || strand==NULL) {
        fprintf(stderr,"[ERROR: failed to create index tables, exiting]\n");
        exit(1);
    }

    fseek (cps_file, 0, SEEK_SET);
    if(verbose) fprintf(stderr,"[<%s, pass 2", cps_file_name);
    while(fgets(buff,MAXBUFFLENGTH,cps_file)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s" , aux);
        i = assign_code(aux);
        m = chridx[i]+chroff[i];
        sscanf(buff,"%*s %i %i" , &position[m], &strand[m]);
        chroff[i]++;
    }
    fclose(cps_file);
    if(verbose) fprintf(stderr,"]\n");


    if(verbose) fprintf(stderr,"[Sort by position (if not done before)");
    for(i=0;i<MAXCHR;i++) {
        k=1;
        while(k) {
            k=0;
            for(j=chridx[i];j<chridx[i+1]-1;j++) {
                if(position[j]>position[j+1]) {
                    k=1;
                    swapi(position+j,position+j+1);
                    swapi(strand+j,strand+j+1);
                }
            }
        }
    }
    if(verbose) fprintf(stderr,"]\n");


/**********************************************************************************************/
    size = (int*) malloc(sizeof(int)*(MAXALN + ARRAY_MARGIN));
    dq   = (int*) malloc(sizeof(int)*(MAXALN + ARRAY_MARGIN));
    dt   = (int*) malloc(sizeof(int)*(MAXALN + ARRAY_MARGIN));

    if(size ==0 || dq ==0 || dt==0) {
        fprintf(stderr,"[ERROR: not enough memory for chains, exiting]\n");
        exit(1);
    }


/**********************************************************************************************/

    chain_file = fopen(chain_file_name,"r");
    if(chain_file==NULL) {
	fprintf(stderr,"[ERROR: cannot access %s, exiting]\n", chain_file_name);
	exit(1);
    }

    fseek(chain_file, 0, SEEK_END);
    unsigned int last_pos = ftell(chain_file);
    fseek(chain_file, 0, SEEK_SET);

    while(fgets(buff,MAXBUFFLENGTH,chain_file)) {
        if(strlen(buff)<2) break;
     	buff[5]=0;
     	if(strcmp(buff,"chain")==0) {
       	    sscanf(buff+6,"%li %s %i %c %i %i %s %i %c %i %i %li",&score, &chr1[0], &len1, &strand1, &start1, &end1, &chr2[0], &len2, &strand2, &start2, &end2, &chain_id);
	    k=0;
	    while(fgets(buff,MAXBUFFLENGTH,chain_file)) {
		if(strlen(buff)<2) break;
		progressbar(ftell(chain_file), last_pos-1, (char*)"Processing ", verbose);
	    	sscanf(buff,"%i %i %i",size + k, dt + k, dq + k);
	    	k++;
	    	if(k>=MAXALN) {
		    fprintf(stderr,"[ERROR: chain too long, exiting]\n");
		    exit(1);
	    	}
	    }

	    x = get_chr_code(chr1);
	    if(x<0) continue;

            a=start1;b=start2;
            j=0;

	    for(i=chridx[x];i<chridx[x+1] && position[i]<start1;i++);
	    for(;i<chridx[x+1]&& position[i]<end1;i++) {
	    	while(position[i]>a+size[j]+dt[j] && j<k){
		    a+=size[j]+dt[j];
		    b+=size[j]+dq[j];
		    j++;
	    	}
	        if(j>=k) break;
	        if(position[i]-a > marginlength && a+size[j]-position[i] >= marginlength) {
                    if(strand1==strand2) {
                    	resstr = strand[i];
                    	rescrd = position[i] - a + b;
                    }
                    else {
                    	resstr = -strand[i];
                    	rescrd = len2 - (position[i] - a + b - 1) ;
                    }
		    fprintf(out_file,"%s\t%i\t%i\t%s\t%i\t%i\t%li\n",chr1, position[i], strand[i], chr2, rescrd, resstr, chain_id);
	    	}
	    }
     	}
    }
    fclose(chain_file);
    fclose(out_file);
    timestamp_report();

    free(size);
    free(dq);
    free(dt);

    free(position);
    free(strand);
    exit(0);
}
