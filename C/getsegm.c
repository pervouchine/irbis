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
#include <sys/dir.h>
#include "genutils.h"

int main(int argc, char* argv[]) {
    char out_file_name[MAXBUFFLENGTH]="";
    char bed_file_name[MAXBUFFLENGTH]="";
    char idx_file_name[MAXBUFFLENGTH];
    char dbx_file_name[MAXBUFFLENGTH];

    FILE *idx_file;
    FILE *dbx_file;
    FILE *bed_file;
    FILE *outfile;


    long offset;
    long seqlen;

    long length_limit = MAXLONGBUFFLENGTH;

    long half_length_limit;
    char name[MAXBUFFLENGTH];

    char chr_name[MAXBUFFLENGTH];
    char buff[MAXBUFFLENGTH];
    char longbuff[MAXLONGBUFFLENGTH + 100];

    int output_type=OUTPUT_TBMAF;

    int b,i,j,q,k,n,a,m,s;
    char c;
    long x,y;

    long** beg;
    long** end;
    int** str;
    char*** ids;
    long pos,p,l;
    long margin=0;
    int warnings=0;
    int cut=0;

    int record_count[MAXCHR];
    int record_idx[MAXCHR];

    if(argc==1) {
	fprintf(stderr,"This routine get sequence segments from a custom compressed FASTA repository  (see transf)\n");
        fprintf(stderr,"Last update by Dmitri Pervouchine (dp@crg.eu) on Mar 22, 2013\n");
        fprintf(stderr," -in <bed_file>\n -dbx <database_file>\n -idx <database_index_file>n -out <output_file>\n");
        fprintf(stderr," -quiet suppress verbose output [default=NO]\n -limit <length_limit> [default=%i]\n -margin <margin_length> [default=%i]\n", length_limit, margin);
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {
        if(strcmp(argv[i],"-in")==0) {
            sscanf(argv[++i], "%s", &bed_file_name[0]);
        }

        if(strcmp(argv[i],"-dbx")==0) {
            sscanf(argv[++i], "%s", &dbx_file_name[0]);
        }

        if(strcmp(argv[i],"-idx")==0) {
            sscanf(argv[++i], "%s", &idx_file_name[0]);
        }

        if(strcmp(argv[i],"-out")==0) {
            sscanf(argv[++i], "%s", &out_file_name[0]);
        }

        if(strcmp(argv[i],"-limit")==0) {
            sscanf(argv[++i], "%li", &length_limit);
        }

        if(strcmp(argv[i],"-margin")==0) {
            sscanf(argv[++i], "%li", &margin);
        }

        if(strcmp(argv[i],"-quiet")==0) {
            verbose = 0;
        }
    }


    half_length_limit = length_limit/2;


    if(out_file_name[0]==0) {
        fprintf(stderr,"[WARNING: output file not specified, redirect to stdout]\n");
        outfile = stdout;
    }
    else {
        outfile = fopen(out_file_name,"w");
        if(outfile == NULL) {
            fprintf(stderr,"[ERROR: output file %s cannot be opened for writing, exiting]\n", out_file_name);
            exit(1);
        }
	if(verbose) fprintf(stderr,"[>%s]\n",out_file_name);
    }


    bed_file = fopen(bed_file_name,"r");
    if(bed_file == NULL) {
        fprintf(stderr,"[ERROR: cannot access %s, exiting]\n", bed_file_name);
        exit(1);
    }

    if(verbose) fprintf(stderr,"[<%s, pass 1", bed_file_name);
    while(fgets(buff,MAXBUFFLENGTH,bed_file)) {
        if(strlen(buff)<2) break;
      	sscanf(buff,"%s" , &chr_name[0]);
      	n = assign_code(chr_name);
      	record_count[n]++;
    }
    if(verbose) fprintf(stderr,"]\n");

    beg = (long**) malloc(sizeof(long*)*(N_CHR_NAMES+1));
    end = (long**) malloc(sizeof(long*)*(N_CHR_NAMES+1));
    str = (int**) malloc(sizeof(int*)*(N_CHR_NAMES+1));
    ids = (char***)malloc(sizeof(char**)*(N_CHR_NAMES+1));

    if(beg==NULL || end ==NULL || str==NULL || ids==NULL) {
        fprintf(stderr,"[ERROR: failed to create index tables, exiting]\n");
        exit(1);
    }

    for(i=0;i<N_CHR_NAMES;i++) {
    	if(record_count[i]>0) {
	    beg[i] = (long*) malloc(sizeof(long)*(record_count[i]+1));
	    end[i] = (long*) malloc(sizeof(long)*(record_count[i]+1));
	    str[i] = (int*) malloc(sizeof(int)*(record_count[i]+1));
	    ids[i] = (char**) malloc(sizeof(char*)*(record_count[i]+1));
	    record_idx[i]=0;
	    if(beg[i]==NULL || end[i]==NULL || str[i]==NULL || ids[i]==NULL) {
        	fprintf(stderr,"[ERROR: failed to create index tables, exiting]\n");
        	exit(1);
    	    }
	}
    }

    fseek(bed_file,  0, SEEK_SET);
    if(verbose) fprintf(stderr,"[<%s, pass 2", bed_file_name);
    while(fgets(buff,MAXBUFFLENGTH,bed_file)) {
        if(strlen(buff)<2) break;
	sscanf(buff,"%s" , &chr_name[0]);
	i = get_chr_code(chr_name);
	j = record_idx[i];
	sscanf(buff,"%*s %li %li %*s %s %i" ,&beg[i][j], &end[i][j], &name[0], &str[i][j]);
        if(beg[i][j]>end[i][j]) {
            pos=beg[i][j];beg[i][j]=end[i][j];end[i][j]=pos;
        }
	ids[i][j] = (char*)malloc(sizeof(char)*(strlen(name)+1));
	strcpy(ids[i][j],name);
	record_idx[i]++;
    }
    fclose(bed_file);
    if(verbose) fprintf(stderr,"]\n");

    if(verbose) fprintf(stderr,"[<%s,%s",idx_file_name,dbx_file_name);
    idx_file = fopen(idx_file_name,"r");
    dbx_file = fopen(dbx_file_name,"r");
    if(idx_file == NULL || dbx_file == NULL) {
        fprintf(stderr,"[ERROR: cannot access %s or %s, exiting]\n", idx_file_name, dbx_file_name);
        exit(1);
    }

    offset = 0;
    while(fgets(buff,MAXBUFFLENGTH,idx_file)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s" , &name[0]);
	while(fgets(buff,MAXBUFFLENGTH,idx_file)) {
            if(strlen(buff)<2) break;
	    sscanf(buff,"%s %li" , &chr_name[0], &seqlen);
	    i = get_chr_code(chr_name);
	    for(k=0;k<record_count[i];k++) {
		if(end[i][k]>seqlen) {
		    warnings++;
		    continue;
		}
		l = end[i][k]-beg[i][k];
		if(l==0) continue;
		beg[i][k]-= margin;
		end[i][k]+= margin;
		l+=2*margin;
		p = beg[i][k] - 1 + (str[i][k]<0 ? 1 : 0);
		if(l<length_limit) {
		    fget_segment(longbuff, dbx_file, offset, p, l);
		}
		else {
		    fget_segment(longbuff, dbx_file, offset, p, half_length_limit);
		    strcat(longbuff + half_length_limit,".....");
                    fget_segment(longbuff+half_length_limit+5, dbx_file, offset, p+l-half_length_limit, half_length_limit);
		    cut++;
		}
		if(str[i][k]<0) {
		    rev1(longbuff);
		}
		if(is_all_n(longbuff)) continue;
		m = strlen(longbuff)-1;
		if(output_type==OUTPUT_FASTA) 
		    fprintf(outfile,">%s\n%s\n",ids[i][k], longbuff);
                if(output_type==OUTPUT_TABSQ) 
		    fprintf(outfile,"%s\t%s\n",ids[i][k], longbuff);
		if(output_type==OUTPUT_TBMAF) 
		    fprintf(outfile,"%s\t%s\t%li\t%li\t%i\t%li\t%s\n",ids[i][k], chr_name, (str[i][k]>0 ? beg[i][k] : seqlen - end[i][k]), l, str[i][k], seqlen, longbuff);
	    }
            offset+= (seqlen % 8 == 0) ? seqlen/8 : (seqlen/8 + 1);
	}
    }
    fclose(outfile);
    fclose(idx_file);
    fclose(dbx_file);
    if(verbose) fprintf(stderr,"]\n");
    if(warnings) fprintf(stderr,"[WARNING: %i segments out of range]\n", warnings);
    if(warnings) fprintf(stderr,"[WARNING: %i segments shortened due to length limit]\n", cut);
    timestamp_report();
    exit(0);
}
