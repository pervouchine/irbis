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
    char input_file_name[MAXBUFFLENGTH]="";
    char idx_file_name[MAXBUFFLENGTH];
    char dbx_file_name[MAXBUFFLENGTH];

    long offset;
    long seqlen;
    long intronic_window = 150;
    long exonic_window = 0;

    char filename[MAXBUFFLENGTH];
    char name[MAXBUFFLENGTH];

    char chr_name[MAXBUFFLENGTH];
    char buff[MAXBUFFLENGTH];
    char longbuff[MAXLONGBUFFLENGTH];
    char longbuffm[MAXLONGBUFFLENGTH];

    FILE *idx_file;
    FILE *dbx_file;
    FILE *input_file;
    FILE *outfile;
   
    int b,i,j,q,k,n,a,m,s;
    char c;
    long position,y;
    int strand;

    long** pos;
    int** str;
    char** typ;
    char*** ids;
    long p,l;

    int record_count[MAXCHR];
    int record_idx[MAXCHR];
    int cis = 0;
    int coord = 0;
    int all=0;

    int warnings = 0;

    char format[][64] = {"%*s %*i %*i %s %li %i %*s %s %c", "%s %li %i %*s %s %c"};

    if(argc==1) {
        fprintf(stderr,"This routine get sequence segments from a custom compressed FASTA repository  (see transf)\n");
        fprintf(stderr,"Last update by Dmitri Pervouchine (dp@crg.eu) on Mar 22, 2013\n");
        fprintf(stderr," -in <aln_file>\n -dbx <database_file>\n -idx <index_file>\n -out <output_file>\n");
        fprintf(stderr," -we <exonic_window> [default=%i]\n -wi <intronic_window> [default=%i]\n -cis [use colunms 1-3] [default=%i]\n", exonic_window, intronic_window, cis);
	fprintf(stderr," -quiet <suppress verbose output> [default=no]\n -all <include all sites>\n -coord <offset for 3'-sites> [default=%i]\n",coord);
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {
        if(strcmp(argv[i],"-in")==0) {
            sscanf(argv[++i], "%s", &input_file_name[0]);
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

        if(strcmp(argv[i],"-we")==0) {
            sscanf(argv[++i], "%li", &exonic_window);
        }

        if(strcmp(argv[i],"-wi")==0) {
            sscanf(argv[++i], "%li", &intronic_window);
        }

        if(strcmp(argv[i],"-coord")==0) {
            sscanf(argv[++i], "%i", &coord);
        }

        if(strcmp(argv[i],"-quiet")==0) {
	    verbose = 0;
	}

        if(strcmp(argv[i],"-cis")==0) {
            cis = 1;
        }

        if(strcmp(argv[i],"-all")==0) {
            all = 1;
        }

    }

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

    input_file = fopen(input_file_name,"r");
    if(input_file == NULL) {
        fprintf(stderr,"[ERROR: cannot access %s, exiting]\n", input_file_name);
        exit(1);
    }

    if(verbose) fprintf(stderr,"[<%s, pass 1",input_file_name);
    while(fgets(buff,MAXBUFFLENGTH,input_file)) {
        if(strlen(buff)<2) break;
      	sscanf(buff, format[cis], &chr_name[0], &position, &strand, &name[0], &c);
      	n = assign_code(chr_name);
      	record_count[n]++;
    }
    if(verbose) fprintf(stderr,"]\n");

    pos = (long**) malloc(sizeof(long*)*(N_CHR_NAMES+1));
    str = (int**)  malloc(sizeof(int*)*(N_CHR_NAMES+1));
    typ = (char**) malloc(sizeof(char*)*(N_CHR_NAMES+1));
    ids = (char***)malloc(sizeof(char**)*(N_CHR_NAMES+1));

    if(pos==NULL || str==NULL || typ==NULL || ids==NULL) {
        fprintf(stderr,"[ERROR: failed to create index tables, exiting]\n");
        exit(1);
    }

    for(i=0;i<N_CHR_NAMES;i++) {
    	if(record_count[i]>0) {
	    pos[i] = (long*)  malloc(sizeof(long)*(record_count[i]+1));
	    str[i] = (int*)   malloc(sizeof(int)*(record_count[i]+1));
            typ[i] = (char*)  malloc(sizeof(char)*(record_count[i]+1));
	    ids[i] = (char**) malloc(sizeof(char*)*(record_count[i]+1));
	    if(pos[i]==NULL || str[i]==NULL || typ[i]==NULL || ids[i]==NULL) {
        	fprintf(stderr,"[ERROR: failed to create index tables, exiting]\n");
        	exit(1);
    	    }
	    record_idx[i]=0;
	}
    }

    if(verbose) fprintf(stderr,"[<%s, pass 2",input_file_name);
    fseek(input_file,  0, SEEK_SET);
    while(fgets(buff,MAXBUFFLENGTH,input_file)) {
        if(strlen(buff)<2) break;
	sscanf(buff, format[cis], &chr_name[0], &position, &strand, &name[0], &c);
	i = get_chr_code(chr_name);
	j = record_idx[i];
	pos[i][j] = position + (c=='D' ? coord*strand : 0);
	str[i][j] = strand;
	typ[i][j] = c;
	ids[i][j] = (char*) malloc(sizeof(char)*(strlen(name)+1));
	strcpy(ids[i][j],name);
	record_idx[i]++;
    }
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
               if(pos[i][k]>seqlen) {
		    warnings++;
                    continue;
                }
		l = exonic_window + intronic_window;
		if(typ[i][k]=='D' || typ[i][k]=='A' || all) { 
		    if(str[i][k]>0) {
		    	p = pos[i][k] - 1 - (typ[i][k]=='D' ? exonic_window : intronic_window);
		    }
		    else {
			p = pos[i][k] - (typ[i][k]=='D' ? intronic_window : exonic_window);
		    }
		    fget_segment(longbuff, dbx_file, offset, p, l);
		    if(str[i][k]<0) {
		    	rev1(longbuff);
		    }
               	    if(is_all_n(longbuff)) continue;
		    fprintf(outfile,"%s\t%s\t%li\t%li\t%i\t%li\t%s\t%c\n",ids[i][k], chr_name, (str[i][k]>0 ? p + 1 : seqlen - (p + l)), l, str[i][k], seqlen, longbuff,typ[i][k]);
		}
	    }
            offset+= (seqlen % 8 == 0) ? seqlen/8 : (seqlen/8 + 1);
	}
    }
    fclose(outfile);
    fclose(idx_file);
    fclose(dbx_file);
    if(verbose) fprintf(stderr,"]\n");
    if(verbose && warnings>0) fprintf(stderr,"[WARNING: %i windows were out of range, they were ignored]\n", warnings);

    timestamp_report();
    exit(0);
}
