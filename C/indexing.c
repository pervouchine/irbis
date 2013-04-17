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

#include "genutils.h"
#include "progressbar.h"


int main(int argc, char* argv[]) {
    char buff[MAXLONGBUFFLENGTH];
    char inp_file_name[MAXLONGBUFFLENGTH]="";
    char out_file_name[MAXLONGBUFFLENGTH]="";
    long pos;

    FILE *inp_file, *out_file;
    char *pc;
    int i, n, key;
    char *ptr;

    if(argc==1) {
	fprintf(stderr,"This script is to create index for sus/suw files\n");
        fprintf(stderr,"Last update by Dmitri Pervouchine (dp@crg.eu) on Mar 22, 2013\n");
        fprintf(stderr,"Usage: %s -in <file_name> -out <index_file_name>\n",argv[0]);
        exit(1);
    }

    for(i=1;i<argc;i++) {
	pc = argv[i];
        if(*pc != '-') continue;
        if(strcmp(pc+1,"in")==0)  sscanf(argv[++i], "%s",  &inp_file_name[0]);
	if(strcmp(pc+1,"out")==0) sscanf(argv[++i], "%s",  &out_file_name[0]);
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

    inp_file = fopen(inp_file_name, "r");
    if(inp_file==NULL) {
        fprintf(stderr,"[ERROR: cannot access %s, exiting]\n", inp_file_name);
        exit(1);
    }

    fprintf(stderr, "[<%s", inp_file_name);
    while(fgets(buff,MAXLONGBUFFLENGTH-ARRAY_MARGIN,inp_file)) {
        sscanf(buff,"%i",&key);
	fprintf(out_file,"%i\t%li\n", key, pos);
	pos = ftell(inp_file);
    }
    fclose(inp_file);
    fclose(out_file);
    fprintf(stderr, "]\n");
    exit(0);
}

