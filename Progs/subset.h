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

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include<genutils.h>

class subset {
public:
    keys_t	*index;
    keys_t	length; // length of the index

    subset() {
	length = 0;
    }
    subset(char* filename) {
	from_cps(filename);
    }

    void from_cps(FILE*);
    void from_cps(char*);
    void info();
};


void subset::from_cps(char *filename) {
    FILE *f;
    f = fopen(filename, "rb");
    if(f==NULL) {
        if(verbose) fprintf(logfile,"Cannot load subset from '%s', exiting\n", filename);
        exit(1);
    }
    if(verbose) fprintf(logfile,"[Load from %s",filename);
    from_cps(f);
    if(verbose) fprintf(logfile,"]\n");
    fclose(f);
}

void subset::from_cps(FILE *f) {
    keys_t key,i;
    char buff[MAXBUFFLENGTH];
    length=0;

    fseek(f,SEEK_SET,0);
    while(fgets(buff,MAXBUFFLENGTH,f)) {
        sscanf(buff,"%*s %*i %*i %*i %i",&key);
        if(key>length) length = key;
    }
    length++;

    index = (keys_t*)malloc(sizeof(keys_t)*(length + ARRAY_MARGIN));

    if(index==NULL) {
        if(verbose) fprintf(logfile,"No space for subset index, exiting\n");
        exit(1);
    }
    for(i=0;i<length;i++) index[i] = 0;

    fseek(f,SEEK_SET,0);
    while(fgets(buff,MAXBUFFLENGTH,f)) {
	sscanf(buff,"%*s %*i %*i %*i %i",&key);
        index[key]++;
    }
}

void subset::info() {
    if(verbose) fprintf(logfile, "(%i)",length);
}
