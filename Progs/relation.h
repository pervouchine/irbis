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
#include "genutils.h"

class relation {
public:
    index_t	*index;
    index_t	*count;

    keys_t	*table;
    keys_t	length; // length of the index
    index_t     size;   //size of the table

    index_t	ptr1, ptr2;

    relation() {
	length = 0;
	size = 0;
    }

    void get_from_bin_file(char*, int);
    void get_from_tab_file(char*, int);

    void get_from_bin_file(FILE*, int);
    void get_from_tab_file(FILE*, int);

    void info();
    void report();

    void init();
    void make_pass1();
    void make_pass2();

/*  void setptr(keys_t);
    bool pass(keys_t);
*/
    bool check(keys_t, keys_t);
    void validate();
};

/*
void relation::setptr(keys_t a) {
    ptr1 = index[a];
    ptr2 = index[a+1];
}

bool relation::pass(keys_t b) {
    while(ptr1<ptr2 && table[ptr1]<b) ptr1++;
    return(ptr1<ptr2 && table[ptr1]==b);
}
*/

bool relation::check(keys_t a, keys_t b) {
    keys_t i;
    int flag = 0;
    if(a>=length) return(0);
    for(i=index[a];i<index[a+1];i++) {
	if(table[i]==b) flag=1;	
    }
    return(flag);
}


void relation::get_from_tab_file(char *filename, int right_end) {
    FILE *f;
    f = fopen(filename, "r");
    if(f==NULL) {
        if(verbose) fprintf(logfile,"Cannot load relation from '%s', exiting\n", filename);
        exit(1);
    }
    if(verbose) fprintf(logfile,"[Load from tabular %s",filename);
    get_from_tab_file(f, right_end);
    if(verbose) fprintf(logfile,"]\n");
    fclose(f);
}

void relation::get_from_bin_file(char *filename, int right_end) {
    FILE *f;
    f = fopen(filename, "rb");
    if(f==NULL) {
        if(verbose) fprintf(logfile,"Cannot load relation from '%s', exiting\n", filename);
        exit(1);
    }
    if(verbose) fprintf(logfile,"[Load from binary %s",filename);
    get_from_bin_file(f, right_end);
    if(verbose) fprintf(logfile,"]\n");
    fclose(f);
}


static void qs(keys_t arr[], int beg, int end) {
    int i = beg, j = end;
    keys_t tmp;
    keys_t pivot = arr[(beg + end) / 2];
    while (i <= j) {
        while (arr[i] < pivot) i++;
        while (pivot < arr[j]) j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    }
    if (beg < j) qs(arr, beg, j);
    if (i < end) qs(arr, i, end);
}

void relation::init() {
    keys_t i;
    index = (index_t*)malloc(sizeof(index_t)*(length + ARRAY_MARGIN));
    count = (index_t*)malloc(sizeof(index_t)*(length + ARRAY_MARGIN));
    if(index==NULL || count==NULL) {
        if(verbose) fprintf(logfile,"No space for relation index, exiting\n");
        exit(1);
    }
    for(i=0;i<length;i++) index[i] = count[i]= 0;
    if(verbose) fprintf(logfile,"(length=%i)",length);
}

void relation::make_pass1() {
    keys_t i;
    index_t a;
    for(size=0,i=0;i<length;i++) {
        a = index[i];
        index[i]=size;
        size+=a;
    }
    index[i]=size;
    if(verbose) fprintf(logfile,"(size=%li)",(long)size);
    table = (keys_t*) malloc(sizeof(keys_t)*(size + ARRAY_MARGIN));
    if(table == NULL) {
        if(verbose) fprintf(logfile,"No space for relation content");
        exit(1);
    }
}

void relation::make_pass2() {
    keys_t i;
    if(verbose) fprintf(logfile,"(sort)");
    for(i=0;i<length;i++) {
        qs(table, index[i], index[i+1] - 1);
    }
    free(count);
}

void relation::get_from_tab_file(FILE *f, int right_end) {
    keys_t left, right;
    keys_t i;
    index_t a;
    char buff[MAXBUFFLENGTH];

    if(verbose) fprintf(logfile,"(init)");
    length = 0;
    fseek(f,SEEK_SET,0);
    while(fgets(buff,MAXBUFFLENGTH,f)) {
        sscanf(buff,"%i %i",&left, &right);
	if(left>length) length = left;
    }
    length++;

    init();    

    if(verbose) fprintf(logfile,"(pass 1)");
    fseek(f,SEEK_SET,0);
    while(fgets(buff,MAXBUFFLENGTH,f)) {
	sscanf(buff,"%i %i",&left, &right);
        index[left]++;
    }
    make_pass1();

    if(verbose) fprintf(logfile,"(pass 2)");
    fseek(f,SEEK_SET,0);
    while(fgets(buff,MAXBUFFLENGTH,f)) {
        sscanf(buff,"%i %i",&left, &right);
        table[index[left]+count[left]]=right - right_end;
	count[left]++;
    }
    make_pass2();
    validate();
}

void relation::get_from_bin_file(FILE *f, int right_end) {
    keys_t left, right;
    keys_t i;
    index_t a;
    char buff[MAXBUFFLENGTH];

    if(verbose) fprintf(logfile,"(init)");
    length = 0;
    fseek(f,SEEK_SET,0);
    while(!feof(f)) {
	left = right = 0;
	fread(&left, 1, sizeof(int), f);
	fread(&right,1, sizeof(int), f);
	if(left<=0 || right<=0) break;
        if(left>length) length = left;
    }
    length++;

    init();

    if(verbose) fprintf(logfile,"(pass 1)");
    fseek(f,SEEK_SET,0);
    while(!feof(f)) {
        left = right = 0;
        fread(&left, 1, sizeof(int), f);
        fread(&right,1, sizeof(int), f);
        if(left<=0 || right<=0) break;
        index[left]++;
    }
    make_pass1();

    if(verbose) fprintf(logfile,"(pass 2)");
    fseek(f,SEEK_SET,0);
    while(!feof(f)) {
        left = right = 0;
        fread(&left, 1, sizeof(int), f);
        fread(&right,1, sizeof(int), f);
        if(left<=0 || right<=0) break;
        table[index[left]+count[left]]=right - right_end;
        count[left]++;
    }
    make_pass2();
    validate();
}



void relation::info() {
    if(verbose) fprintf(logfile, "(%i/%li)",length,(long)size);
}

void relation::validate() {
    keys_t i,j;
    for(i=0;i<length;i++) {
	for(j=index[i];j+1<index[i+1];j++) {
	    if(table[j]>table[j+1]) {
	    	fprintf(logfile,"Found inversion, exiting\n");
	    	exit(1);
	    }
	}
    }
    if(verbose) fprintf(logfile,"(x)");
}

void relation::report() {
    keys_t i,j;
    for(i=0;i<length;i++) {
	if(index[i]==index[i+1]) continue;
	printf("%i:\t",i);
	for(j=index[i];j<index[i+1];j++) {
	    printf("%i ",table[j]);
	}
	printf("\n");
    }
}
