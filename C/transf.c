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

void compress(FILE *fastafile, FILE *maskfastafile, FILE *idxfile, FILE *dbxfile);
void uncompress(FILE *fastafile, FILE *maskfastafile, FILE *idxfile, FILE *dbxfile);

int main(int argc, char* argv[]) {
    int UNCOMPRESS=0;
    int maskdatafound=0;

    char buff[MAXBUFFLENGTH]="";
    char filename[MAXBUFFLENGTH]="";

    char dirname[MAXBUFFLENGTH]="";
    char maskdirname[MAXBUFFLENGTH]="";
    char maskextention[MAXBUFFLENGTH]="";

    char idxfilename[MAXBUFFLENGTH]="";
    char dbxfilename[MAXBUFFLENGTH]="";

    char fastafilename[MAXBUFFLENGTH];
    char maskfastafilename[MAXBUFFLENGTH];

    int remove_source=0;
    struct direct **files;
    int i,n,m;
    char *ptr;

    FILE* dbxfile=NULL;
    FILE* idxfile=NULL;
    FILE* fastafile=NULL;
    FILE* maskfastafile=NULL;


    if(argc==1) {
        fprintf(stderr,"Custom compress/decompress utility for storing masked and unmasked fasta files\n");
        fprintf(stderr,"Last update by Dmitri Pervouchine (dp@crg.eu) on Mar 22, 2013\n");
        fprintf(stderr," -dir name of the directory containing FASTA files\n");
	fprintf(stderr," -maskdir name of the directory containing masked FASTA files (optional)\n");
	fprintf(stderr," -maskext masked FASTA file extension\n");
	fprintf(stderr," -dbx, -idx sequence database names\n");
        fprintf(stderr," -remove remove source FASTA files [default=NO]\n");
	fprintf(stderr," -uncompress [default=NO]\n");
	fprintf(stderr," -quiet suppress verbose output [default=NO]\n");
        exit(1);
    }

    timestamp_set();

    for(i=1;i<argc;i++) {
        if(strcmp(argv[i],"-dir")==0) {
            sscanf(argv[++i], "%s", &dirname[0]);
        }
        if(strcmp(argv[i],"-dbx")==0) {
            sscanf(argv[++i], "%s", &dbxfilename[0]);
        }
        if(strcmp(argv[i],"-idx")==0) {
            sscanf(argv[++i], "%s", &idxfilename[0]);
        }
        if(strcmp(argv[i],"-remove")==0) {
	    remove_source = 1;
        }
        if(strcmp(argv[i],"-quiet")==0) {
            verbose = 0;
        }
	if(strcmp(argv[i],"-maskdir")==0) {
            sscanf(argv[++i], "%s", &maskdirname[0]);
        }
        if(strcmp(argv[i],"-maskext")==0) {
            sscanf(argv[++i], "%s", &maskextention[0]);
        }
        if(strcmp(argv[i],"-uncompress")==0) {
	    UNCOMPRESS = 1;
        }   

    }

    if(maskdirname[0]==0) strcpy(maskdirname, dirname);

    if(dbxfilename[0]==0 || idxfilename[0]==0) {
        if(verbose) fprintf(stderr,"Index or database file name not provided, exiting\n");
        exit(1);
    }

    endslashdir(dirname);
    endslashdir(maskdirname);

    if(UNCOMPRESS) {
        idxfile = fopen(idxfilename,"r");
        dbxfile = fopen(dbxfilename,"rb");
        if(idxfile == NULL || dbxfile == NULL) {
            if(verbose) fprintf(stderr,"DB cannot be open for reading, exiting\n");
            exit(1);
        }
        if(verbose) fprintf(stderr,"Input set to:\t%s,%s\n",idxfilename,dbxfilename);
	while(fgets(buff,MAXBUFFLENGTH,idxfile)) {
	    if(strlen(buff)<2) break;
	    maskdatafound = 0;
            sscanf(buff,"%s %i",&filename[0],&maskdatafound);
            strcpy(fastafilename,dirname);
            strcat(fastafilename,filename);
	    if(verbose) fprintf(stderr,"[%s",fastafilename);
	    fastafile = fopen(fastafilename,"w");
	    if(fastafile==NULL) {
		if(verbose) fprintf(stderr,"] - error, exiting\n");
		exit(1);
	    }
	    if(maskdatafound) {
            	strcpy(maskfastafilename,dirname);
            	strcat(maskfastafilename,filename);
		strcat(maskfastafilename,".masked");
                if(verbose) fprintf(stderr,", %s",maskfastafilename);
		maskfastafile = fopen(maskfastafilename,"w");
		if(maskfastafile == NULL) {
		    if(verbose) fprintf(stderr," can't open, exiting\n");
		    exit(1);
		}
	    }
	    uncompress(fastafile,maskfastafile,idxfile,dbxfile);
	    if(verbose) fprintf(stderr,"]\n");
	}
        if(remove_source) {
            if(verbose) fprintf(stderr,"Removing source databases\n");
	    remove(idxfilename);
            remove(dbxfilename);
	}
    }
    else {
        n=scandir(dirname,&files,0,alphasort);
	for(m=i=0;i<n;i++) {
	    ptr = rindex(files[i]->d_name,'.');
            if(ptr != NULL && strcmp(ptr, ".fa")==0) m++;
	}
	if(m==0) {
	    if(verbose) fprintf(stderr,"The directory doesn't cointain any fasta files, exiting\n");
	    exit(0);
	}

        idxfile = fopen(idxfilename,"w");
    	dbxfile = fopen(dbxfilename,"wb");
    	if(idxfile == NULL || dbxfile == NULL) {
	    if(verbose) fprintf(stderr,"DB cannot be open for writing, exiting\n");
            exit(1);
        }
    	if(verbose) fprintf(stderr,"Output set to:\t%s,%s\n",idxfilename,dbxfilename);

    	for(i=0;i<n;i++) {
            ptr = rindex(files[i]->d_name,'.');
            if(ptr != NULL && strcmp(ptr, ".fa")==0) {
            	strcpy(fastafilename,dirname);
            	strcat(fastafilename,files[i]->d_name);
	    	fastafile = fopen(fastafilename,"r");
	    	if(fastafile==NULL) {
		    if(verbose) fprintf(stderr," can't open, skipping\n");
		    continue;
	    	}
	    	if(verbose) fprintf(stderr,"[%s",fastafilename);
		strcpy(maskfastafilename,maskdirname);
	    	strcat(maskfastafilename,files[i]->d_name);
		strcat(maskfastafilename, (strcmp(dirname,maskdirname)==0) ? ".masked" : maskextention);
		maskfastafile = fopen(maskfastafilename,"r");
		if(maskfastafile!=NULL) {
		    if(verbose) fprintf(stderr,", %s",maskfastafilename);
            	}
		fprintf(idxfile,"%s %i\n",files[i]->d_name,(maskfastafile==NULL)?0:1);
	        compress(fastafile,maskfastafile,idxfile,dbxfile);
                fprintf(idxfile,"\n");
	    	fclose(fastafile);
	    	if(maskfastafile!=NULL) fclose(maskfastafile);
		if(remove_source) {
		    if(verbose) fprintf(stderr," remove sources");
		    remove(fastafilename);
		    if(maskfastafile!=NULL) remove(maskfastafilename);
		}
	    	if(verbose) fprintf(stderr,"]\n");
	    }
    	}
    	fclose(idxfile);
        fclose(dbxfile);
    }
    timestamp_report();
    exit(0);
}

void compress(FILE *fastafile, FILE *maskfastafile, FILE *idxfile, FILE *dbxfile) {
    char chr_name[MAXBUFFLENGTH];
    char buff[MAXBUFFLENGTH];
    char maskbuff[MAXBUFFLENGTH];
    char c,cm;
    long ipos;
    int m;
    long filepos;
    int n,b,t;
    while(fgets(buff,MAXBUFFLENGTH,fastafile)) {
	if(maskfastafile!=NULL) fgets(maskbuff,MAXBUFFLENGTH,maskfastafile);
	if(strlen(buff)<1) break;
        if(buff[0]=='>') {
	    buff[0]=0;
            sscanf(buff+1,"%s",chr_name);
	    ipos=0;
            filepos = ftell(fastafile);
            while(!feof(fastafile)) {
                c = fgetc(fastafile);
                if(c=='>') break;
                if(c>' ') ipos++;
            }
            fseek(fastafile,filepos,SEEK_SET);
	    fprintf(idxfile,"%s %li\n",chr_name,ipos);
	    fcode_start(dbxfile);
	    while(ipos>0) {
                c = fgetc(fastafile);
		if(maskfastafile!=NULL) cm = fgetc(maskfastafile);
                if(c>' ') {
		    fcode(c,cm,dbxfile);
		    ipos--;
		}
            }
	    fcode_stop(dbxfile);
	}
    }
}
	                
void uncompress(FILE *fastafile, FILE *maskfastafile, FILE *idxfile, FILE *dbxfile) {
    char c;
    char chr_name[MAXBUFFLENGTH];
    char buff[MAXBUFFLENGTH];
    int m;
    long ipos;
    long filepos;
    int n,b,t;

    while(fgets(buff,MAXBUFFLENGTH,idxfile)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s %li",&chr_name[0],&ipos);
	fprintf(fastafile,">%s\n",chr_name);
	if(maskfastafile!=NULL) fprintf(maskfastafile,">%s\n",chr_name);
	while(ipos!=0) {
	    fread(&t,sizeof(int),1,dbxfile);
	    fdecode(t,&ipos,fastafile,maskfastafile);
	}
	fprintf(fastafile,"\n");
        if(maskfastafile!=NULL) fprintf(maskfastafile,"\n");
    }
}
