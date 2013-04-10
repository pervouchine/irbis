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

FILE *alnfile;
FILE *cpsfile;
FILE *outfile;

int main(int argc, char* argv[]) {
    char alnfilename[MAXBUFFLENGTH];
    char cpsfilename[MAXBUFFLENGTH];
    char outfilename[MAXBUFFLENGTH]="";
 
    char c;

    char buff[MAXBUFFLENGTH];
    char aux[MAXBUFFLENGTH];

    int start1,end1,len1,start2,end2,len2;
    char strand1, strand2;

    int a,b,k,l,i,j,s,m,q;
    int d,dmin,lmin,score_max;
    int x,y;

    int qbest, kbest, jbest;
    int kprev;

    int *gene_idx;	//index
    int *gene_off;	//offset  
    int *gene_site;	//site number
    char *gene_styp;	//site type
    int *gene_pos;	//site position
    int *gene_chr;	//chromosome
    char *gene_str;	//strand

    int max_genes;
    int max_sites;

    int *site_idx;	//index
    int *site_off;	//offset
    int *site_chr;	//matching chromosome
    int *site_pos;	//matching pos
    int *site_str;	//matching strand
    int *site_score;	//--- optimal score
    int *site_lbest;    //--- where it came from
    int *site_qbest;    //--- where it came from

    int specific_site=0;

    double dthreshold = 0.50;
    int dlimit = 5000;
    int max_depth = 4;

    int pos, strand;

    if(argc==1) {
	fprintf(stderr,"Select best unique maping from the ALN file created by map_single\n");
        fprintf(stderr,"Last updated by Dmitri Pervouchine (dp@crg.eu) on Jan 28, 2013\n");
	fprintf(stderr,"Keys:\n -in <cps file>\n -aln <aln file>\n -out <output file>\n");
 	//fprintf(stderr," -l length difference limit [%i]\n -t percentage difference threshold [%1.2lf] (ONE OR THE OTHER THRESHOLD IS USED)\n -h max_depth [%i]\n",dlimit, dthreshold,max_depth);
	//fprintf(stderr," -v suppress verbose output [NO]\n");
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {
        if(strcmp(argv[i],"-in")==0) {
            sscanf(argv[++i], "%s", &cpsfilename[0]);
        }
        if(strcmp(argv[i],"-aln")==0) {
            sscanf(argv[++i], "%s", &alnfilename[0]);
        }
        if(strcmp(argv[i],"-out")==0) {
            sscanf(argv[++i], "%s", &outfilename[0]);
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
        if(strcmp(argv[i],"-s")==0) {
            sscanf(argv[++i], "%i", &specific_site);
        }
    }

    if(outfilename[0]==0) {
        fprintf(stderr,"[WARNING: output file not specified, redirect to stdout]\n");
        outfile = stdout;
    }
    else {
        outfile = fopen(outfilename,"w");
        if(outfile == NULL) {
            fprintf(stderr,"[ERROR: output file (%s) cannot be opened, exiting]\n", outfilename);
            exit(1);
        }
    }

/*******************************************************************************************************/
    cpsfile= fopen(cpsfilename,"r");
    if(cpsfile==NULL) {
	fprintf(stderr,"Can't access CPS file. Exiting\n");
	exit(1);
    }

    if(verbose) fprintf(stderr,"[Reading CPS, pass 0");
    max_sites = max_genes = 0;
    while(fgets(buff,MAXBUFFLENGTH,cpsfile)) {
      	if(strlen(buff)<2) break;
        sscanf(buff,"%*s %*i %*i %i %i", &i, &j);
        if(i>max_genes) max_genes = i;
	if(j>max_sites) max_sites = j;
    }

    max_genes++;
    max_sites++;

    gene_idx = (int*)  malloc(sizeof(int)*(max_genes+1));
    gene_off = (int*)  malloc(sizeof(int)*(max_genes+1));
    gene_chr = (int*)  malloc(sizeof(int)*(max_genes+1));
    gene_str = (char*) malloc(sizeof(char)*(max_genes+1));

    for(i=0;i<max_genes;i++) gene_idx[i]=gene_off[i]=0;

    if(verbose) fprintf(stderr,", max_genes = %i, max_sites = %i]\n", max_genes, max_sites);

    if(verbose) fprintf(stderr,"[Reading CPS, pass 1");
    fseek (cpsfile, 0, SEEK_SET);
    while(fgets(buff,MAXBUFFLENGTH,cpsfile)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s %*i %i %i %*i" , aux, &strand, &i);
        gene_idx[i]++;
	gene_chr[i] = assign_code(aux);
	gene_str[i] = strand;
    }

    for(s=i=0;i<max_genes;i++) {
        x = gene_idx[i];
        gene_idx[i] =s;
        s+=x;
    }
    gene_idx[i] = s;

    gene_site = (int*)  malloc(sizeof(int)*(s+1));
    gene_styp = (char*) malloc(sizeof(char)*(s+1));
    gene_pos  = (int*)  malloc(sizeof(int)*(s+1));

    if(verbose) fprintf(stderr,", records = %i]\n", s);

    if(verbose) fprintf(stderr,"[Reading CPS, pass 2");
    fseek (cpsfile, 0, SEEK_SET);
    while(fgets(buff,MAXBUFFLENGTH,cpsfile)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%*s %i %*i %i %i %c" , &x, &i, &j, &c);
	gene_site[gene_idx[i]+gene_off[i]]=j;
        gene_styp[gene_idx[i]+gene_off[i]]=c;
        gene_pos[gene_idx[i]+gene_off[i]]=x;
	gene_off[i]++;
    }
    fclose(cpsfile);
    if(verbose) fprintf(stderr,"]\n");
/**********************************************************************************************/

    alnfile = fopen(alnfilename,"r");
    if(alnfile == NULL) {
	fprintf(stderr, "Cant open alignment file, exiting\n");
	exit(1);
    }

    site_idx = (int*) malloc(sizeof(int)*(max_sites+1));
    site_off = (int*) malloc(sizeof(int)*(max_sites+1));

    for(i=0;i<max_sites;i++) {
        site_idx[i]=site_off[i]=0;
    }

    if(verbose) fprintf(stderr,"[Reading alignment file, pass 1");

    while(fgets(buff,MAXBUFFLENGTH,alnfile)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%*s %*i %*i %*s %*i %*i %*i %i %*c" , &i);
	site_idx[i]++;
    }

    for(s=i=0;i<max_sites;i++) {
        x = site_idx[i];
        site_idx[i] =s;
        s+=x;
    }
    site_idx[i] = s;

    site_chr   = (int*) malloc(sizeof(int)*(s+1));
    site_pos   = (int*) malloc(sizeof(int)*(s+1));
    site_str   = (int*) malloc(sizeof(int)*(s+1));
    site_score = (int*) malloc(sizeof(int)*(s+1));
    site_lbest = (int*) malloc(sizeof(int)*(s+1));
    site_qbest = (int*) malloc(sizeof(int)*(s+1));

    if(verbose) fprintf(stderr,", records = %i]\n",s);

    if(verbose) fprintf(stderr,"[Reading alignment file, pass 2");
    fseek (alnfile, 0, SEEK_SET);
    while(fgets(buff,MAXBUFFLENGTH,alnfile)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%*s %*i %*i %s %i %i %*i %i %*c" , &aux, &pos, &strand, &i);
	site_chr[site_idx[i]+site_off[i]] = assign_code(aux);
        site_pos[site_idx[i]+site_off[i]] = pos;
        site_str[site_idx[i]+site_off[i]] = strand;
	site_score[site_idx[i]+site_off[i]] = 0;
        site_lbest[site_idx[i]+site_off[i]] = -1;
        site_qbest[site_idx[i]+site_off[i]] = -1;
	site_off[i]++;
    }
    fclose(alnfile);
    if(verbose) fprintf(stderr,"]\n");

    for(i=0;i<max_genes;i++) {
	progressbar(i,max_genes-1, (char*)"Processing");
	score_max=0;
	kbest = -1;
	for(j=gene_idx[i];j<gene_idx[i+1];j++) {
	    x = gene_site[j];
	    for(k=site_idx[x];k<site_idx[x+1];k++) {
		site_score[k] = 0;
		site_lbest[k] = site_qbest[k] = -1;
	    }
	    for(q=1;q<=max_depth;q++) {
            	if(j-q>=gene_idx[i]) {
                    y = gene_site[j-q];
	   	    a = abs(gene_pos[j]-gene_pos[j-q]);
		    for(k=site_idx[x];k<site_idx[x+1];k++) {
			dmin = INFTY;
			lmin = -1;
			for(l=site_idx[y];l<site_idx[y+1];l++) {
			    if(site_chr[k] == site_chr[l] && site_str[k] == site_str[l]) {
			    	b = abs(site_pos[k]-site_pos[l]);
				d = abs(b-a);
				if(d<dmin) {
				    dmin = d;
				    lmin = l;
				}
				if(x==specific_site) fprintf(stderr,"[prev=%i curr=%i pos_p=%i pos_c=%i d=%i]\n",y,x,site_pos[l],site_pos[k],d);
			    }
			}
			m  = (lmin>=0 && (((double)dmin/a)<dthreshold || dmin<dlimit)) ? site_score[lmin] + a : 0;			
			if(m>site_score[k]) {
			    site_score[k] = m;
			    site_lbest[k] = lmin;
                            site_qbest[k] = q;
			}
			if(site_score[k]>score_max) {
			    score_max = site_score[k];
			    kbest = k; jbest = j;
			}
			if(x==specific_site) fprintf(stderr,"[curr=%i score=%i]\n",x,site_score[k]);
		    }
		}
	    }
	}
	j = jbest;
	k = kbest;
	if(k>=0 && site_score[k]>0) {
	    fprintf(outfile,"%s\t%i\t%i\t%s\t%i\t%i\t",get_chr_name(gene_chr[i]),gene_pos[j],gene_str[i],get_chr_name(site_chr[k]),site_pos[k],site_str[k]);
	    fprintf(outfile,"%i\t%i\t%c\t%i\t%i\n",i,gene_site[j],gene_styp[j],site_pos[site_lbest[k]],0);
	    while(site_score[k]>0 && site_lbest[k]>=0 && site_qbest[k]>=0) {
		kprev = k;
            	j = j - site_qbest[k];
                k = site_lbest[k];
 		fprintf(outfile,"%s\t%i\t%i\t%s\t%i\t%i\t",get_chr_name(gene_chr[i]),gene_pos[j],gene_str[i],get_chr_name(site_chr[k]),site_pos[k],site_str[k]);
		fprintf(outfile,"%i\t%i\t%c\t%i\t%i\n",i,gene_site[j],gene_styp[j],(site_lbest[k]>=0?site_pos[site_lbest[k]]:0),(kprev>=0?site_pos[kprev]:0));
	    }
 
	}
    }
    timestamp_report();
    exit(0);
}
