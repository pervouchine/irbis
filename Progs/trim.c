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

#include<dictionary.h>
#include<wordexp.h>
 
extern int lowercase_allowed;
extern int samestrandcontrol;

int main(int argc, char** argv) {
    dictionary<ATOMIC>* dict[MAXSPECIES];
    dictionary<ATOMIC>* filter;
    conservation_table* cons[MAXSPECIES];

    int i,j;
    FILE *configfile, *metacalc, *consfile;

    wordexp_t exp_result;

    char buff[MAXBUFFLENGTH];
    char tag[MAXBUFFLENGTH];
    char aux[MAXBUFFLENGTH];

    char path[MAXBUFFLENGTH]="";
    char suffix[MAXBUFFLENGTH]="";

    char configfilename[MAXBUFFLENGTH];
    char restrictionfilename[MAXBUFFLENGTH];
    char metacalcfilename[MAXBUFFLENGTH]="";
    char consfilename[MAXBUFFLENGTH]="";
    char maffilename[MAXBUFFLENGTH]="";

    int n_species=0;
    char* file_name[MAXSPECIES];
    double weight[MAXSPECIES];


    double threshold=0;
    double tot_weight=0;
    int halfsize=4;
    int gap=2;
    int min_gc=1;
    int repeat_length=2;
    int length_limit=0;

    subset *restriction=NULL;

    index_t p,q,r,s;
    word_t a;
    char *pc;
    char *p_file_name;

    int reverse_complement=0;


/***************************************************************************************************************************************/

    if(argc==1) {
	fprintf(stderr, "TRIM - this program takes sequences, creates corresponding hash tables, and intersects them by trimmming\n");
        fprintf(stderr, "Usage: %s -i config_file [-maf maf_file_name] -r restriction_file -o output_metafile\n",argv[0]);
	fprintf(stderr, "Other parameters are:\n -verbose [suppress verbose output]\n -rc [reverse complement sequences before adding]\n");
	exit(1);
    }

    timestamp_set();
    maffilename[0]=0;
    for(i=1;i<argc;i++) {
        pc = argv[i];
        if(*pc == '-') {
            if(strcmp(pc+1,"i") == 0) {
		configfile = fopen(argv[++i],"r");
		if(configfile==NULL) {
		    fprintf(logfile,"Input file (%s) cannot be opened, exiting", argv[i]);
		    exit(1);
		}
		n_species = 0;
		tot_weight= 0;
    		while(fgets(buff,MAXBUFFLENGTH,configfile)) {
        	    tag[0]=0;
        	    sscanf(buff,"%s %*s",tag);
        	    if(strcmp(tag,"path")==0) sscanf(buff,"%*s %s", &path[0]);
        	    if(strcmp(tag,"extention")==0)  sscanf(buff,"%*s %s", &suffix[0]);
        	    if(strcmp(tag,"halfsize")==0)   sscanf(buff,"%*s %i", &halfsize);
        	    if(strcmp(tag,"gap")==0)        sscanf(buff,"%*s %i", &gap);
        	    if(strcmp(tag,"repeatmask")==0) sscanf(buff,"%*s %i", &repeat_length);
        	    if(strcmp(tag,"gccontent")==0)  sscanf(buff,"%*s %i", &min_gc);
        	    if(strcmp(tag,"lowercase")==0)  sscanf(buff,"%*s %i", &lowercase_allowed);
        	    if(strcmp(tag,"threshold")==0)  sscanf(buff,"%*s %lf",&threshold);
		    if(strcmp(tag,"samewindow")==0) sscanf(buff,"%*s %i", &SAMEWINDOW);
		    if(strcmp(tag,"limit")==0)	    sscanf(buff,"%*s %i", &length_limit);
        	    if(strcmp(tag,"species")==0) {
			file_name[n_species] = (char*)malloc(sizeof(char)*MAXBUFFLENGTH);
			sscanf(buff,"%*s %s %lf", &file_name[n_species][0], &weight[n_species]);
			tot_weight+=weight[n_species];
            		n_species++;
		    }
		}
		fclose(configfile);
	    }
	    if(strcmp(pc+1,"maf") == 0) sscanf(argv[++i],"%s",&maffilename[0]);
	    if(strcmp(pc+1,"r") == 0)   restriction = new subset(argv[++i]);
	    if(strcmp(pc+1,"o") == 0) 	sscanf(argv[++i], "%s", &metacalcfilename[0]);
	    if(strcmp(pc+1,"v") == 0)	verbose = 0;
	    if(strcmp(pc+1,"rc") == 0)		reverse_complement = 1;

            if(strcmp(pc+1,"log")==0) {
            	logfile = fopen(argv[++i], "w");
            	if(logfile==NULL) logfile = stderr;
            }
	    
        }
    }

    if(halfsize<2 || gap>3) {
        if(verbose) fprintf(logfile,"Word size or gap are incorrect, exiting\n");
        exit(1);
    }

    if(metacalcfilename[0]==0) {
	fprintf(logfile,"Output file not specified, exiting\n");
        exit(1);
    }

    metacalc = fopen(metacalcfilename,"wb");
    if(metacalc==NULL) {
        if(verbose) fprintf(logfile,"Output metacalc file cannot be opened, exiting\n");
        exit(1);
    }

    strcpy(consfilename, metacalcfilename);
    strcat(consfilename, (char*)".cns");

    consfile = fopen(consfilename,"wb");
    if(consfile==NULL) {
        if(verbose) fprintf(logfile,"Output conservation file cannot be opened, exiting\n");
        exit(1);
    }




/***************************************************************************************************************************************/

    if(threshold==0) threshold = 0.5;

    if(verbose) fprintf(logfile, "[Numeric parameters are: pattern=%i-%i-%i, minGC=%i, repeat=%i, threshold=%2.1lf%%]\n", halfsize, gap, halfsize, min_gc, repeat_length, threshold*100);
    if(verbose) fprintf(logfile, "[Boolean parameters are: lowercase=%s, reverse_complement=%s]\n",yesno(lowercase_allowed), yesno(reverse_complement)); 
    if(verbose) fprintf(logfile, "[Extra parameters are: SAMEWINDOW=%i, sequence length limit=%i]\n",SAMEWINDOW,length_limit); 

    wordexp(path, &exp_result, 0);
    sprintf(&path[0],"%s",exp_result.we_wordv[0]);


/***************************************************************************************************************************************/

    filter = new dictionary<ATOMIC>(halfsize,gap);
    filter->mask_low_complexity(repeat_length);
    if(min_gc>0) filter->mask_low_GCcontent(min_gc);
    for(i=0;i<n_species;i++) {
        dict[i] = new dictionary<ATOMIC>(halfsize,gap);
        dict[i]->filter = filter;
    }

    if(maffilename[0]==0) {
    	for(i=0;i<n_species;i++) {
	    p_file_name = getfullname(path,file_name[i],suffix);
	    if(verbose) fprintf(logfile,"[Reading %s: ",p_file_name);
	    dict[i]->read_from_suf(p_file_name, restriction, reverse_complement,length_limit);
	    dict[i]->check();
	    cons[i] = new conservation_table(dict[i]->max_key);
	    cons[i]->before = dict[i]->fill_cons();
            if(verbose) fprintf(logfile,", weight=%2.1lf]\n",weight[i]);
	}
    }
    else {
    	if(verbose) fprintf(logfile,"[Reading %s: ", maffilename);
    	dictionary<ATOMIC>::read_from_muf(maffilename, dict, restriction, file_name, n_species, reverse_complement,length_limit);
    	if(verbose) fprintf(logfile,"]\n");
    	for(i=0;i<n_species;i++) {
            if(verbose) fprintf(logfile,"[%s",file_name[i]);
            dict[i]->info();
            dict[i]->check();
            cons[i] = new conservation_table(dict[i]->max_key);
            cons[i]->before = dict[i]->fill_cons();
            if(verbose) fprintf(logfile,", weight=%2.1lf]\n",weight[i]);
    	}
    }

    dictionary<ATOMIC>::intersect_many(dict, n_species, weight, threshold*tot_weight, (char*)"Intersect ");

    if(verbose) fprintf(logfile,"[Compute conservation rates:");
    for(i=0;i<n_species;i++) {
	if(verbose) fprintf(logfile," %s",file_name[i]);
	cons[i]->after = dict[i]->fill_cons();
    }
    if(verbose) fprintf(logfile,"]\n");

    if(verbose) fprintf(logfile,"[Metacalc: saving data to %s]\n",metacalcfilename);
    fprintf(metacalc,"%i %i %i %lf\n", halfsize, gap, n_species, threshold);
    for(i=0;i<n_species;i++) {
    	if(verbose) fprintf(logfile,"[Saving %s",file_name[i]);
	dict[i]->check();
	dict[i]->info();
	fprintf(metacalc,"%s %lf\n", file_name[i], weight[i]);
	dict[i]->save(metacalc);
	if(i==0) cons[i]->save(consfile);
	if(verbose) fprintf(logfile,"]\n");
    }

    fclose(metacalc);
    fclose(consfile);
    timestamp_report();
    return(0);
}
