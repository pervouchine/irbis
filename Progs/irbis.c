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

#include <dictionary.h>
#include <relation.h>

extern int lowercase_allowed;
extern int SAMEWINDOW;


int main(int argc, char** argv) {
    char buff[MAXBUFFLENGTH];
    char aux[MAXBUFFLENGTH];

    dictionary<ATOMIC>* dict_f[MAXSPECIES];
    dictionary<ATOMIC>* dict_r[MAXSPECIES];
    dictionary<PAIR<ATOMIC> >* dict_bp[MAXSPECIES];
    dictionary<ATOMIC>* filter;

    conservation_table* cons_f[MAXSPECIES];
    conservation_table* cons_r[MAXSPECIES];

    char right_file_name[MAXBUFFLENGTH]="";
    char left_file_name[MAXBUFFLENGTH]="";
    char right_cons_file_name[MAXBUFFLENGTH]="";
    char left_cons_file_name[MAXBUFFLENGTH]="";
    char right_sbs_file_name[MAXBUFFLENGTH]="";
    char left_sbs_file_name[MAXBUFFLENGTH]="";

    char rel_binary_file_name[MAXBUFFLENGTH]="";
    char rel_tabdel_file_name[MAXBUFFLENGTH]="";

    relation* rel=NULL;

    FILE *right_file, *left_file, *rel_file, *outfile;
    FILE *right_cons_file, *left_cons_file; 

    char* specie_name[MAXSPECIES];
    double weight[MAXSPECIES];

    char filename[MAXBUFFLENGTH];
    char outfilename[MAXBUFFLENGTH];

    double threshold, threshold1, threshold2;
    int halfsize,gap,n_species;
    int halfsize1,gap1,n_species1;
    double tot_weight=0;
    int right_end=0;

    int redundant_pairs = 2;
    const char RED_P[3][24]={"strictly not allowed", "not allowed", "allowed"};

    double min_cons = 0;
    int min_length   = 0;
    int joint_fold = 1;
    int number_of_gt_basepairs=0;

    index_t p,q,r,s;
    word_t a,b;
    int i,j;
    int repeat_length=1;
    int min_gc=0;
    char *pc;
    word_t aux_array[GT_ARRAY_CAPACITY];
    int aux_array_length;
    int n,m;

    subset *left_restriction=NULL;
    subset *right_restriction=NULL;

    if(argc==1) {
	fprintf(stderr, "BInary Search for COnserved RNa Structures\n");
	fprintf(stderr, "Usage: %s -l [left_metacalc] -r [right_metacalc] -o [output_tabular file]\n",argv[0]);
	fprintf(stderr, " -v suppress verbose output\n -t threshold for intersection\n -u redundant pairs (0=strict no, 1=no, 2=allowed)\n -L min length\n");
	fprintf(stderr, " -B binary relation file\n -b tab-delimited relation file\n -j suppress joint fold (fold species separately)\n");
	exit(0);
    }

    timestamp_set();
    rel_binary_file_name[0] = rel_tabdel_file_name[0] = 0;
    for(i=1;i<argc;i++) {
        pc = argv[i];
        if(*pc == '-') {
            if(strcmp(pc+1,"v") == 0) verbose=0;
	    if(strcmp(pc+1,"l") == 0)  sscanf(argv[++i], "%s",  &left_file_name[0]);
	    if(strcmp(pc+1,"r") == 0)  sscanf(argv[++i], "%s",  &right_file_name[0]);
            if(strcmp(pc+1,"lr")== 0)  sscanf(argv[++i], "%s",  &left_sbs_file_name[0]);
            if(strcmp(pc+1,"rr")== 0)  sscanf(argv[++i], "%s",  &right_sbs_file_name[0]);
            if(strcmp(pc+1,"B") == 0)  sscanf(argv[++i], "%s",  &rel_binary_file_name[0]);
            if(strcmp(pc+1,"b") == 0)  sscanf(argv[++i], "%s",  &rel_tabdel_file_name[0]);
	    if(strcmp(pc+1,"o") == 0)  sscanf(argv[++i], "%s",  &outfilename[0]);
	    if(strcmp(pc+1,"t") == 0)  sscanf(argv[++i], "%lf", &threshold);
	    if(strcmp(pc+1,"L") == 0)  sscanf(argv[++i], "%i",  &min_length);
	    if(strcmp(pc+1,"C") == 0)  sscanf(argv[++i], "%lf", &min_cons);
	    if(strcmp(pc+1,"u") == 0)  sscanf(argv[++i], "%i",  &redundant_pairs);
	    if(strcmp(pc+1,"g") == 0)  sscanf(argv[++i], "%i",  &number_of_gt_basepairs);
	    if(strcmp(pc+1,"E") == 0)  right_end=1;
	    if(strcmp(pc+1,"j") == 0)  joint_fold = 0;
        }
    }

    outfile = fopen(outfilename,"w");
    if(outfile==NULL) {
        if(verbose) fprintf(logfile,"Cannot open outfile, exiting\n");
        exit(1);
    }

    if(left_sbs_file_name[0]!=0)  left_restriction  = new subset(left_sbs_file_name);
    if(right_sbs_file_name[0]!=0) right_restriction = new subset(right_sbs_file_name);

    if(redundant_pairs>2) redundant_pairs=2;
    if(verbose) fprintf(logfile,"[Params: max_gt=%i, threshold=%2.1lf%%, min_length=%i, min_cons_score=%2.1lf, redundant_pairs=%s, joint_fold=%s]\n",
            number_of_gt_basepairs, threshold*100, min_length, min_cons, RED_P[redundant_pairs], yesno(joint_fold));

    if(rel_binary_file_name[0]!=0) {
	rel = new relation();
	rel->get_from_bin_file(rel_binary_file_name, right_end);
    }

    if(rel_tabdel_file_name[0]!=0) {
        rel = new relation();
        rel->get_from_tab_file(rel_tabdel_file_name, right_end);
    }

    if(verbose) fprintf(logfile,"Reading left %s", left_file_name);
    left_file = fopen(left_file_name,"r");
    if(left_file==NULL) {
	fprintf(stderr,"cannot be found, exiting\n");
	exit(1);
    }
    strcpy(left_cons_file_name, left_file_name);
    strcat(left_cons_file_name, (char*)".cns");
    left_cons_file = fopen(left_cons_file_name,"r");
    if(left_cons_file == NULL) {
        fprintf(stderr, "cannot be found, exiting\n");
        exit(1);
    }

    fgets(buff, MAXBUFFLENGTH, left_file);
    sscanf(buff,"%i %i %i %lf",&halfsize, &gap, &n_species, &threshold1);
    if(verbose) fprintf(logfile,"\n");
    for(i=0;i<n_species;i++) {
        dict_f[i] = new dictionary<ATOMIC>(halfsize,gap);
	fgets(buff, MAXBUFFLENGTH, left_file);
	specie_name[i] = (char*)malloc(sizeof(char)*MAXBUFFLENGTH);
	sscanf(buff,"%s %lf",&specie_name[i][0],&weight[i]);
 	tot_weight+=weight[i];
        if(verbose) fprintf(logfile,"[Loading %s ",specie_name[i]);
        dict_f[i]->load(left_file);
        cons_f[i] = new conservation_table();
        if(i==0) cons_f[i]->load(left_cons_file);
        dict_f[i]->info();
	if(left_restriction!=NULL) {
	    dict_f[i]->subset_on(left_restriction);
	}
        dict_f[i]->check();
        if(verbose) fprintf(logfile,"]\n");
    }

    if(verbose) fprintf(logfile,"Reading right %s", right_file_name);
    right_file = fopen(right_file_name,"r");
    if(right_file == NULL) {
	fprintf(stderr, "cannot be found, exiting\n");
	exit(1);
    }
    strcpy(right_cons_file_name, right_file_name);
    strcat(right_cons_file_name, (char*)".cns");
    right_cons_file = fopen(right_cons_file_name, "r");
    if(right_cons_file == NULL) {
	fprintf(stderr," conservation file cannot be found, exiting\n");
        exit(1);
    }

    if(verbose) fprintf(logfile,"\n");
    fgets(buff, MAXBUFFLENGTH, right_file);
    sscanf(buff,"%i %i %i %lf",&halfsize1, &gap1, &n_species1, &threshold2);
    if(halfsize!=halfsize1 || gap!=gap1) {
	fprintf(logfile,"Dictionaries not compartible (%i %i) and (%i %i), exiting\n", halfsize, gap, halfsize1, gap1);
	exit(1);
    }
    if(n_species!=n_species1) {
	fprintf(logfile,"Number of species are different");
	exit(1);
    }
    for(i=0;i<n_species;i++) {
        dict_r[i] = new dictionary<ATOMIC>(halfsize,gap);
        fgets(buff, MAXBUFFLENGTH, right_file);
        sscanf(buff,"%s",&aux[0]);
	if(strcmp(aux,specie_name[i])!=0) {
	    fprintf(logfile,"Species are not the same (%s vs %s), exiting\n", aux, specie_name[i]);
	    exit(1);
	}
        if(verbose) fprintf(logfile,"[Loading %s ",specie_name[i]);
        dict_r[i]->load(right_file);
        cons_r[i] = new conservation_table();
        if(i==0) cons_r[i]->load(right_cons_file);
        dict_r[i]->info();
	if(right_restriction!=NULL) {
	    dict_r[i]->subset_on(right_restriction);
	}
        dict_r[i]->check();
        if(verbose) fprintf(logfile,"]\n");
    }

    if(threshold<threshold1 || threshold<threshold2) {
        threshold=MAX2(threshold1, threshold2);
        if(verbose) fprintf(logfile,"[WARNING: product threshold set to %2.2lf%%]\n",threshold*100);
    }

    if((threshold1<threshold || threshold2<threshold) && verbose) fprintf(logfile,"Threshold rose to %2.2lf%%, need to trim again\n",threshold*100);
    if(threshold1<threshold) dictionary<ATOMIC>::intersect_many(dict_f, n_species, weight, threshold*tot_weight, (char*)"Forward\t");
    if(threshold2<threshold) dictionary<ATOMIC>::intersect_many(dict_r, n_species, weight, threshold*tot_weight, (char*)"Reverse\t");

    for(i=0;i<n_species;i++) {
	if(verbose) fprintf(logfile,"[%s forward ",specie_name[i]);
	dict_f[i]->info();
	if(verbose) fprintf(logfile,", reverse: ");
        dict_r[i]->info();
        if(verbose) fprintf(logfile,", product: ");

	/** times procedure **/
    	dict_bp[i] = new dictionary<PAIR<ATOMIC> >(dict_f[i]->half_size,dict_f[i]->gap_size);
    	for(s=0,a=0;a<=dict_f[i]->last_word;a++) {
	    dict_bp[i]->revcomp_gt(a, number_of_gt_basepairs, aux_array, &aux_array_length);
	    for(;aux_array_length>0;aux_array_length--) { 
	    	b = aux_array[aux_array_length-1];
	    	if(rel) {
	    	    for(p=dict_f[i]->index[a];p<dict_f[i]->index[a+1];p++) {
		    	for(q=dict_r[i]->index[b];q<dict_r[i]->index[b+1];q++) {
		    	    if(rel->check(dict_f[i]->table[p].getid(), dict_r[i]->table[q].getid())) s++;
		    	}
		    }
	    	}
	    	else {
	    	    s+=(dict_f[i]->index[a+1]-dict_f[i]->index[a])*(dict_r[i]->index[b+1]-dict_r[i]->index[b]);
	    	}
	    }
	}

    	dict_bp[i]->table = (PAIR<ATOMIC>*)malloc((size_t)sizeof(PAIR<ATOMIC>)*(s + ARRAY_MARGIN));
    	r=0;
    	for(a=0;a<=dict_f[i]->last_word;a++) {
	    dict_bp[i]->index[a]=r;
	    dict_bp[i]->revcomp_gt(a, number_of_gt_basepairs, aux_array, &aux_array_length);
            for(;aux_array_length>0;aux_array_length--) {
	    	b = aux_array[aux_array_length-1];
            	for(p=dict_f[i]->index[a];p<dict_f[i]->index[a+1];p++) {
            	    for(q=dict_r[i]->index[b];q<dict_r[i]->index[b+1];q++) {
		    	if(rel) if(!rel->check(dict_f[i]->table[p].getid(), dict_r[i]->table[q].getid())) continue;
                    	if(redundant_pairs==2 || 
			   redundant_pairs==1 && (dict_f[i]->table[p].getid()<dict_r[i]->table[q].getid() || 
						  dict_f[i]->table[p].getid()==dict_r[i]->table[q].getid() && dict_f[i]->table[p].getpos() + dict_f[i]->full_size<dict_r[i]->table[q].getpos()) || 
			   redundant_pairs==0 && dict_f[i]->table[p].getid()<dict_r[i]->table[q].getid()
			  ) 
			    dict_bp[i]->table[r++].set(dict_f[i]->table[p], dict_r[i]->table[q]);
            	    }
                }
	    }
	    PAIR<ATOMIC>::quicksort_less(dict_bp[i]->table, dict_bp[i]->index[a], r-1); 
    	}
    	dict_bp[i]->index[a]=r;

    	if(verbose) fprintf(logfile,"%li/%li/%li ",(long)dict_bp[i]->last_word,(long)r,(long)s);
        free(dict_bp[i]->count);
        dict_bp[i]->phase = PHASE_COMPLETED;
	/** end of times procedure **/

	dict_bp[i]->info();
	dict_bp[i]->check();
	dict_f[i]->drop();
        dict_r[i]->drop();
	if(verbose) fprintf(logfile,"]\n");
    }


    dictionary<PAIR<ATOMIC> >::intersect_many(dict_bp, n_species, weight, threshold*tot_weight, (char*)"Intersection of products \t");

    if(joint_fold) {
    	dictionary<PAIR<ATOMIC> >::foldall(outfile, dict_bp, specie_name, n_species, weight, threshold*tot_weight, cons_f[0], cons_r[0], min_cons, min_length);
    }
    else {
    	if(verbose) fprintf(logfile,"[Folding and saving to %s...]\n",outfilename);
    	for(i=0;i<n_species;i++) {
	    if(verbose) fprintf(logfile,"[%s ",specie_name[i]);
	    dict_bp[i]->info();
	    n = dict_bp[i]->fold_and_save_standalone(outfile, i, cons_f[0], cons_r[0], min_length);
            if(verbose) fprintf(logfile,", %i structures]\n",n);
    	}
    }
    fclose(outfile);
    timestamp_report();
    return(0);
}
