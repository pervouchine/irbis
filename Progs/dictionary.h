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
#include <orderedset.h>
#include <subset.h>
#include <linearfold.h>
#include <conservation.h>

#define PHASE_CONTAINERONLY 0
#define PHASE_READY1 1
#define PHASE_READY2 2
#define PHASE_COMPLETED 3
#define DICT_VERSION 301

#define MAXFILEBUFFLENGTH 600000

using namespace std;

static char PHASE[][32] = {"not ready", "ready for pass 1", "ready for pass 2", "completed"};

template <class ordered_set>
class dictionary {
public:
    int word_size;
    int gap_size;
    int half_size;
    int full_size;
    int phase;
    word_t upper_mask;
    word_t lower_mask;
    word_t last_word;
    index_t	*index;
    index_t	*count;
    ordered_set *table;
    keys_t	max_key;
    dictionary  *filter;

public:
    dictionary() {
        phase = PHASE_CONTAINERONLY;
    }

    dictionary(int n) {
	init(n, 0);
    }

    dictionary(int n, int g) {
	init(n, g);
    }

    void init(int, int);
    void pass1_add(char *seq, int length_limit);
    void pass1_make();
    void pass2_add(char *seq, keys_t key, int length_limit);
    void pass2_make();
    void drop();  
    void read_from_suf(char*, int reverse_complement,int length_limit);
    void read_from_suf(FILE*, int reverse_complement,int length_limit);
    void read_from_suf(char*, subset*, int reverse_complement, int length_limit);
    void read_from_suf(FILE*, subset*, int reverse_complement, int length_limit);
    void subset_on(subset*);

    static void read_from_muf(char*, dictionary**, subset*, char**, int, int, int length_limit);
    static void read_from_muf(FILE*, dictionary**, subset*, char**, int, int, int length_limit);

    static void intersect_many(dictionary**, int, double*, double);
    static void intersect_many(dictionary**, int, double*, double,char*);
    void save(FILE *f);
    void load(FILE *f);
    void save(char*);
    void load(char*);
    void fprint(FILE *f);
    void print();
    word_t reverse_word(word_t w);
    void print_word(word_t w);
    void fprint_word(FILE*, word_t);
    void describe();
    void info();
    void check();

    void finalize();

    static void foldall(FILE *f, dictionary **dict, char ** name, int n, double *weight, double threshold, conservation_table *cons_f, conservation_table *cons_r, double min_cons, int report_length);

    int  fold_and_save_standalone(FILE*, int species, conservation_table*, conservation_table*, int report_length);
    void fold_simple(char *buff, word_t* word, int species, index_t i, index_t j, conservation_table* cons_f, conservation_table *cons_r, int *length);

    int* fill_cons();

    void mask_low_complexity(int);
    void mask_low_GCcontent(int);
    bool passfilter(word_t);

    void var_ac_sub(word_t w, word_t a, word_t c, word_t g, word_t t, int pos, int max_mutations, word_t *result, int *n_results);
    void revcomp_gt(word_t w, int n_gt_pairs, word_t *result, int *n_results);
};


/***********************************************************************************************************/

template<class ordered_set>
void dictionary<ordered_set>::init(int n, int g) {
    word_t i;
    int k;

    word_size = n*2;
    last_word = 0;

    if(word_size>=sizeof(word_t)*4) {
        if(verbose) fprintf(logfile,"Max word length exceeded (%i), exiting\n",word_size);
        exit(1);
    }

    for(k=0;k<word_size;k++) {
        last_word <<= 2;
        last_word |= 3;
    }

    half_size = n;
    gap_size = g;
    full_size = 2*n+g;

    lower_mask = 1;
    for(int i=0;i<half_size;i++) lower_mask<<=2;
    lower_mask--;

    upper_mask = last_word ^ lower_mask;

    index = (index_t*) malloc((size_t)sizeof(index_t)*(last_word + ARRAY_MARGIN));
    count = (index_t*) malloc((size_t)sizeof(index_t)*(last_word + ARRAY_MARGIN));

    if(index==0 ||count==0) {
        if(verbose) fprintf(logfile,"Not enough memory to create index (last_word=%lx), exiting\n",(long)last_word);
        exit(1);
    }

    for(i=0;i<=last_word;i++) index[i] = count[i] = 0;

    max_key = 0;

    phase = PHASE_READY1;

}

template<class ordered_set>
void dictionary<ordered_set>::describe() {
    if(!verbose) return;
    fprintf(logfile,"-- gapped dictionary --\n");
    fprintf(logfile,"phase:\t%s\n",PHASE[phase]);
    fprintf(logfile,"word:\t%i\ngap:\t%i\nunit:\t%i\n",word_size,gap_size,(int)sizeof(ordered_set));
    fprintf(logfile,"index:\t%li entries\n\t%li bytes\n",(long)last_word+1,(long)sizeof(index_t)*(last_word+1));
    fprintf(logfile,"size:\t%li entries\n\t%li bytes\n",(long)index[last_word+1],(long)(sizeof(ordered_set)*index[last_word+1]));
    fprintf(logfile,"-----------------------\n\n");
}

template<class ordered_set>
void dictionary<ordered_set>::info() {
    if(verbose) fprintf(logfile," %li/%li",(long)last_word,(long)index[last_word+1]);
}

template<class ordered_set>
void dictionary<ordered_set>::check() {
    bool flag = true;
    word_t a;
    index_t i;

    if(phase<PHASE_COMPLETED) if(verbose) fprintf(logfile,"Warning: too early to check!");

    for(a=0;a<=last_word;a++) {
	for(i=index[a];i<index[a+1]-1;i++) {
	    if(table[i+1] < table[i]) {
		if(verbose) {
		    fprintf(logfile,"\nFound inverstion at ");
           	    fprint_word(logfile,a);
		    fprintf(logfile,"\t");
            	    table[i].fprint(logfile);
            	    fprintf(logfile," ");
            	    table[i+1].fprint(logfile);
		    fprintf(logfile,", exiting\n");
            	    exit(1);
		}
	    }
	}
    }
    if(verbose) fprintf(logfile," (c)");
}

/***********************************************************************************************************/
template<class ordered_set>
void dictionary<ordered_set>::read_from_suf(char *filename, int reverse_complement, int length_limit) {
    read_from_suf(filename, NULL, reverse_complement, length_limit);
}

template<class ordered_set>
void dictionary<ordered_set>::read_from_suf(FILE *f, int reverse_complement, int length_limit) {
    read_from_suf(f, NULL, reverse_complement, length_limit);
}

template<class ordered_set>
void dictionary<ordered_set>::read_from_suf(char *filename, subset *sbs, int reverse_complement, int length_limit) {
    FILE *f;
    f = fopen(filename, "rb");
    if(f==NULL) {
        fprintf(logfile,"Cannot read from '%s', exiting\n", filename);
	exit(1);
    }
    read_from_suf(f, sbs, reverse_complement, length_limit);
    fclose(f);
}

template<class ordered_set>
void dictionary<ordered_set>::read_from_suf(FILE *f, subset *sbs, int reverse_complement, int length_limit) {
    char buff[MAXFILEBUFFLENGTH];
    char seq[MAXFILEBUFFLENGTH];
    keys_t key;

    fseek(f,SEEK_SET,0);
    if(verbose) fprintf(logfile," pass 1");
    while(fgets(buff,MAXFILEBUFFLENGTH,f)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%i %*s %*i %*i %*i %*i %s",&key,seq);
	if(sbs!=NULL) if(key>=sbs->length) continue;
	if(sbs!=NULL) if(sbs->index[key]==0) continue;
	if(reverse_complement) rev1(seq);
        pass1_add(seq, length_limit);
    }
    pass1_make();
    if(verbose) fprintf(logfile," ok");

    fseek(f,SEEK_SET,0);
    if(verbose) fprintf(logfile," pass 2");
    while(fgets(buff,MAXFILEBUFFLENGTH,f)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%i %*s %*i %*i %*i %*i %s",&key,seq);
        if(sbs!=NULL) if(key>=sbs->length) continue;
        if(sbs!=NULL) if(sbs->index[key]==0) continue;
	if(reverse_complement) rev1(seq);
        pass2_add(seq, key, length_limit);
    	if(key>max_key) max_key=key;
    }
    pass2_make();
    if(verbose) fprintf(logfile," ok");
    if(verbose) info();
}

/****/
template<class ordered_set>
void dictionary<ordered_set>::read_from_muf(char *filename, dictionary **dict, subset *sbs, char **species_name, int n, int reverse_complement, int length_limit) {
    FILE *f;
    f = fopen(filename, "rb");
    if(f==NULL) {
        fprintf(logfile,"Cannot read from '%s', exiting\n", filename);
        exit(1);
    }
    read_from_muf(f, dict, sbs, species_name, n, reverse_complement, length_limit);
    fclose(f);
}


template<class ordered_set>
void dictionary<ordered_set>::read_from_muf(FILE *f, dictionary **dict, subset *sbs, char **species_name, int n, int reverse_complement, int length_limit) {
    char buff[MAXFILEBUFFLENGTH];
    char seq[MAXFILEBUFFLENGTH];
    char name[MAXFILEBUFFLENGTH];
    keys_t key;
    int i,j;

    fseek(f,SEEK_SET,0);
    if(verbose) fprintf(logfile," pass 1");
    while(fgets(buff,MAXFILEBUFFLENGTH,f)) {
	if(strlen(buff)<2) break;
        if(strcmp(buff,"a id=")>0) {
            sscanf(buff+5,"%li",(long*)&key);
	    if(sbs!=NULL) if(key>=sbs->length) continue;
	    if(sbs!=NULL) if(sbs->index[key]==0) continue;
            while(fgets(buff,MAXFILEBUFFLENGTH,f)) {
		if(strlen(buff)<2) break;
                sscanf(buff,"%*c %s %*i %*i %*i %*i %s",&name[0],&seq[0]);
		for(j=0;name[j]!=0 && name[j]!='.';j++); name[j]=0;
                for(i=0;i<n;i++) {
                    if(strcmp(name,species_name[i])==0) {
			if(reverse_complement) rev1(seq);
                        dict[i]->pass1_add(seq, length_limit);
                    }
                }
            }
        }
    }

    for(i=0;i<n;i++) dict[i]->pass1_make();
    if(verbose) fprintf(logfile," ok");

    fseek(f,SEEK_SET,0);
    if(verbose) fprintf(logfile," pass 2");
    while(fgets(buff,MAXFILEBUFFLENGTH,f)) {
	if(strlen(buff)<2) break;
        if(strcmp(buff,"a id=")>0) {
            sscanf(buff+5,"%li",(long*)&key);
            if(sbs!=NULL) if(key>=sbs->length) continue;
            if(sbs!=NULL) if(sbs->index[key]==0) continue;
            while(fgets(buff,MAXFILEBUFFLENGTH,f)) {
		if(strlen(buff)<2) break;
                sscanf(buff,"%*c %s %*i %*i %*i %*i %s",&name[0],&seq[0]);
	   	for(j=0;name[j]!=0 && name[j]!='.';j++); name[j]=0;
                for(i=0;i<n;i++) {
                    if(strcmp(name,species_name[i])==0) {
			if(reverse_complement) rev1(seq);
                  	dict[i]->pass2_add(seq,key, length_limit);
                  	if(key > dict[i]->max_key) dict[i]->max_key=key;
                    }
                }
            }
        }
    }
    for(i=0;i<n;i++) dict[i]->pass2_make();
    if(verbose) fprintf(logfile," ok");
}

template<class ordered_set>
void dictionary<ordered_set>::subset_on(subset* restriction) {
    word_t a;
    index_t i,j,s;
    keys_t k;
    j=s=0;
    for(a=0;a<=last_word;a++) {
	for(i=index[a]; i<index[a+1]; i++) {
	    k = table[i].getid();
	    table[j] = table[i];
	    if(restriction->index[k]>0) j++;
	}
	index[a]=s;
	s=j;
    }
    index[a]=s;
    if(verbose) fprintf(logfile,"->");
    info();
}

/************************************************************************************************************/

template<class ordered_set>
void dictionary<ordered_set>::intersect_many(dictionary **dict, int n, double *weight, double threshold) {
    intersect_many(dict, n, weight, threshold,(char*)"");
}

template<class ordered_set>
void dictionary<ordered_set>::intersect_many(dictionary **dict, int n, double *weight, double threshold, char *message) {
    int i;
    word_t a;
    index_t *ptr,*qtr,ppp;
    bool flag;
    double score;
    int k;

    ptr = (index_t*) malloc((size_t)sizeof(index_t)*(n + ARRAY_MARGIN));
    qtr = (index_t*) malloc((size_t)sizeof(index_t)*(n + ARRAY_MARGIN));

    for(i=0;i<n;i++) qtr[i]=ptr[i]=0;

    for(a=0;a<=dict[0]->last_word;a++) {
	progressbar(a,dict[0]->last_word,message);
        flag=true;
        while(flag) {
            flag = false;
	    /** check if at least one pointer is within range **/
            for(i=0;i<n;i++) {
		if(ptr[i] < dict[i]->index[a+1]) { 
		    k = i;
		    flag = true;
		}
	    }
	    if(flag) {
                for(i=0;i<n;i++) {
		    if(ptr[i] < dict[i]->index[a+1]) 
			if(dict[i]->table[ptr[i]] < dict[k]->table[ptr[k]]) k = i;
		}
		score=0;
		for(i=0;i<n;i++) {
		    if(ptr[i]<dict[i]->index[a+1]) 
			if(dict[i]->table[ptr[i]].equiv(dict[k]->table[ptr[k]])) 
			    score+=weight[i];
		}
                if(score>=threshold) {
                    for(i=0;i<n;i++) {
			if(ptr[i]<dict[i]->index[a+1]) 
			    if(dict[i]->table[ptr[i]].equiv(dict[k]->table[ptr[k]])) 
				dict[i]->table[qtr[i]++] = dict[i]->table[ptr[i]];
		    }
                }
		ppp = ptr[k];
                for(i=0;i<n;i++) {
		    if(ptr[i]<dict[i]->index[a+1])
			if(dict[i]->table[ptr[i]].equiv(dict[k]->table[ppp])) ptr[i]++;
		}
            }
        }
        for(i=0;i<n;i++) dict[i]->index[a+1]=qtr[i];
    }
    for(i=0;i<n;i++) {
	dict[i]->index[a+1]=qtr[i];
	dict[i]->table = (ordered_set*)realloc(dict[i]->table,(size_t)sizeof(ordered_set)*(qtr[i] + ARRAY_MARGIN));
    }
    free(ptr);
    free(qtr);
}

template<class ordered_set>
void dictionary<ordered_set>::foldall(FILE *f, dictionary **dict, char ** name, int n, double *weight, double threshold, conservation_table *cons_f, conservation_table *cons_r, double min_cons, int report_length) {
    char *buff;
    char aux[MAXBUFFLENGTH];
    int i;
    word_t a;
    bool flag=true;
    index_t *ptr, *qtr, ppp;
    index_t j;
    word_t **word;
    double score,cons;
    int length;
    int count=0;
    int k;

    buff = (char*) malloc(sizeof(char)*(MAXBUFFLENGTH * n + ARRAY_MARGIN));

    word = (word_t**)malloc((size_t)sizeof(word_t*)*(n + ARRAY_MARGIN));
    ptr = (index_t*) malloc((size_t)sizeof(index_t)*(n + ARRAY_MARGIN));
    qtr = (index_t*) malloc((size_t)sizeof(index_t)*(n + ARRAY_MARGIN));

    for(i=0;i<n;i++) {
	if(verbose) fprintf(logfile, "[Pre-fold %s", name[i]);
	word[i] = (word_t*)malloc(sizeof(word_t)*(dict[i]->index[dict[i]->last_word + 1] + ARRAY_MARGIN));
	if(word[i] == NULL) {
	    fprintf(stderr,"ERROR: not enough memory, exiting\n");
	    exit(1);
	}
	for(a=0;a<=dict[i]->last_word;a++) for(j=dict[i]->index[a];j<dict[i]->index[a+1];j++) word[i][j] = a;
	ordered_set::quicksort_leid(dict[i]->table, word[i], 0, dict[i]->index[dict[i]->last_word + 1] - 1); // in each specie the table with pairs is sorted by (left.id, right.id) 
	ptr[i] = 0; // pointer to the current element in species i
	if(verbose) dict[i]->info();
    	if(verbose) fprintf(logfile, "]\n");
    }

    flag=true;
    while(flag) {
	flag = false; 
        for(i=0;i<n;i++) {
	    a = dict[i]->last_word;
            if(ptr[i] < dict[i]->index[a+1]) {
                k = i;		// k is the species where ptr points to an existing element 
                flag = true;	// flag is up if there is once such species
            }
        }
	if(flag) { // if at least one ptr points to an existing element
	    progressbar(ptr[k],dict[k]->index[dict[k]->last_word + 1], (char*) "Joint fold\t");
            for(i=0;i<n;i++) { // by the end of this loop k will be the specie with the smallest current element pointed by ptr
		a = dict[i]->last_word;
                if(ptr[i] < dict[i]->index[a+1])
                    if(dict[i]->table[ptr[i]].leid(dict[k]->table[ptr[k]])) k = i;
            }
	    ppp = ptr[k]; // keep the pointer because ptr will change later

            score = 0;
	    cons  = 0;
	    buff[0] = 0;
            for(i=0;i<n;i++) { // identify species in which the pairid of the current pointer is the same as in specie k and ppp; look up same segment by qtr and fold
		a = dict[i]->last_word;
                if(ptr[i]<dict[i]->index[a+1]) {
                    if(dict[i]->table[ptr[i]].left.getid() == dict[k]->table[ppp].left.getid() && dict[i]->table[ptr[i]].right.getid() == dict[k]->table[ppp].right.getid()) {
			for(qtr[i]=ptr[i];qtr[i]<dict[i]->index[a+1] && dict[i]->table[ptr[i]].left.getid()  == dict[i]->table[qtr[i]].left.getid() &&
				dict[i]->table[ptr[i]].right.getid() == dict[i]->table[qtr[i]].right.getid();qtr[i]++); // this loop is to find all entries ptr[i] to qtr[i]
			dict[i]->fold_simple(aux, word[i], i, ptr[i], qtr[i], cons_f, cons_r, &length);
			if(length>=report_length) {
			    strcat(buff,aux);
			    score += weight[i];
			    cons  += weight[i]*length*(cons_f->rep_log(dict[i]->table[ptr[i]].left.getid()) + cons_r->rep_log(dict[i]->table[ptr[i]].right.getid()))/dict[i]->word_size;
			}
			ptr[i] = qtr[i];
		    }
		}
	    }
	    if(score>=threshold && cons>=min_cons) {
		fprintf(f, "%s", buff);
		count++;
	    } 
	}
	else {
	   progressbar(1,1,(char*) "Joint fold\t");
	}
    }
    fprintf(stdout, "%i structures found\n", count);

    for(i=0;i<n;i++) free(word[i]);
    free(word);
    free(ptr);
    free(qtr);
    free(buff);

}

template<class ordered_set>
int dictionary<ordered_set>::fold_and_save_standalone(FILE *f, int species, conservation_table *cons_f, conservation_table *cons_r, int report_length) {
    char buff[MAXBUFFLENGTH];
    index_t i,j;
    word_t a;
    int count=0;
    word_t *word;
    int length;

    if(verbose) fprintf(logfile, " sort");
    word = (word_t*)malloc(sizeof(word_t)*(index[last_word+1] + ARRAY_MARGIN));
    for(a=0;a<=last_word;a++) for(i=index[a];i<index[a+1];i++) word[i] = a;
    ordered_set::quicksort_leid(table, word, 0, index[last_word+1] - 1); // in each specie the table with pairs is sorted by (left.id, right.id)

    if(verbose) fprintf(logfile, ", fold");
    for(i=0;i<index[last_word+1];i=j) {
        for(j=i;j<index[last_word+1] && table[i].left.id==table[j].left.id && table[i].right.id==table[j].right.id;j++);
	fold_simple(buff, word, species, i, j, cons_f, cons_r, &length);
	if(length>=report_length) {
	    fprintf(f,"%s",buff);
	    count++;
	}
    }
    free(word);
    return(count);
}

template<class ordered_set>
void dictionary<ordered_set>::finalize() {
    index_t i=0;
}

/**********************************************************************************************************************************************************/

#define FOLD 100
/** This part is to be modified to get the structure from WORDS not base pairs **/

template<class ordered_set>
void dictionary<ordered_set>::fold_simple(char *buff, word_t* word, int species, index_t i, index_t j, conservation_table* cons_f, conservation_table *cons_r, int *length) {
/*  Save output into string buffer 'buff'
 *  the sorted dictionary is the paernt object; the corresponding word array is 'word'
 *  look between indices [i,j) in the table object of the dictionary
 *  conservation tables, weights of species are given
 *  threshold is the weight threshold; report_length is the threshold on length; min_cons is the threshold on sum cons score */

    index_t d;
    char arm5[MAXFILEBUFFLENGTH];
    char arm3[MAXFILEBUFFLENGTH];
    word_t a;
    int max_m,z;
    int k,l,m,n,q;
    int max_s,p;
    pos_t lp,rp;
    basepair *bp;

    bp = (basepair*)malloc(sizeof(basepair)*((j - i + 1) * half_size * 2 + ARRAY_MARGIN));
    if(bp==NULL) {
        if(verbose) fprintf(logfile, "Not enough memory to store base pairs, exiting\n");
        exit(1);
    }

    /** melt by single basepair **/
    n = 0;
    for(d=i;d<j;d++) {
        a = word[d];
        lp = table[d].left.getpos();
        rp = table[d].right.getpos() - 2*half_size - table[d].right.getgap() + 1;
        for(z=0;z<half_size;z++,a>>=2) bp[n++].set(lp--, rp++, a & 3);
        lp-=table[d].left.getgap();
        rp+=table[d].right.getgap();
        for(z=0;z<half_size;z++,a>>=2) bp[n++].set(lp--, rp++, a & 3);
    }

    /** sort to prepare dp table **/
    basepair::qs(bp,0,n-1);

    /** remove redundant base pairs **/
    basepair::compress(bp,(index_t*)&n);

    /** perform Nussinov-dp; can possibly be modified to use stacking energies + average special internal loops **/
    max_s = -INFTY;
    for(m = n-1; m>=0; m--) {
        bp[m].score = 1;
        bp[m].next = -INFTY;
        for(k = m + 1;k < n;k++) {
            if(bp[m].pos5 == bp[k].pos5) continue;
            if(bp[k].pos5 - bp[m].pos5 - 1 > gap_size) break;
            p = 0; //abs(bp[k].pos5 - bp[m].pos5 - 1) + abs(bp[k].pos3 - bp[m].pos3 - 1);
            if(bp[m].pos3 > bp[k].pos3 && bp[m].pos3 - bp[k].pos3 - 1 <= gap_size) {
                if(bp[k].score + FOLD - p > bp[m].score) {
                    bp[m].score = bp[k].score + FOLD - p;
                    bp[m].next = k;
                }
            }
        }
        if(bp[m].score >= max_s) {
            max_s = bp[m].score;
            max_m = m;
        }
    }

    /** Traceback the best structure **/
    if(max_s<2*half_size) {
        if(verbose) fprintf(logfile,"WARNING: Suspiciously few base pairings\n");
    }
    k = l = 0;
    for(m = max_m; m>=0; m = bp[m].next) {
        arm5[k++] = decode1(bp[m].base5);
        arm3[l++] = decode1(3 - bp[m].base5);
        if(bp[m].next>=0) {
            for(q = bp[m].pos5; q < bp[bp[m].next].pos5 - 1; q++) arm5[k++] = '.';
            for(q = bp[m].pos3; q > bp[bp[m].next].pos3 + 1; q--) arm3[l++] = '.';
        }
        n = m;
    }
    arm5[k] = arm3[l] = 0;
    mirr(arm3);

    *length = (int)ceil((double)max_s/FOLD);
    //sprintf(buff, "%i\t%li\t%li\t%i\t%i\t%i\t%2.2lf\t%s\t%i\t%i\t%2.2lf\t%s\n", species, (long)table[i].left.id, (long)table[i].right.id, (*length),
    //	bp[max_m].pos5, bp[n].pos5, cons_f->rep_log(table[i].left.id), arm5, bp[n].pos3, bp[max_m].pos3, cons_r->rep_log(table[i].right.id), arm3);
    sprintf(buff, "%i\t%li\t%li\t%i\t%i\t%i\t%i\t%i\t%2.2lf\t%2.2lf\t%s\t%s\n",species, (long)table[i].left.id, (long)table[i].right.id, (*length),
	bp[max_m].pos5, bp[n].pos5, bp[n].pos3, bp[max_m].pos3, cons_f->rep_log(table[i].left.id), cons_r->rep_log(table[i].right.id), arm5, arm3);
    free(bp);
}


/***********************************************************************************************************/
template<class ordered_set>
void dictionary<ordered_set>::var_ac_sub(word_t w, word_t a, word_t c, word_t g, word_t t, int pos, int max_mutations, word_t *result, int *n_results) {
    word_t a0 = a; //nucleotide 'a' at position pos
    word_t c0 = c; //nucleotide 'c' at position pos
    word_t g0 = g; //nucleotide 'g' at position pos
    word_t t0 = t; //nucleotide 't' at position pos
    if(pos>=word_size || max_mutations==0) {
        if(*n_results < GT_ARRAY_CAPACITY) result[(*n_results)++]=w;
        return;
    }
    a<<=2;
    c<<=2;
    g<<=2;
    t<<=2;
    if((w&t0)==a0) { var_ac_sub(w&(~t0)|g0, a, c, g, t, pos + 1, max_mutations - 1, result, n_results); } // if at position pos of w it is 'a' then change 'a' to 'g' and proceed recursively
    if((w&t0)==c0) { var_ac_sub(w&(~t0)|t0, a, c, g, t, pos + 1, max_mutations - 1, result, n_results); } // if at position pos of w it is 'c' then change 'c' to 't' and proceed recursively
    var_ac_sub(w, a, c, g, t, pos+1, max_mutations, result, n_results);
    return;
}

template<class ordered_set>
void dictionary<ordered_set>::revcomp_gt(word_t w, int n_gt_pairs, word_t *result, int *n_results) {
// takes the word w and tabulates all the complementary words with not more than n_gt_pairs GT base pairs
    word_t rc = reverse_word(w);
    rc = rc ^ last_word;
    *n_results=0;
    var_ac_sub(rc, 0, 1, 2, 3, 0, n_gt_pairs, result, n_results);
}




/***********************************************************************************************************/

template<class ordered_set>
void dictionary<ordered_set>::pass1_add(char* str, int length_limit) {
    int b;
    int i,j,l,k;
    bool flag=false;

    word_t y;
    word_t a=0;

    l = strlen(str);
    i = INFTY;
    length_limit/=2;

    for(j=0;j<l;j++) {
	if(j<length_limit && j<=l-length_limit) {
	    j = l-length_limit;
	    flag = false;
	}
	b = code1(str[j]);
        if(b<4) {
            if(!flag) {
                i=j;			// i is the first position of the continuous (uninterrupted) stretch
                flag=true;		// j is the current position
            }
	    a<<=2;
            a|=b;
	    k = j - i + 1 - word_size;	// k is the max possible gap size given the values of i and j
	    if(k>gap_size) k=gap_size;
	    for(;k>=0;k--) {
		y = (a & lower_mask) | ( (a >> k >> k) & upper_mask);
		if(passfilter(y)) { 
		    index[y]++;
		}
	    }
	}
	else {
	    flag=false;
	}
    }
}

template<class ordered_set>
void dictionary<ordered_set>::pass2_add(char* str, keys_t key, int length_limit) {
    int b;
    int i,j,l,k;
    bool flag=false;

    word_t y;
    word_t a=0;

    l = strlen(str);
    i = INFTY;
    length_limit/=2;

    for(j=0;j<l;j++) {
        if(j<length_limit && j<=l-length_limit) {
            j = l-length_limit;
            flag = false;
        }
        b = code1(str[j]);
        if(b<4) {
            if(!flag) {
                i=j;
                flag=true;
            }
            a<<=2;
            a|=b;
            k = j - i + 1 - word_size;
            if(k>gap_size) k=gap_size;
            for(;k>=0;k--) {
                y = (a & lower_mask) | ( (a >> k >> k) & upper_mask);
                if(passfilter(y)) { 
		    table[index[y]+count[y]].set(key, j, k);
                    count[y]++;
		}
            }
        }
        else {
            flag=false;
        }
    }
}

template<class ordered_set>
void dictionary<ordered_set>::pass1_make() {
    word_t i;
    index_t a;
    index_t s=0;

    for(i=0;i<=last_word;i++) {
        a = index[i];
        index[i]=s;
        s+=a;
    }
    index[i]=s;

    table = (ordered_set*)malloc((size_t)sizeof(ordered_set)*(s + ARRAY_MARGIN));
    if(table==NULL) {
        if(verbose) fprintf(logfile,"Not enough memory to create content tables (last entry=%li), exiting\n",(long)s);
        exit(1);
    }
    phase = PHASE_READY2;
}

template<class ordered_set>
void dictionary<ordered_set>::pass2_make() {
    free(count);
    count = NULL;
    phase = PHASE_COMPLETED;
}

template<class ordered_set>
void dictionary<ordered_set>::drop() {
    free(index);
    if(phase>PHASE_READY2) {
    	free(table);
    }
    word_size = 0;
    last_word = 0;
    phase = PHASE_CONTAINERONLY;
}

/****************************************************************************************************************************/

template<class ordered_set>
word_t dictionary<ordered_set>::reverse_word(word_t w) {
    word_t x=w;
    word_t y=0;
    for(int k=0;k<word_size;k++) {
        y<<=2;
        y|=(x & 3);
        x>>=2;
    }
    return(y);
}

template<class ordered_set>
void dictionary<ordered_set>::print_word(word_t w) {
    fprint_word(stdout,w);
}

template<class ordered_set>
void dictionary<ordered_set>::fprint_word(FILE *f, word_t w) {
    word_t x = reverse_word(w);
    for(int k=0;k<word_size;k++) {
        fprintf(f,"%c",decode1(x&3));
        x>>=2;
    }
}

template<class ordered_set>
void dictionary<ordered_set>::fprint(FILE *f) {
    word_t i;
    index_t k;
    for(i=0;i<=last_word;i++) {
        if(index[i]<index[i+1]) {
            fprint_word(f,i);
            for(k=index[i];k<index[i+1];k++) {
		table[k].fprint(f);
		fprintf(f," ");
            }
            fprintf(f,"\n");           
        }
    }
}

template<class ordered_set>
void dictionary<ordered_set>::print() {
    fprint(logfile);
}

/************************************************************************************/
template<class ordered_set>
void dictionary<ordered_set>::save(FILE *f) {
    int dict_type = DICT_VERSION;
    fwrite(&dict_type, sizeof(int), 1, f);
    fwrite(&word_size, sizeof(int), 1, f);
    fwrite(&half_size, sizeof(int), 1, f);
    fwrite(&full_size, sizeof(int), 1, f);
    fwrite(&gap_size,  sizeof(int), 1, f);
    fwrite(&upper_mask, sizeof(word_t), 1, f);
    fwrite(&lower_mask, sizeof(word_t), 1, f);
    fwrite(&last_word, sizeof(word_t), 1, f);
    fwrite(index, sizeof(index_t), last_word + FILE_MARGIN,f);
    fwrite(table, sizeof(ordered_set),  index[last_word+1] + FILE_MARGIN, f);
}

template<class ordered_set>
void dictionary<ordered_set>::load(FILE *f) {
    int dict_type;
    fread(&dict_type, sizeof(int), 1, f);
    if(dict_type!=DICT_VERSION) {
        if(verbose) fprintf(logfile,"Dictionary cannot be loaded within gapped class, exiting\n");
        exit(1);
    }
    fread(&word_size, sizeof(int), 1, f);
    fread(&half_size, sizeof(int), 1, f);
    fread(&full_size, sizeof(int), 1, f);
    fread(&gap_size,  sizeof(int), 1, f);
    fread(&upper_mask, sizeof(word_t), 1, f);
    fread(&lower_mask, sizeof(word_t), 1, f);
    fread(&last_word, sizeof(word_t), 1, f);
    index = (index_t*) malloc((size_t)sizeof(index_t)*(last_word + ARRAY_MARGIN));
    if(index==NULL) {
        if(verbose) fprintf(logfile,"Not enough memory to create index (last_word=%lx), exiting\n",(long)last_word);
        exit(1);
    }
    fread(index, sizeof(index_t), last_word + FILE_MARGIN,f);
    table = (ordered_set*) malloc((size_t)sizeof(ordered_set)*(index[last_word+1] + ARRAY_MARGIN));
    if(table==NULL) {
        if(verbose) fprintf(logfile,"Not enough memory to create content tables (last entry=%li), exiting\n",(long)index[last_word+1]);
        exit(1);
    }
    fread(table, sizeof(ordered_set), index[last_word+1] + FILE_MARGIN, f);
    phase = PHASE_COMPLETED;
}

template<class ordered_set>
void dictionary<ordered_set>::save(char *filename) {
    FILE *f;
    f = fopen(filename, "wb");
    if(f==NULL) {
        if(verbose) fprintf(logfile,"Cannot write to '%s', exiting\n", filename);
	exit(1);
    }
    save(f);
    fclose(f);
}

template<class ordered_set>
void dictionary<ordered_set>::load(char *filename) {
    FILE *f;
    f = fopen(filename, "rb");
    if(f==NULL) {
        if(verbose) fprintf(logfile,"Cannot read from '%s', exiting\n", filename);
	exit(1);
    }
    load(f);
    fclose(f);
}

/*************************************************************************************************/
template<class ordered_set>
int* dictionary<ordered_set>::fill_cons() {
    index_t i;
    keys_t k;
    int *res;
    res = (int*)malloc(sizeof(int)*(max_key + ARRAY_MARGIN));
    for(k=0;k<=max_key;k++) res[k] = 0;
    for(i=0;i<index[last_word+1];i++) {
	k = table[i].getid();
	if(k<=max_key) res[k]++; else fprintf(stderr,"!");
    }
    return(res);
}
/*************************************************************************************************/
template<class ordered_set>
void dictionary<ordered_set>::mask_low_complexity(int n) {
    int n_times = full_size / n;
    word_t a, b, w;
    int i,k;

    if(full_size % n != 0) {
        if(verbose) fprintf(logfile,"[WARNING: low complexity incomplete or redundant: %i is not a multiple of %i]\n",full_size, n);
        n_times++;
    }

    w = 1<<(n<<1);
    for(a=0;a<w;a++) {
        for(b=0,i=0;i<n_times;i++) b=(b<<(n<<1))|a;
        for(k=0;k<gap_size;k++) {
            index[((b>>(2*k)) & upper_mask) | (b & lower_mask)]=1;

        }
    }
    free(count);
    if(verbose) fprintf(logfile,"[Created low complexity mask: %i-nt repeats excluded]\n", n);
}

template<class ordered_set>
void dictionary<ordered_set>::mask_low_GCcontent(int l) {
    word_t a,b;
    int k,s;

    for(a=0;a<=last_word;a++) {
        b=a;s=0;
        for(k=0;k<word_size;k++) {
            if((b&3)==1||(b&3)==2) s++;
            b>>=2;
        }
        if(s<l) index[a]=1;
    }
    if(verbose) fprintf(logfile,"[Created min GC content mask: at least %i GC base pair per %i-mer]\n", l, word_size);
}

template<class ordered_set>
bool dictionary<ordered_set>::passfilter(word_t y) {
    return(filter->index[y]==0);
}

