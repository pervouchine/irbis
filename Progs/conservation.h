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

#define CTFACTOR 10000
#define BUFFLEN 1000

class conservation_table {
public:
    keys_t max_key; 
    int *before;
    int *after;

    conservation_table() {
    }

    void printf() {
	int i;
	for(i=0;i<=max_key;i++) {
	    fprintf(stdout,"%i\t%lf\n",i,rep_lin(i));
	}
    };

    conservation_table(keys_t the_max_key) {
	max_key = the_max_key;
	before = NULL;
	after = NULL;
    };
    void save(FILE *f) {
	int i,j;
	for(i=0;i<=max_key;i++) after[i] = before[i]>0 ? CTFACTOR*after[i]/before[i] : 0;
        /*fwrite(&max_key, sizeof(int), 1, f);
    	fwrite(after,   sizeof(int), max_key + FILE_MARGIN, f);*/
	fprintf(f, "L\t%i\n", max_key);
	for(i=0;i<=max_key;i++) fprintf(f, "%i\t%i\n",i, after[i]);
    }
    void load(FILE *f) {
	before = NULL;
	char buff[BUFFLEN];
	int i, v;
	/*fread(&max_key, sizeof(int), 1, f);
	after = (int*)malloc(sizeof(int)*(max_key + ARRAY_MARGIN));
	fread(after,   sizeof(int), max_key + FILE_MARGIN, f); */
	fgets(buff,BUFFLEN,f);
	sscanf(buff,"%*s %i", &max_key);
	after = (int*)malloc(sizeof(int)*(max_key + ARRAY_MARGIN));
	while(fgets(buff,BUFFLEN,f)) {
            if(strlen(buff)<2) break;
            sscanf(buff,"%i %i",&i, &v);
	    after[i] = v;    
	}
    }
    double rep_lin(keys_t key) {
	//if(key>max_key) return(0);
	return((double)after[key]/CTFACTOR);
    }
    
    double rep_log(keys_t key) {
	//if(key>max_key) return(0);
	if(after[key]==0) return(10);
	return((double)log(CTFACTOR/after[key])/log(10));
    }
};
