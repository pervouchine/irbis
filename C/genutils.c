#include "genutils.h"
#include <time.h>

int verbose = 1;
FILE* logfile = stderr;
int lowercase_allowed = 1;
char chr_names[MAXCHR][MAXCHRNAMELENGTH];
int N_CHR_NAMES=0;


/**********************************************************************/

int assign_code(char* name) {
    int i;

    for(i=0;i<N_CHR_NAMES;i++) if(strcmp(chr_names[i],name)==0) break;
    if(i<N_CHR_NAMES) return(i);
    strcpy(chr_names[i],name);
    N_CHR_NAMES++;
    if(N_CHR_NAMES>=MAXCHR) {
        if(verbose) fprintf(stderr,"Too many chromosome names, exiting\n");
        exit(1);
    }
    return(i);
}

int get_chr_code(char* name) {
    int i;

    for(i=0;i<N_CHR_NAMES;i++) if(strcmp(chr_names[i],name)==0) break;
    if(i<N_CHR_NAMES) return(i);
    return(-1);
}

char* get_chr_name(int code) {
    return(chr_names[code]);
}

/**********************************************************************/

void swapi(int *i, int *j) {
    int x;
    x = *i, *i = *j; *j = x;
}

void swapc(char *c, char *d) {
    char x;
    x = *c; *c = *d; *d = x;
}

void swaps(char **c, char **d) {
    char* x;
    x = *c; *c = *d; *d = x;
}

/**********************************************************************************************************************/
int code1(char c) {
    int res=8;
    if(c>='a' && c<='z') c = c + 'A' - 'a';

    if(c=='A') res=0;
    if(c=='C') res=1;
    if(c=='G') res=2;
    if(c=='T') res=3;
    return(res);
}

int code2(char c, char cm) {
    int res;
    res = code1(c);
    if(cm=='n' || cm=='N') res|=4;
    if(c!=cm && (res&4)==0 && cm>0) fprintf(stderr,"[%c%c]",c,cm);
    return(res);
}

char decode1(int i) {
    if(i&8) return('N');
    if((i&3)==0) return('A');
    if((i&3)==1) return('C');
    if((i&3)==2) return('G');
    if((i&3)==3) return('T');
}

void decode2(int i, char *c, char *cm) {
    *cm = *c = decode1(i);
    if(i&4) {
        *cm = 'N';
        *c=tolower(*c);
    }
}

/**********************************************************************************************************************/

char compl1(char c) {
    if(c=='A') return('T');
    if(c=='C') return('G');
    if(c=='G') return('C');
    if(c=='T') return('A');

    if(c=='a') return('t');
    if(c=='c') return('g');
    if(c=='g') return('c');
    if(c=='t') return('a');

    if(c=='n'||c=='N'||c=='.') return(c);
    return('?');
}

void rev1(char *str) {
   int i,n,m;
   char a,b;
   char cm;

   n = strlen(str);
   if(n<2) return;
   for(i=0;i<(n+1)/2;i++) {
        a=compl1(str[i]);
        b=compl1(str[n-1-i]);
        str[i]=b;
        str[n-1-i]=a;
   }
}

void rev2(char *str, char *strm) {
   int i,n,m;
   char a,b,am,bm;

   n = strlen(str);
   if(n<2) return;
   for(i=0;i<(n+1)/2;i++) {
        a = compl1(str[i]);
        b = compl1(str[n-1-i]);
        str[i]=b;
        str[n-1-i]=a;
        am = compl1(strm[i]);
        bm = compl1(strm[n-1-i]);
        strm[i]=bm;
        strm[n-1-i]=am;
   }
}

int is_all_n(char *str) {
    int i,s,n;
    n = strlen(str);
    for(s=0,i=0;i<n;i++) if(str[i]=='n' || str[i]=='N') s++;
    return(s==n);
}


/**********************************************************************************************************************/

void mirr(char *str) {
    int i,n;
    char a;
    n = strlen(str);
    if(n<2) return;
    for(i=0;i<(n+1)/2;i++) {
	a=str[i];
	str[i]=str[n-1-i];
	str[n-1-i]=a;
    }
}	

void endslashdir(char *name) {
    int n;
    n = strlen(name);
    if(n==0) return;
    while(n>0 && name[n-1]!='/') {
	name[n-1]=0;
	n--;
    }
}

/**********************************************************************************************************************/


int fcode_buffer;
int fcode_count;

void fcode_start(FILE *f) {
    fcode_buffer=0;
    fcode_count=0;
}

void fcode(char c, char cm, FILE *f) {
    int b;
    b = code2(c,cm);
    fcode_buffer = (fcode_buffer<<4) | b;
    fcode_count++;
    if(fcode_count==8) {
        fwrite(&fcode_buffer,sizeof(int),1,f);
	fcode_count=0;
    }
}

void fcode_stop(FILE *f) {
    if(fcode_count>0) {
	for(;fcode_count<8;fcode_count++) fcode_buffer <<= 4;
	fcode_count = 0;
    	fwrite(&fcode_buffer,sizeof(int),1,f);
    }
}

char fdecode(int d, long *counter, FILE *g, FILE *gm) {
    int n=0;
    char c,cm;
    int b;
    while((*counter)>0) {
        b = (d >> 28) & 15;
        d <<= 4;
        decode2(b,&c,&cm);
        fprintf(g,"%c",c);
	if(gm!=NULL) fprintf(gm,"%c",cm);
        n++;
	(*counter)--;
        if(n==8) break;
    }
}

void fget_segment(char *dest, FILE *f, long offset, long a, long count) {
    long i,n;
    int b, d;
    char c,cm;

    i = a / 8;
    n = a % 8;

    fseek(f,sizeof(int)*(offset+i),SEEK_SET);
    fread(&d,sizeof(int),1,f);

    d <<= (4*n);

    while(count>0) {
        b = (d >> 28) & 15;
        d <<= 4;
        decode2(b,&c,&cm);
        *(dest++) = c;
        n++;
        count--;
        if(n==8) {
            fread(&d,sizeof(int),1,f);
            n=0;
        }
    }
    *dest=0;
}

/***************************************************************************************/


void progressbar(unsigned long current, unsigned long last, char *inscription) {
    int i,l,k;
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    l = (int)(log(w.ws_col)/log(2));
    unsigned long m = (1 << l);

    if(last==0) return;

    if(logfile!=stderr || !verbose) return;

    if((m*(current-1))/last < (m*current)/last) {
	k = (m*current)/last;
	fprintf(stderr,"%c%s[",13,inscription);
	for(i=0;i<k;i++) fprintf(stderr,"=");
	if(k<m) fprintf(stderr,">");
	for(i++;i<m;i++) fprintf(stderr," ");
	fprintf(stderr,"] %2.1lf%% ",(double)100*current/last);
    }
    if(current==last) fprintf(stderr,"\n");

}

/***************************************/
time_t timestamp;
void timestamp_set() {
    timestamp = time(NULL);
}

void timestamp_report() {
    time_t current;
    current = time(NULL);
    if(verbose && current>timestamp) fprintf(logfile,"[Completed in %2.2lf seconds]\n",difftime(current,timestamp));
}

void replace_extention(char *name, char *ext) {
    int i;
    for(i=strlen(name)-1;i>=0 && name[i]!='.';i--);
    if(i>=0) {
        name[i]=0;
        strcat(name,ext);
    }
    else {
        fprintf(logfile,"Erroneous file name specified: %s, exiting\n",name);
        exit(1);
    }
}

char* getfullname(char *path, char *name, char *ext) {
    char *result;
    result = (char*)malloc(sizeof(char)*(strlen(path)+strlen(name)+strlen(ext)+3));
    strcpy(result, path);
    strcat(result, name);
    strcat(result, ext);
    return(result);
}

const char ynreply[2][4]={"no","yes"};
char* yesno(int i) {
    return((char*)ynreply[i>0? 1 : 0]);
}
