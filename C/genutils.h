#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/ioctl.h>
#include <string.h>
#include <iostream>
#include <locale>

#define MAXCHR 100000
#define MAXCHRNAMELENGTH 64
#define MAX2(A,B) ((A)<(B) ? (B) : (A))
#define MIN2(A,B) ((A)<(B) ? (A) : (B))

#define ARRAY_MARGIN 2 
#define FILE_MARGIN 2
#define INFTY 65535

#define MAXLONGBUFFLENGTH 400000
#define MAXBUFFLENGTH 1024

#define OUTPUT_TABSQ 0
#define OUTPUT_FASTA 1
#define OUTPUT_TBMAF 2

#define MAXSPECIES 64
#define ATOMIC LOS3A
#define GT_ARRAY_CAPACITY 256

extern int verbose;
extern FILE* logfile;
extern int N_CHR_NAMES;

int assign_code(char*);
int get_chr_code(char*);
char* get_chr_name(int);

void endslashdir(char*);
void mirr(char*);
void rev1(char*);
void rev2(char*, char*);
int is_all_n(char*);

int code1(char);
char decode1(int);

void fcode_start(FILE *f);
void fcode(char c, char cm, FILE *f);
void fcode_stop(FILE *f);
char fdecode(int i, long *counter, FILE *g, FILE *gm);

void swapi(int*, int*);
void swapc(char*, char*);
void swaps(char**, char**);

void progressbar(unsigned long, unsigned long, char*);
void timestamp_set();
void timestamp_report();

void replace_extention(char*, char*);
char* getfullname(char*, char*, char*);
char* yesno(int i);

void fget_segment(char *dest, FILE *f, long offset, long a, long count);

void swapi(int *i, int *j);
void swapc(char *c, char *d);
