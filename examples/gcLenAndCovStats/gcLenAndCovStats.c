/**
 * qseqToFastq
 *  Author: John St. John
 *  Date: 1/14/2010
 *  
 *  Qseq format; Tab seperated on a single line:
    0. Machine name (hopefully) unique
    1. Run number: (hopefully) unique
    2. lane number: [1..8]
    3. Tile number: positive integer
    4. X: x coordinate of the spot, (can be negative)
    5. Y: y coordinate of the spot, can be negative)
    6. Index: positive integer, should be greater than 0 (the files I see have it == to 0)
    7. Read number, 1 for single read, 1 or 2 for paired end
    8. sequence
    9. quality, the calibrated quality string.
    10. filter, did the read pass qc 0 - no, 1 - yes

    Also note that '.' in the sequence is equivalent to 'N'.
 *
 *
 */
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <limits.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "common.h"
#include "options.h"
#include "linefile.h"
#include "dnautil.h"
#include "fa.h"
#include "uthash.h"

/*
Global variables
 */
bool verboseOut = false;
int sampleCount = 0;
struct cov_gc_stats_hash {
    char name[25];                /* key (structure POINTS TO string */
    unsigned length;
    float gc;
    unsigned short *cov; /* initialize to the sequence length */
    UT_hash_handle hh;         /* makes this structure hashable */
};


struct cov_gc_stats_hash *chrInfo = NULL;

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "gcLenAndCovStats -- Outputs some coverage and gc statistics given a fasta and the output of samtools depth.\n"
      "usage: samtools depth -q [min bq] -Q [min mq] | gcLenAndCovStats [required options]\n"
      "\n**required** options:\n"
      "\t-fasta=FILE\tFile name holding the vcf file to parse.\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"fasta",OPTION_STRING},
    {NULL, 0}
}; //end options()


int initInfoHashFromFasta(char *fasta)
/**
 * Fills the hash table of chrom sizes
 *
 *
 *
 */
{
  DNA *seq;
  int length;
  char *name;
  struct cov_gc_stats_hash *s;
  FILE * fafile = fopen(fasta,"r");

  while(faFastReadNext(fafile, &seq, &length, &name)){
    s = malloc(sizeof(struct cov_gc_stats_hash));
    strcpy(s->name, name);
    s->length = length;
    s->cov = malloc(sizeof(unsigned short)*s->length);
    zeroBytes(s->cov,sizeof(unsigned short)*s->length);
    int histogram[4]; //a,c,g,t
    dnaBaseHistogram(seq, length, histogram);
    /*
               0 for t
               1 for c
               2 for a
               3 for g
               */
    int A = histogram[ntVal['A']];
    int C = histogram[ntVal['C']];
    int G = histogram[ntVal['G']];
    int T = histogram[ntVal['T']];
    s->gc = (((float)(G+C))/((float)(A+C+G+T)));
    HASH_ADD_STR(chrInfo,name,s);
  }
  return 0;
}

int fillCoverage(void)
/**
    seqid\tposition(1based)\tcoverage
 */
{
  struct lineFile *lf;
  lf = lineFileStdin(TRUE);
  char *line;
  struct cov_gc_stats_hash *s;
  char lastSeqId[200];
  strcpy(lastSeqId,"23121NoTASeqUeNCE22");


  while(lineFileNextReal(lf,&line)){
    // if line starts with "#" skip it
    char *split[3];
    int numRecords = chopByWhite(line,split,3);
    char *seqId = split[0];
    int pos0 = atoi(split[1])-1;
    int cov = atoi(split[2]);

    if(strcmp(seqId,lastSeqId)!= 0){
      //update the pointer to our chrom info struct if
      //we have a new chrom name in the line
      strcpy(lastSeqId,seqId);
      HASH_FIND_STR(chrInfo,seqId,s);
    }

    s->cov[pos0] = (unsigned short int)cov;

  }
  return 0;
}

int isort(const void *x, const void *y) {
  return (*(int*)x - *(int*)y);
}


float mean(unsigned short *arr, int length){
  int i;
  unsigned int sum;
  for(i=0;i<length;i++){
    sum += arr[i];
  }
  return (((float)sum)/((float)length));
}

float meanShorth(unsigned short *arr, int length){
  int i = 0;
  //get the starting position of the shorth
  int halfLen = length/2; //floor
  int minDiff = INT_MAX;
  int shorthStart = -1;

  //get the shorth info
  for(i=0;i<halfLen;i++){
    int begin = arr[i];
    int end = arr[i+halfLen-1];
    int diff = abs(end - begin);
    if(diff < minDiff){
      shorthStart = i;
      minDiff = diff;
    }
  }

  int upperBounds = shorthStart+halfLen;

  //return the mean of the half interval starting at
  //the shorth start
  return mean(arr+shorthStart,halfLen);
}

struct minMax{
  float min;
  float max;
};

struct minMax minMaxMeanWindow(unsigned short *arr, int windowSize, int length){
  struct minMax m;
  m.min = 0;
  m.max = 0;
  int i;
  int first = 0;
  int upper = length-(windowSize*2);
  float mymean;
  for(i=windowSize;i<upper;i++){
    mymean = mean(arr+i,windowSize);
    if(first == 0){
      //initialize
      m.max = mymean;
      m.min = mymean;
    }else{
    if(mymean > m.max)
      m.max = mymean;
    if(mymean < m.min)
      m.min = mymean;
    }
  }

  return m;
}


int printChromInfo(FILE *out){
  struct cov_gc_stats_hash *s,*tmp;
  HASH_ITER(hh, chrInfo, s, tmp) {
//      HASH_DEL(chrSizes,s);  /* delete; users advances to next */
//      free(s);            /* optional- if you want to free  */
    int i;

    //sort the array
    qsort(s->cov,s->length,sizeof(unsigned short),isort);
    float mymean = mean(s->cov,s->length);
    float mymeanshorth = meanShorth(s->cov,s->length);
    struct minMax mm = minMaxMeanWindow(s->cov, 50, s->length);

    fprintf(out,"%s\t%d\t%f\t%f\t%f\t%f\n",s->name,s->length,mymeanshorth,mymean,mm.min,mm.max);
  }
  return 0;
}


int main(int argc, char *argv[])
/* Process command line. */
{
  char *fasta = NULL;
  optionInit(&argc, argv, options);

  verboseOut = optionExists("verbose");
  fasta = optionVal("fasta",NULL);

  if (fasta == NULL){
    fprintf(stderr,"Error: Must provide fasta file\n");
    usage();
  }
  initInfoHashFromFasta(fasta);
  fillCoverage();
  printChromInfo(stdout);
  return 0;
} //end main()


