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
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "common.h"
#include "options.h"
#include "linefile.h"
#include "fastq.h"



/*
Global variables
 */
#define DEF_MIN_LEN (30)


void usage()
/* Explain usage and exit. */
{
  errAbort(
      "fastqCleanup -- Makes sure sequence and quality strings are the same length, trims to the shorter of the two.\n"
      "usage:\n"
      "\tfastqCleanup [options] inread1.fq[.gz] inread2.fq[.gz] outread1.fq.gz outread2.fq.gz \n"
      "Options:\n"
      "\t-help\tPrints this message.\n"
      "\t-minLen=NUM\tIf at least one sequence is shorter than NUM after cleanup, don't print either (-1 to disable removal) (DEFAULT 30).\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"help",OPTION_STRING},
    {"minLen",OPTION_INT},
    {NULL, 0}
}; //end options()


/**
 * Main function, takes filenames for paired qseq reads
 * and outputs three files.
 */
int fastqCleanup(char *inread1, char *inread2, char *outread1, char *outread2, const int minlen){
  struct lineFile *lf1;
  struct lineFile *lf2;
  lf1 = lineFileOpen(inread1,TRUE);
  lf2 = lineFileOpen(inread2,TRUE);
  gzFile *zout1 = gzopen(outread1,"wb");
  gzFile *zout2 = gzopen(outread2,"wb");

  struct fastqItem *fq1 = allocFastqItem();
  struct fastqItem *fq2 = allocFastqItem();

  //loop over pairs of reads and clean up
  while(fastqItemNext(lf1,fq1) && fastqItemNext(lf2,fq2)){

    boolean trash1 = cleanAndThrowOutFastqItem(fq1,minlen);
    boolean trash2 = cleanAndThrowOutFastqItem(fq2,minlen);

    if((trash1 == TRUE) || (trash2 == TRUE))
      continue;

    gzPrintFastqItem(zout1,fq1);
    gzPrintFastqItem(zout2,fq2);
  }

  //close and clean up
  gzclose(zout1);
  gzclose(zout2);
  freeFastqItem(fq1);
  freeFastqItem(fq2);
  lineFileClose(&lf1);
  lineFileClose(&lf2);
  return(0);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  char *input1 = NULL;
  char *input2 = NULL;
  char *output1 = NULL;
  char *output2 = NULL;
  optionInit(&argc, argv, options);
  if(optionExists("help") || argc != 5){
    usage();
  }

  int minLen = optionInt("minLen", DEF_MIN_LEN);

  input1 = argv[1];
  input2 = argv[2];
  output1 = argv[3];
  output2 = argv[4];

  //make sure .gz is on the outfile name
  char newOut1[1000];
  char newOut2[1000];
  int outlen1 = strlen(output1);
  int outlen2 = strlen(output2);

  //put .gz at the end of output1 if needed
  if(tolower(output1[outlen1-2]) == 'g' && tolower(output1[outlen1-1]) == 'z'){
    strcpy(newOut1,output1);
  }else{
    sprintf(newOut1,"%s.gz",output1);
  }

  //put .gz at the end of output2 if needed
  if(tolower(output2[outlen2-2]) == 'g' && tolower(output2[outlen2-1]) == 'z'){
    strcpy(newOut2,output2);
  }else{
    sprintf(newOut2,"%s.gz",output2);
  }

  fastqCleanup(input1, input2, newOut1, newOut2, minLen);

  return(0);
} //end main()


