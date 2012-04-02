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
bool verboseOut = false;

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "fastq64to33 -- Converts fastq reads of a file in phred+64 format to phred+33 format.\n"
      "usage:\n"
      "\tfastq64to33 [required options] \n"
      "\n**required** options:\n"
      "\t-input=FILE\tFile name holding the fastq reads to be converted\n"
      "\t-output=NAME\tFilename to output\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"input",OPTION_STRING},
    {"output",OPTION_STRING},
    {NULL, 0}
}; //end options()

/**
 * Main function, takes filenames for paired qseq reads
 * and outputs three files.
 */
int fastq64To33(char *inread, char *outread){
  struct lineFile *lf;
  lf = lineFileOpen(inread,TRUE);
  gzFile *zout = gzopen(outread,"wb");
  struct fastqItem *fq = allocFastqItem();
  while(fastqItemNext(lf,fq)){
    convPhred64ToPhred33(fq);
    gzPrintFastqItem(zout,fq);
  }
  gzclose(zout);
  freeFastqItem(fq);
  return(0);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  char *input = NULL;
  char *output = NULL;
  optionInit(&argc, argv, options);

  input = optionVal("input",NULL);
  output = optionVal("output",NULL);

  if (input == NULL || output == NULL){
    if(input == NULL) fprintf(stderr,"Must provide -input=FILE\n");
    if(output == NULL) fprintf(stderr,"Must provide -output=FILE\n");
    usage();
  }


  //make sure .gz is on the outfile name
  char newOut[1000];
  int outlen = strlen(output);
  if(tolower(output[outlen-2]) == 'g' && tolower(output[outlen-1]) == 'z'){
    strcpy(newOut,output);
  }else{
    sprintf(newOut,"%s.gz",output);
  }

  fastq64To33(input, newOut);

  return(0);
} //end main()


