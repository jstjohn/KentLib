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



/*
Global variables
 */
bool verboseOut = false;

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "qseqToFastq -- Converts illumina qseq format to fastq format, outputs a pair 1, 2 and singleton file for the reads where only one pair passes the Illumina QC filter\n"
      "usage:\n"
      "\tqseqToFastq [required options] [options] \n"
      "\n**required** options:\n"
      "\t-read1=FILE\tFile name holding the first qseq mate\n"
      "\t-read2=FILE\tFile name holding the second qseq mate.\n"
      "\t-prefix=NAME\tPrefix for output, example 's_5_', would result in s_5_1.fastq.gz (read 1) s_5_2.fastq.gz (read 2) and s_5_s.fastq.gz (singleton reads)\n"
      "\noptions:\n"
      "\t-verbose\tOutput verbose debug messages to stderr.\n\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"read1",OPTION_STRING},
    {"read2",OPTION_STRING},
    {"prefix",OPTION_STRING},
    {"verbose",OPTION_BOOLEAN},
    {NULL, 0}
}; //end options()



int printFastqFromSplitQseq(gzFile *zfp, char *qseq[11]){
  int i;
  for(i=0;i<strlen(qseq[8]);i++){
    if(qseq[8][i] == '.') qseq[8][i] = 'N';
  }
  gzprintf(zfp,"@%s_%s:%s:%s:%s:%s#%s/%s\n%s\n+\n%s\n",
      qseq[0],qseq[1], // machine ID_Run number:
      qseq[2],qseq[3], // lane number: tile number
      qseq[4],qseq[5], //x coordinate: y coordinate
      qseq[6],         // multiplex index (0 for none)
      qseq[7],         // read number (1 or 2)
      qseq[8],         // sequence
      qseq[9]          // quality string
      );

}


/**
 * Main function, takes filenames for paired qseq reads
 * and outputs three files.
 */
int qseqToFastq(char *inread1, char *inread2, char *outread1, char *outread2, char *outreads){
  struct lineFile *lf1;
  struct lineFile *lf2;
  lf1 = lineFileOpen(inread1,TRUE);
  lf2 = lineFileOpen(inread2,TRUE);
  gzFile *zout1 = gzopen(outread1,"wb");
  gzFile *zout2 = gzopen(outread2,"wb");
  gzFile *zouts = gzopen(outreads,"wb");
  char *line1;
  char *line2;
  while(lineFileNextReal(lf1,&line1) && lineFileNextReal(lf2,&line2)){
    char *split1[11];
    char *split2[11];
    chopByWhite(line1,split1,11);
    chopByWhite(line2,split2,11);
    int qc1 = atoi(split1[10]);
    int qc2 = atoi(split2[10]);
    if(qc1 == 1 && qc2 == 1){
      printFastqFromSplitQseq(zout1, split1);
      printFastqFromSplitQseq(zout2, split2);
    }else if(qc1 == 1){
      printFastqFromSplitQseq(zouts, split1);
    }else if(qc2 == 1){
      printFastqFromSplitQseq(zouts, split2);
    }
  }
  gzclose(zout1);
  gzclose(zout2);
  gzclose(zouts);
  return 0;
}


int main(int argc, char *argv[])
/* Process command line. */
{
  char *read1 = NULL;
  char *read2 = NULL;
  char *prefix= NULL;
  optionInit(&argc, argv, options);
  verboseOut = optionExists("verbose");
  read1 = optionVal("read1",NULL);
  read2 = optionVal("read2",NULL);
  prefix = optionVal("prefix",NULL);

  if (read1 == NULL || read2 == NULL || prefix == NULL){
    if(read1 == NULL) fprintf(stderr,"Must provide -read1\n");
    if(read2 == NULL) fprintf(stderr,"Must provide -read2\n");
    if(prefix == NULL) fprintf(stderr,"Must provide -prefix\n");
    usage();
  }
  char out1[255];
  char out2[255];
  char outs[255];
  sprintf(out1,"%s1.fastq.gz",prefix);
  sprintf(out2,"%s2.fastq.gz",prefix);
  sprintf(outs,"%ss.fastq.gz",prefix);
  qseqToFastq(read1,read2,out1,out2,outs);

  return 0;
} //end main()


