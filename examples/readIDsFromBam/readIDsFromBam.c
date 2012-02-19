/*
 * readIDsFromBam.c
 *
 *  Created on: Feb 9, 2012
 *      Author: jstjohn
 */


/**
 * Includes
 *
 */

#include "common.h"
#include "options.h"
#include "linefile.h"
#include "uthash.h"
#include "sam.h"


/**
 * Global Variables and definitions
 *
 */
struct nameSet{
  char name[50];
  UT_hash_handle hh;
};


struct nameSet *badSequences = NULL;





/**
 * Command line options
 */

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "readIDsFromBam -- Given a sorted and indexed BAM alignment, and a blacklist (one per line) of sequence names, prints the IDs of all of the reads that map to these sequences.\n"
      "\tUsage: readIDsFromBam reads.sorted.bam blacklist.txt \n"
      );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {NULL, 0}
}; //end options()










/**
 * Program's functions
 *
 */


int initBlackListSet(char* infile){
  /**
   * Read in the blacklist file, and store it in the
   * global set of badSequences
   */
  struct lineFile *lf;
  lf = lineFileOpen(infile,TRUE);
  char *line;
  struct nameSet *s;
  while(lineFileNextReal(lf,&line)){
    char *seqid = trimSpaces(line);
    if (strlen(seqid) < 1) continue; //skip blank lines
    s = malloc(sizeof(struct nameSet));
    strcpy(s->name, seqid);
    HASH_ADD_STR(badSequences,name,s);
  }
  return 0;
}


int getFastqIDsFromBam(char *bamfile){
 /**
  * Grab the sequences aligning to each of our blacklist regions,
  * and print out their sequence IDs.
  */
  return 0;
}



int main(int argc, char *argv[])
/* Process command line. */
{

  optionInit(&argc, argv, options);
  boolean help = optionExists("help");
  if(help) usage();
  if(argc != 3) usage();


  initBlackListSet(argv[1]);

  getFastqIDsFromBam(argv[2]);


  return 0;
} //end main()


