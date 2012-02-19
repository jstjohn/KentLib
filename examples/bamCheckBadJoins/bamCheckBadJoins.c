/**
 * bamCheckBadJoins
 *  Author: John St. John
 *  Date: 2/8/2012
 *
 *  This program looks for examples of weird mapping
 *  properties of a mate-pair file. Specifically it
 *  looks for cases where an inner scaffold mate strongly
 *  supports a join to a different genomic fragment.
 *
 */



/**
 * Includes
 *
 */

#include "common.h"
#include "options.h"
#include "linefile.h"
#include "bamFile.h"
#include "sam.h"



/**
 * Global Variables and definitions
 *
 */ 

#define DEFAULT_WINDOW 500
#define DEFAULT_EDGES 4000
#define DEFAULT_MQ 30



/**
 * Command line options
 */ 

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "bamCheckBadJoins -- check for support of different contig joins within scaffolds, reports summary.\n"
      "\tUsage: bamCheckBadJoins [options] bamFile.sorted.bam \n"
      "\noptions:\n"
      "\t-window\tthe window size for binning reads (default: 500)\n"
      "\t-edges\tskip this many bases at beginning and end of scaffolds (default: 4000)\n"
      "\t-minq\tonly consider alignments where both reads have at least mapq of (default: 30)\n"
      "\t-verbose\twrite some program status information to stderr.\n"
      "\t-help\twrite this help to the screen.\n"
      );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"window",OPTION_INT},
    {"edges",OPTION_INT},
    {"minq",OPTION_INT},
    {"help",OPTION_BOOLEAN},
    {"verbose",OPTION_BOOLEAN},
    {NULL, 0}
}; //end options()










/**
 * Program's functions
 *
 */




int bamPrintInfo(samfile_t *bamFile, FILE* out, int edges, int window, unsigned short minmq)
/* iterate through bam alignments, storing */
{
  int lastTID=-1; //real TIDs are never negative
  bam_header_t *header = bamFile->header;
  unsigned short *bucket_counts_to_other = NULL;
  unsigned short *bucket_counts = NULL;
  boolean skipTID = TRUE;
  unsigned int numBuckets = 0;
  unsigned length = 0;
  int i;

  bam1_t *b = bam_init1();
  while(bam_read1(bamFile->x.bam,b)){
    if(b->core.tid != lastTID){
      //we have an alignment to something
      //new
      if (bucket_counts_to_other != NULL && bucket_counts != NULL){
        char *name = header->target_name[lastTID];

        for(i=0;i<numBuckets;i++)
          fprintf(out, "%s\t%d\t%hu\t%hu\n", name, (window * i) + edges, bucket_counts_to_other[i], bucket_counts[i] );
        

        
        //free old buckets
        free(bucket_counts_to_other);
        bucket_counts_to_other = NULL;
        free(bucket_counts);
        bucket_counts = NULL;
      }

      lastTID = b->core.tid;
      
      length = header->target_len[lastTID];

      //nothing to see here
      if(length < 2 * edges){
        skipTID = TRUE;
        numBuckets = 0;
        continue;
      }
      
      //number of buckets to look at
      numBuckets = (length - (2 * edges) - window) / window;
      
      //calloc should 0 out arrays
      bucket_counts_to_other = (unsigned short *) calloc(numBuckets, sizeof(unsigned short));
      bucket_counts = (unsigned short *) calloc(numBuckets, sizeof(unsigned short));



    }
    if(skipTID == TRUE)
      continue;


    //now increment counts of our buckets
    if(b->core.qual >= minmq 
        && !(b->core.flag & BAM_FUNMAP) 
        && !(b->core.flag & BAM_FMUNMAP)
        && b->core.pos >= edges
        && b->core.pos <= length - edges)
    {
      int bucket_idx = (b->core.pos - edges) / window;
      bucket_counts[bucket_idx] += 1;

      if(b->core.mtid != b->core.tid){
        //the mate is mapped to a different sequece
        bucket_counts_to_other[bucket_idx] += 1;
      }

    }

    
  }
  bam_destroy1(b);
}







int main(int argc, char *argv[])
/* Process command line. */
{

  optionInit(&argc, argv, options);
  boolean help = optionExists("help");
  boolean verbose = optionExists("verbose");
  if(help || argc != 2) usage();
  if(!bamFileExists(argv[1])){
    fprintf(stderr,"ERROR: Bam file %s could not be found, or bam index not found.\n",argv[1]);
    usage();
  }

  int windowsize = optionInt("window",DEFAULT_WINDOW);
  int edgelen = optionInt("edges",DEFAULT_EDGES);
  int minq = optionInt("minq",DEFAULT_MQ);

  char *bamFileName;
  samfile_t *bamFile = bamOpen(argv[1],&bamFileName);
  
  bamPrintInfo(bamFile, stdout, edgelen, windowsize, minq);

  bamClose(&bamFile);  


  return 0;
} //end main()


