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
#define MAX_OCOUNT_LEN 301



/**
 * Command line options
 */ 

void usage()
  /* Explain usage and exit. */
{
  errAbort(
      "bamCheckBadJoins -- WARNING: CURRENTLY HAS BUGS, USE V2: check for support of different contig joins within scaffolds, reports summary.\n"
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

int compare_unsigned (const unsigned *a, const unsigned *b)
{
  unsigned temp = *a - *b;
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

void setAdd(unsigned *lst, unsigned mtid){
  if(lst[0] > 1){
    //search after the first element of the list for mtid
    unsigned *rlst = bsearch(&mtid, lst+1, lst[0], sizeof(unsigned), compare_unsigned );
    if(rlst == NULL && lst[0] < MAX_OCOUNT_LEN){
      //insert the new element into the list at the current end of the list
      lst[lst[0]+1] = mtid;
      lst[0]++;
      //now re-sort the list after the first element
      qsort(lst+1, lst[0], sizeof(unsigned), compare_unsigned );

    }
  }else{ //base case
    if (lst[0] == 0){
      lst[1] = mtid;
      lst[0] = 1;
    }else{ //exactly one element
      if(mtid != lst[1]){
        if(mtid < lst[1]){
          unsigned tmp = lst[1];
          lst[1] = mtid;
          lst[2] = tmp;
        }else{
          lst[2] = mtid;
        }
        lst[0]++;
      }
    }

  }
}


unsigned ** allocSetLst(int length){
  int i = 0;
  unsigned ** setList = (unsigned **) malloc(length * sizeof(unsigned *));
  for(i=0;i<length;i++){
    setList[i] = (unsigned *) calloc(MAX_OCOUNT_LEN, sizeof(unsigned));
  }
}

void freeSetLst(unsigned **setLst, int length){
  int i = 0;
  for(i=0;i<length;i++){
    free(setLst[i]);
  }
  free(setLst);
}


int bamPrintInfo(samfile_t *bamFile, FILE* out, int edges, int window, unsigned short minmq)
  /* iterate through bam alignments, storing */
{
  int lastTID=-1; //real TIDs are never negative
  bam_header_t *header = bamFile->header;
  unsigned short *bucket_counts_to_other = NULL;
  unsigned short *bucket_counts = NULL;
  unsigned **bucket_num_other = NULL;
  boolean skipTID = TRUE;
  int numBuckets = 0;
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
          fprintf(out, "%s\t%d\t%hu\t%u\t%hu\n", name, (window * i) + edges, bucket_counts_to_other[i], bucket_num_other[i][0], bucket_counts[i] );



        //free old buckets
        free(bucket_counts_to_other);
        bucket_counts_to_other = NULL;
        free(bucket_counts);
        bucket_counts = NULL;
        freeSetLst(bucket_num_other, numBuckets);
        bucket_num_other = NULL;
      }

      lastTID = b->core.tid;

      length = header->target_len[lastTID];

      //nothing to see here
      if(length < (2 * edges + window)){
        skipTID = TRUE;
        numBuckets = 0;
        continue;
      }

      skipTID = FALSE;

      //number of buckets to look at
      numBuckets = (length - (2 * edges) - window) / window;

      //calloc should 0 out arrays
      bucket_counts_to_other = (unsigned short *) calloc(numBuckets, sizeof(unsigned short));
      bucket_counts = (unsigned short *) calloc(numBuckets, sizeof(unsigned short));
      bucket_num_other = allocSetLst(numBuckets);
   
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
      
      if(bucket_idx >= numBuckets)
        continue;

      bucket_counts[bucket_idx]++;

      if(b->core.mtid != b->core.tid){
        //the mate is mapped to a different sequece
        bucket_counts_to_other[bucket_idx]++;

        //add this mtid to the set of mtids for this bucket
        setAdd(bucket_num_other[bucket_idx], b->core.mtid);

      }

    }//end bucket incrementation


  }//end loop over reads

  
  //now that loop is done, take care of last bucket if it is there.
  if (bucket_counts_to_other != NULL && bucket_counts != NULL){
    char *name = header->target_name[lastTID];

    for(i=0;i<numBuckets;i++)
     fprintf(out, "%s\t%d\t%hu\t%u\t%hu\n", name, (window * i) + edges, bucket_counts_to_other[i], bucket_num_other[i][0], bucket_counts[i] );


    //free old buckets
    free(bucket_counts_to_other);
    bucket_counts_to_other = NULL;
    free(bucket_counts);
    bucket_counts = NULL;
    freeSetLst(bucket_num_other, numBuckets);
    bucket_num_other = NULL;

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


