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
#include "sam.h"



/**
 * Global Variables and definitions
 *
 */ 

#define DEFAULT_MQ 30
#define MAX_OCOUNT_LEN 301
#define DEFAULT_MIN_INSERT 10
#define DEFAULT_MAX_INSERT 5000
#define MAX_COUNT (1000)


/**
 * Command line options
 */ 

void usage()
  /* Explain usage and exit. */
{
  errAbort(
      "bamFragmentCoverage -- make list of per-base regions of unexpected mate-pair insert sizes.\n"
      "Usage: bamFragmentCoverage [options] bamFile.sorted.bam \n"
      "\noptions:\n"
      "\t-minInsert=INT\tread pairs with inserts less than this value aren't considered (default: 10)\n"
      "\t-maxInsert=INT\tread pairs with inserts greater than this value are marked as bad (default: 5000)\n"
      "\t-minq=INT\tonly consider alignments where the left read has at least a mapq of (default: 30)\n"
      "\t-verbose\twrite some program status information to stderr.\n"
      "\t-help\twrite this help to the screen.\n"
      );
}//end usage()


static struct optionSpec options[] = {
  /* Structure holding command line options */
  {"minInsert",OPTION_INT},
  {"maxInsert",OPTION_INT},
  {"minq",OPTION_INT},
  {"help",OPTION_BOOLEAN},
  {"verbose",OPTION_BOOLEAN},
  {NULL, 0}
}; //end options()



/**
 * Program's functions
 *
 */

int addReadCovToCovLst(bam1_t *b, unsigned short *insert_coverage_counts){
  int i,k,l,op,end;
  bam1_core_t *c = &(b->core);
  unsigned int tpos = c->pos;
  unsigned int tpos_ori = tpos;
  k = 0;
  unsigned int *cigar;
  cigar = bam1_cigar(b);
  i = 0;


  for(k=0;k<c->n_cigar;k++){
    op = cigar[k] & BAM_CIGAR_MASK;
    l = cigar[k] >> BAM_CIGAR_SHIFT;

    switch(op)
    {
    case BAM_CSOFT_CLIP:
      i+=l;
      break;

    case BAM_CHARD_CLIP:
      //do nothing, the sequence was physically removed
      //i need not be incremented
      break;

    case BAM_CINS: //insertion relative to the reference
      i += l;
      break;

    case BAM_CDEL: //deletion relative to the reference
      end = tpos+l;
      for(;tpos<end;tpos++){
        if(insert_coverage_counts[tpos] < MAX_COUNT)
          insert_coverage_counts[tpos]++;
      }
      break;

    case BAM_CREF_SKIP: //skip stuff in the reference
      tpos += l;
      break;

    case BAM_CPAD: //silent deletion from padded reference
      i += l;
      break;


      //treat match, equal or diff the same, we still need to know the actual
      //underlying bases
    case BAM_CMATCH:
    case BAM_CEQUAL:
    case BAM_CDIFF:
      end = i+l;
      for(;i<end;i++,tpos++){
        if(insert_coverage_counts[tpos] < MAX_COUNT)
           insert_coverage_counts[tpos]++;
      }
      break;

    default:
      break;
    }//done switching over alignment

  }//done loop through cigar alignment things
  return(tpos - tpos_ori); //return the length on the genome that this read covers
}




void bamPrintInfo(samfile_t *bamFile, FILE* out, int minInsert, int maxInsert, int minmq, boolean verbose)
  /* iterate through bam alignments, storing */
{
  int lastTID=-1; //real TIDs are never negative
  bam_header_t *header = bamFile->header;
  unsigned short *insert_coverage_counts = NULL;
  boolean skipTID = TRUE;
  int length = 0;
  int i;

  bam1_t *b = bam_init1();
  fprintf(out, "#seq_name\tposition(0-based)\tisize_out_of_range\tdiscontiguous_in_avg_insert_window\tok_looking_inserts\n");
  while(samread(bamFile, b)>=0)
  {
    if(b->core.tid != lastTID){
      //we have an alignment to something
      //new
      if (insert_coverage_counts != NULL){
        char *name = header->target_name[lastTID];

        for(i=0;i<length;i++)
          fprintf(out, "%s\t%d\t%hu\n", name, i, insert_coverage_counts[i] );

        //free old count structures
        free(insert_coverage_counts);
        insert_coverage_counts = NULL;

      }

      lastTID = b->core.tid;

      length = header->target_len[lastTID];

      skipTID = FALSE;

      //calloc should 0 out arrays

      insert_coverage_counts = (unsigned short *) calloc(length, sizeof(unsigned short));
   
    }
    if(skipTID == TRUE)
      continue;


    //now increment counts of our 
    if(b->core.qual >= minmq && (0 == (b->core.flag & BAM_FUNMAP)))
    { //this read aligns somewhere

      int seqlen = addReadCovToCovLst(b,insert_coverage_counts);

      if((!(b->core.flag & BAM_FUNMAP)) && (b->core.tid == b->core.mtid)){
        //read and mate both map
        int mypos = b->core.pos;
        int opos = b->core.mpos;
        int endpos = mypos+seqlen;
        int gaplen = opos-endpos;

        if( gaplen > 0 && gaplen >= minInsert && gaplen <= maxInsert){
          for(i=endpos;i<opos;i++){
            if(insert_coverage_counts[i] < MAX_COUNT)
              insert_coverage_counts[i]++;
          }
        }
      }
    }//end coverage incrementation


  }//end loop over reads

  
  //now that loop is done, take care of last bucket if it is there.
  if (insert_coverage_counts != NULL){
    char *name = header->target_name[lastTID];

    for(i=0;i<length;i++)
      fprintf(out, "%s\t%d\t%hu\n", name, i, insert_coverage_counts[i] );



    //free old count structures
    free(insert_coverage_counts);
    insert_coverage_counts = NULL;

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

  int maxInsert = optionInt("maxInsert",DEFAULT_MAX_INSERT);
  int minInsert = optionInt("minInsert",DEFAULT_MIN_INSERT);
  int minq = optionInt("minq",DEFAULT_MQ);

  samfile_t *bamFile = samopen(argv[1],"rb",NULL);
  if(!bamFile){
      fprintf(stderr,"ERROR: Bam file %s could not be opened.\n",argv[1]);
      usage();
    }


  bamPrintInfo(bamFile, stdout, minInsert, maxInsert, minq, verbose);

  samclose(bamFile);


  return(0);
} //end main()


