/**
 * helloWorld
 *  Author: John St. John
 *  Date: 2/8/2012
 *
 *  This program is an example shell of a C program
 *  for working with the kent source library.
 *
 */



/**
 * Includes
 *
 */

#include "common.h"
#include "options.h"
#include "linefile.h"
#include "fa.h"
#include <math.h>
#include <ctype.h>



/**
 * Global Variables and definitions
 *
 */

#define MIN_SEQ_LEN 200




/**
 * Command line options
 */ 

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "faRemoveShortContigsFromScaffolds -- Given a fasta file, and the filtered set of bad looking positions from bamCheckBadJoins2, split the sequence on any 'N' characters that appear in the marked bad regions and output the resulting fragmented assembly.\n"
      "\tUsage: faRemoveShortContigsFromScaffolds [options] genome.scaf.fa out.genome.scaf.fa\n"
      "\noptions:\n"
      "\t-minLen=NUM\tMinimum contig length in scaffold file to keep. (default: 200)\n"
      "\t-help\tWrites this help to the screen, and exits.\n"
      );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"help",OPTION_BOOLEAN},
    {"minLen",OPTION_INT},
    {NULL, 0}
}; //end options()










/**
 * Program's functions
 *
 */

inline void maskRange(DNA *seq, int start, int end){
  int i;
  for(i=start;i<end;i++){
    seq[i] = 'N';
  }
}

inline boolean allN(DNA *seq, int len){
  int i;
  for(i=0;i<len;i++){
    if(toupper(seq[i]) != 'N'){
      return(TRUE);
    }
  }
  return(FALSE);
}

inline int beginNTrimPos(DNA *seq, int len){
  int i = 0;
  for(i=0;i<len;i++){
    if(toupper(seq[i]) != 'N')
      return(i);
  }
  return(i);
}

inline int endNTrimPos(DNA *seq, int len){
  int i = 0;
  for(i=len-1; i >= 0; i--){
    if(toupper(seq[i]) != 'N')
      return(i+1);
  }
  return(i+1);
}

void faRemoveShortContigsFromScaffolds(char *faFile, FILE *outstream, int minLen){
  struct lineFile *lf = lineFileOpen(faFile,TRUE);
  DNA *seq;
  int seqLen , i = 0;
  char *seqName;

  while(faMixedSpeedReadNext(lf, &seq, &seqLen, &seqName)){

    int ctgLen = 0;
    for(i=0;i<seqLen;i++){
      if(toupper(seq[i]) == 'N'){
        if(ctgLen > 0 && ctgLen < minLen){
          //need to mask
          maskRange(seq,i-ctgLen,i);
        }
        ctgLen = 0; //make sure it is 0 since we have seen a gap
      }else{
        ctgLen++;
      }
    }

    int beginPos = beginNTrimPos(seq,seqLen);
    int endPos = endNTrimPos(seq,seqLen);


    //if the sequence is longer than our minLen
    //and not all gaps, then print it!
    if(endPos <= beginPos)
      continue; //discard sequence

    int afterTrimLen = endPos - beginPos;
    if(afterTrimLen < minLen)
      continue; //discard if too short

    if(allN(seq+beginPos, afterTrimLen) == FALSE)
      faWriteNext(outstream, seqName, seq+beginPos, afterTrimLen);

  }//end loop over fasta sequences
  faFreeFastBuf();
  lineFileClose(&lf);
}



int main(int argc, char *argv[])
/* Process command line. */
{

  optionInit(&argc, argv, options);
  boolean help = optionExists("help");
  int minLen = optionInt("minLen",MIN_SEQ_LEN);
  if(help) usage();
  if(argc != 3) usage();
  FILE *outstream = fopen(argv[2],"w");
  faRemoveShortContigsFromScaffolds(argv[1], outstream, minLen);
  fclose(outstream);
  return(0);
} //end main()


