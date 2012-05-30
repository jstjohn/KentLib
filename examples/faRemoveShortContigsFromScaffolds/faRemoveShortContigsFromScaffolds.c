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
#include "dnaLoad.h"
#include <math.h>
#include <ctype.h>



/**
 * Global Variables and definitions
 *
 */

#define MIN_SEQ_LEN 200
#define MIN_GAP_LEN 10
#define MIN_TO_CALL_GAP_LEN 2
#define LINE_WRAP_LEN 50




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
      "\t-lineWrapLen=NUM\tLength of line wrap in output fasta file. (default: 50)\n"
      "\t-minToCallGapLen=NUM\tMinimum length of a stretch of 'N' to be considered a gap. (default: 2)\n"
      "\t-minGapLen=NUM\tIf a gap is found, minimum length of that gap in output file. (default: 10)\n"
      "\t-help\tWrites this help to the screen, and exits.\n"
      );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"help",OPTION_BOOLEAN},
    {"minLen",OPTION_INT},
    {"lineWrapLen",OPTION_INT},
    {"minToCallGapLen",OPTION_INT},
    {"minGapLen",OPTION_INT},
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
      return(FALSE);
    }
  }
  return(TRUE);
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

inline void faWriteCWithMinGaps(FILE *f, char c, const int maxPerLine, int *leftPerLine){
  fputc(c,f);
  if((--(*leftPerLine)) == 0){
    fputc('\n',f);
    (*leftPerLine) = maxPerLine;
  }
}

void faWriteWithMinGaps(FILE *f, char *startLine,
    const DNA *const letters, const int dnaSize,
    const int maxPerLine,
    const int minToCallGapLen,
    const int minGapLen){
  int lettersLeft = dnaSize;
  int leftPerLine = maxPerLine;
  int gapLen = 0;
  if (startLine != NULL)
      fprintf(f, ">%s\n", startLine);
  int i;
  for(i=0;i<dnaSize;i++){
    char c = letters[i];
    char uc = toupper(c);

    if(uc == 'N'){
      gapLen++;
      faWriteCWithMinGaps(f,'N',maxPerLine,&leftPerLine);
    }else{

      //deal with previous gap
      if(gapLen > 0){
        if(gapLen >= minToCallGapLen && gapLen < minGapLen){
          int j;
          for(j=gapLen; j<minGapLen;j++)
            faWriteCWithMinGaps(f,'N',maxPerLine,&leftPerLine);
        }
        //reset gapLen
        gapLen = 0;
      }
      //now write out the char
      faWriteCWithMinGaps(f,c,maxPerLine,&leftPerLine);
    }
  }
  //newline after sequence
  fputc('\n',f);
}

void faRemoveShortContigsFromScaffolds(char *faFile, FILE *outstream,
    const int minLen, const int minToCallGapLen, const int minGapLen, const int lineWrapLen){

  struct lineFile *lf = lineFileOpen(faFile,TRUE);
  DNA *seq;
  int seqLen , i = 0;
  char *seqName;


  while(faMixedSpeedReadNext(lf, &seq, &seqLen, &seqName)){
    int ctgLen = 0;
    int gapLen = 0;
    for(i=0;i<seqLen;i++){
      if(toupper(seq[i]) == 'N'){
        gapLen++;
        if(ctgLen-gapLen+1 > 0 && ctgLen-gapLen+1 < minLen && gapLen == minToCallGapLen){
          //need to mask
          //we have seen minGapLen N's already
          //and ctgLen has been incremented by minGapLen-1 N's
          //into this gap
          maskRange(seq,i-ctgLen,i-gapLen+1);
        }
        if(gapLen >= minToCallGapLen){
          ctgLen = 0; //make sure it is 0 since we have seen a real gap
        }else{
          ctgLen++;
        }
      }else{
        ctgLen++;
        gapLen = 0;
      }
    }
    if(ctgLen-gapLen+1 > 0 && ctgLen-gapLen+1 < minLen){
      //need to mask
      //we have seen minGapLen N's already
      //and ctgLen has been incremented by minGapLen-1 N's
      //into this gap
      maskRange(seq,i-ctgLen,i-gapLen+1);
    }

    int beginPos = beginNTrimPos(seq,seqLen);
    int endPos = endNTrimPos(seq,seqLen);


    //if the sequence is longer than our minLen
    //and not all gaps, then print it!
    if(endPos <= beginPos)
      continue;

    int afterTrimLen = endPos - beginPos;
    if(afterTrimLen < minLen)
      continue;

    if(allN(seq+beginPos, afterTrimLen) == FALSE){
      //faWriteNext(outstream, seqName, seq+beginPos, afterTrimLen);
      faWriteWithMinGaps(outstream, seqName, seq+beginPos, afterTrimLen, lineWrapLen, minToCallGapLen, minGapLen);
    }
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
  int minGapLen = optionInt("minGapLen",MIN_GAP_LEN);
  int minToCallGapLen = optionInt("minToCallGapLen",MIN_TO_CALL_GAP_LEN);
  int lineWrapLen = optionInt("lineWrapLen",LINE_WRAP_LEN);

  if(help) usage();
  if(argc != 3) usage();
  FILE *outstream = fopen(argv[2],"w");
  faRemoveShortContigsFromScaffolds(argv[1], outstream, minLen, minToCallGapLen, minGapLen, lineWrapLen);
  fclose(outstream);
  return(0);
} //end main()


