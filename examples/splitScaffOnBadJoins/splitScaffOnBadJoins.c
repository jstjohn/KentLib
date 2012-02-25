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

#include "uthash.h"




/**
 * Global Variables and definitions
 *
 */

#define MAX_NAME_LEN 50
#define MAX_SAME_RANGE 1000
#define MIN_SEQ_LEN 100

struct begin_end_lst{
  int begin;
  int end;
  struct begin_end_lst *next;
};

struct chrom_to_begin_end_lst{
  char name[MAX_NAME_LEN];
  struct begin_end_lst *head;
  struct begin_end_lst *tail;
  UT_hash_handle hh;
};


struct chrom_to_begin_end_lst *chromToBeginEndLstDict = NULL;





/**
 * Command line options
 */ 

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "splitScaffOnBadJoins -- Given a fasta file, and the filtered set of bad looking positions from bamCheckBadJoins2, split the sequence on any 'N' characters that appear in the marked bad regions and output the resulting fragmented assembly.\n"
      "\tUsage: splitScaffOnBadJoins [options] genome.fa badJoins.filtered.tab \n"
      "\noptions:\n"
      "\t-help\tWrites this help to the screen, and exits.\n"
      );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"help",OPTION_BOOLEAN},
    {NULL, 0}
}; //end options()










/**
 * Program's functions
 *
 */

void addAndOrCreateToBegEndLstDict(char *chrom, int begin, int end){
  struct chrom_to_begin_end_lst *s;
  HASH_FIND_STR(chromToBeginEndLstDict, chrom, s);
  if(s == NULL){
    s = (struct chrom_to_begin_end_lst *) malloc(sizeof(struct chrom_to_begin_end_lst));
    strcpy(s->name, chrom);
    s->head = (struct begin_end_lst *)malloc(sizeof(struct begin_end_lst));
    s->tail = s->head; //only one element
    s->tail->next = NULL;
    HASH_ADD_STR(chromToBeginEndLstDict, name, s);//add this item to the hashtable
  }else{
    s->tail->next = (struct begin_end_lst *)malloc(sizeof(struct begin_end_lst));
    s->tail = s->tail->next;
    s->tail->next = NULL;
  }
  s->tail->begin = begin;
  s->tail->end = end;
}


void parseBegEndLstFile(char *chromBegEndLstFile){
  struct lineFile *lf = lineFileOpen(chromBegEndLstFile,TRUE);
  char lastChrom[MAX_NAME_LEN]; 
  strcpy(lastChrom,"NoTAcHrOmNaMe__223");
  char *line;
  int lastPos = -1;
  int lastStart = -1;


  while(lineFileNextReal(lf,&line)){
    char *split[3];
    chopByWhite(line,split,3);
    int pos = atoi(split[1]);
    char *chrom = split[0];
    if(strcmp(chrom,lastChrom) != 0){
      if(lastPos != -1 && lastStart != -1){
        //need to add this to the dict
        addAndOrCreateToBegEndLstDict(lastChrom, lastStart, lastPos +1);  
      }
      //set the last chrom to this chrom
      strcpy(lastChrom, chrom);
      lastStart = pos;
    }else{
      //same sid
      if(abs(pos - lastPos) > MAX_SAME_RANGE){
        addAndOrCreateToBegEndLstDict(lastChrom, lastStart, lastPos +1);
        lastStart = pos;
      }
    }
    
    lastPos = pos;

  }//end loop over lines
  //take care of last range, if it exists
  if(lastPos != -1 && lastStart != -1){
    addAndOrCreateToBegEndLstDict(lastChrom, lastStart, lastPos +1);
  }
}


void splitFaOnNsInBeginEndRegions(char *faFile){
  struct lineFile *lf = lineFileOpen(faFile,TRUE);
  DNA *seq;
  int seqLen , i = 0;
  char *seqName;
  struct chrom_to_begin_end_lst *s;
  char tmpName[MAX_NAME_LEN+10];

  while(faSpeedReadNext(lf, &seq, &seqLen, &seqName)){

    HASH_FIND_STR(chromToBeginEndLstDict, seqName, s);

    //case 1: there aren't any splits
    if(s == NULL)
      faWriteNext(stdout, seqName, seq, seqLen);
    //case 2: deal with splits
    else{
      int numSplits = 0;
      struct begin_end_lst *bel;

      //replace all N regions overlapping start->end regions of bad bases
      for(bel = s->head; bel != NULL; bel = bel->next){
        //get rid of N's of region starting before bel->begin but going to or past it
        for(i=bel->begin; i>=0; i--){
          if(seq[i] != '\0' && toupper(seq[i]) == 'N'){
            seq[i] = '\0'; //swap all N's in the bad region with null terminators
          }else{
            break;
          }
        }
        
        //now swap out termination char for all N's within the start->end block
        for(i = bel->begin; i < bel->end; i++){
          if(seq[i] == '\0')
            continue;
          else if(toupper(seq[i]) == 'N'){
            seq[i] = '\0';
          }
        }

        if(seq[bel->end-1] == '\0'){
          //we need to look past bel->end until we see a non-N
          for(i=bel->end;i<seqLen;i++){
            if(seq[i] != '\0' && toupper(seq[i]) == 'N')
              seq[i] = '\0';
          }
        }
      }//end loop over bad regions

      //Now all of the N's in and around bad regions are '\0'
      //first skip any initial Ns
      int count = 0;
      int start;
      for(start = 0; start < seqLen; start++){
        if((! seq[i] == '\0') && toupper(seq[i]) != 'N' ){
          int length = strlen(seq+start); //finds length until next '\0' char
          if(length >= MIN_SEQ_LEN){
            sprintf(tmpName,"%s_%d",seqName,++count);
            faWriteNext(stdout, tmpName, seq+start, length);
          }
          start += length-1; //incrememted by one by the loop
        }
      }//end loop over seq for printing
    }//end deal with breaks
  }//end loop over fasta sequences

  faFreeFastBuf();
}



int main(int argc, char *argv[])
/* Process command line. */
{

  optionInit(&argc, argv, options);
  boolean help = optionExists("help");
  if(help) usage();
  if(argc != 3) usage();
  
  parseBegEndLstFile(argv[2]);
  splitFaOnNsInBeginEndRegions(argv[1]);
  //WARNING, if we do anything else the above has memory leaks. OK since program is
  //now done though.
  
  return 0;
} //end main()


