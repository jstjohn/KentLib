/*
	Author John St. John
	07/12/2010

	Define a few structures and functions for
	reading and writing fastq format.

	Also defines a few functions for converting
	between various phred and ascii phred
	format.
*/

#include <ctype.h>
#include "common.h"
#include "linefile.h"
#include "dnautil.h"
#include <math.h>
#include <string.h>
#include "fastq.h"


/* Function Prototypes */
/*
void convPhred33ToPhred64( struct fastqItem * );
void convPhred64ToPhred33( struct fastqItem * );
void phred33ToPhred64( char *, int );
void phred64ToPhred33( char *, int );
int  phred64ToPhred( char );
int  phred33ToPhred( char );
double phredToDouble( int );
double phred33ToDouble( char );
double phred64ToDouble( char );
int doubleToPhred( double );
char phredToPhred33( int );
char phredToPhred64( int );
boolean fastqItemNext(struct lineFile *, struct fastqItem *);
struct fastqItem * allocFastqItem();
void freeFastqItem(struct fastqItem * fq);
void printFastqItem(FILE *fp, struct fastqItem *fq);
*/



struct fastqItem * allocFastqItem()
{
  struct fastqItem *fq = (struct fastqItem *) malloc(sizeof(struct fastqItem ));
  fq->id = (char *)malloc(sizeof(char) * MAX_ID_LENGTH);
  fq->seq = (char *)malloc(sizeof(char) * MAX_SEQ_LENGTH);
  fq->score = (char*)malloc(sizeof(char) * MAX_SEQ_LENGTH);
  return fq;
}

void freeFastqItem(struct fastqItem * fq)
{
  if(fq)
  {
    if(fq->id)
    {
      free(fq->id);
      fq->id = NULL;
    }
    if(fq->seq)
    {
      free(fq->seq);
      fq->seq = NULL;
    }
    if(fq->score)
    {
      free(fq->score);
      fq->score = NULL;
    }
    free(fq);
    fq = NULL;
  }
}


inline void printFastqItem(FILE *fp, struct fastqItem *fq)
{
  int i;
  fprintf(fp,"@%s\n",fq->id);
  for(i=0; i < fq->len; i++)
  {
    fprintf(fp,"%c",fq->seq[i]);
  }
  fprintf(fp,"\n+\n");
  for(i=0; i< fq->len; i++)
  {
    fprintf(fp,"%c",fq->score[i]);
  }
  fprintf(fp,"\n");
}

inline void gzPrintFastqItem(gzFile *fp, struct fastqItem *fq){
  int i;

  gzprintf(fp,"@%s\n",fq->id);
  for(i=0; i < fq->len; i++)
  {
    gzprintf(fp,"%c",fq->seq[i]);
  }
  gzprintf(fp,"\n+\n");
  for(i=0; i< fq->len; i++)
  {
    gzprintf(fp,"%c",fq->score[i]);
  }
  gzprintf(fp,"\n");
}


inline boolean cleanAndThrowOutFastqItem(struct fastqItem *fq, const int minlen){
  boolean throwOut = FALSE;
  int seql = strlen(fq->seq);
  int quall = strlen(fq->score);

  if(seql != quall){
    int minl = min(seql,quall);
    if(minl < minlen)
      return(TRUE);

    fq->seq[minl] = '\0';
    fq->score[minl] = '\0';
    fq->len = minl;
  }

  return(throwOut);
}


boolean fastqItemNext(struct lineFile *lf, struct fastqItem *fq)
{
  int len = 0;
  char *line = NULL;
  boolean gotId = FALSE;
  boolean gotSeq = FALSE;
  boolean gotScoreHead = FALSE;
  boolean gotScore = FALSE;



  //Get the ID line
  gotId = lineFileNextReal(lf,&line); //grab next non-blank non-comment line, should have @
  if (! gotId) return(FALSE);
  line = skipLeadingSpaces(line);
  eraseTrailingSpaces(line);
  if (line[0] != '@')
    errAbort("ERROR: %s doesn't seem to be fastq format.  "
        "Expecting '@' start of line %d, got %c.",
        lf->fileName, lf->lineIx, line[0]);
  //copyID(line,fq->id,1,strlen(line));
  //sprintf(fq->id, "%s",line+1);
  strcpy(fq->id,line+1);
  //get the Sequence

  //get the sequence
  gotSeq = lineFileNext(lf, &line, &len);
  if (! gotSeq) lineFileUnexpectedEnd(lf);
  line = skipLeadingSpaces(line);
  eraseTrailingSpaces(line);
  strcpy((fq->seq),line);
  fq->len = strlen(line);


  //check the score header line
  gotScoreHead = lineFileNext(lf, &line, &len);
  if(!gotScoreHead)
    lineFileUnexpectedEnd(lf);
  line = skipLeadingSpaces(line);
  eraseTrailingSpaces(line);
  if(line[0] != '+')
    errAbort("ERROR: Expected '+' on line %d of file %s, got %s",lf->lineIx,lf->fileName,line);

  //now get the score line
  gotScore = lineFileNext(lf,&line,&len);
  if(! gotScore)
    lineFileUnexpectedEnd(lf);
  line = skipLeadingSpaces(line);
  eraseTrailingSpaces(line);
  strcpy((fq->score),line);
  len = strlen(fq->score);
  if(len != fq->len){
    warn("WARNING: %s has sequences and score strings that "
            "are not the same length. Problem at line %d.\n"
        "Seq: %s\n"
        "Score: %s\n",
            lf->fileName, lf->lineIx,
            fq->seq, fq->score);
  }

  return(TRUE);
}

inline void convPhred33ToPhred64( struct fastqItem *fq )
{
  phred33ToPhred64(fq->score,fq->len);
}//end phred33To64

inline void convPhred64ToPhred33(struct fastqItem *fq)
{
  phred64ToPhred33(fq->score,fq->len);
}

inline void phred33ToPhred64( char * p33, int l )
{
  int i;
  for(i=0;i<l;i++)
  {
    p33[i] = phredToPhred64(phred33ToPhred(p33[i]));
  }
}

inline void phred64ToPhred33( char * p64, int l)
{
  int i;
  for(i=0;i<l;i++)
  {
    p64[i] = phredToPhred33(phred64ToPhred(p64[i]));
  }
}

inline char phredToPhred33( int p )
{
  if (p > MAX_PHRED) p=MAX_PHRED;
  else if (p < MIN_PHRED) p=MIN_PHRED;
  return ((char) (p + 33));
}

inline int phred33ToPhred( char p )
{
  return ((int)p) - 33;
}

inline int phred64ToPhred( char p )
{
  return ((int)p) - 64;
}

inline char phredToPhred64( int p )
{
  if (p > MAX_PHRED) p=MAX_PHRED;
  else if (p < MIN_PHRED) p=MIN_PHRED;
  return ((char) (p + 64));
}

inline int doubleToPhred( double p )
/* formula: -10 log10(p) */
{
  double res = -10.0 * log10(p);
  if (res > MAX_PHRED) res=MAX_PHRED;
  else if (res < MIN_PHRED) res=MIN_PHRED;
  return ((int) (res +0.5));  //guarenteed >= 0
}

inline double phredToDouble( int p )
/* formula: 10^(-p/10)  */
{
  return 1.0/pow(10,((double)p)/10.0);
}

/* Some Functions that can now be made from  a combination of existing functions */

inline double phred33ToDouble( char p )
{
  return phredToDouble(phred33ToPhred(p));
}

inline double phred64ToDouble( char p )
{
  return phredToDouble(phred64ToPhred(p));
}

void strrevi(char *s,int n)
{
  int i=0;
  while (i<n/2)
  {
    *(s+n) = *(s+i);//uses the null character as the temporary storage.
    *(s+i) = *(s + n - i -1);
    *(s+n-i-1) = *(s+n);
    i++;
  }
  *(s+n) = '\0';
}

inline void reverseComplementFastqItem(struct fastqItem *fq){
  //reverse complement the seq, and reverse the quality string
  strrevi(fq->score,fq->len);
  reverseComplement((DNA *)fq->seq,fq->len);
}


