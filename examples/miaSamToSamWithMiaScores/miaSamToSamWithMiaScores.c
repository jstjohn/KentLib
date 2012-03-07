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
#include "sam.h"
#include "bamFile.h"
#include "params.h"
#include "pssm.h"
#include "dnaseq.h"
#include "dnaLoad.h"
#include "hash.h"




/**
 * Global Variables and definitions
 *
 */ 





/**
 * Command line options
 */ 

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "miaSamToSamWithMiaScores -- Takes a Sam file, and adds in the Mia alignment score into the AS auxiliary tag.\n"
      "\nUsage: miaSamToSamWithMiaScores [options] scoring_matrix.txt reference.[fa,2bit,nib] file_in.[s,b]am file_out.[b,s]am\n"
      "\noptions:\n"
      "\t-fai=FILE\tIf a sam file is provided that does not have a header.\n"
      "\t-help\tWrites this help to the screen.\n"
      );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"help",OPTION_BOOLEAN},
    {"fai",OPTION_STRING},
    {NULL, 0}
}; //end options()



/**
 * Program's functions
 *
 */



inline unsigned short base2inx(const char base) {
  switch (base) {
  case 'A':
  case 'a':
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g':
    return 2;
  case 'T':
  case 't':
    return 3;
  default:
    return 4;
  }
}







void addAlignScoreToBamAux(samfile_t *bamFileIn, samfile_t *bamFileOut, PSSMP fpssmp, PSSMP rpssmp, struct hash* dnaHash){

  bam_header_t *header = bamFileIn->header;
  bam1_t *b;
  bam1_core_t *c;
  unsigned int *cigar;
  int ***currPSSMP;
  uint8_t *packedQSeq;
  int lastTID = -272823; //this number should never happen
  struct dnaSeq *refSeq = NULL;
  DNA *refDna;
  const int bamAuxAppendIntSize = sizeof(int)/sizeof(uint8_t);

  b = bam_init1();
  c = &b->core;

  while(samread(bamFileIn, b)>=0){
    if(c->tid < 0 || (c->flag & BAM_FUNMAP))
      continue; //skip unmapped reads


    int alnScore = 0;


    if(c->tid != lastTID){
      //get the new reference fasta!
      refSeq = (struct dnaSeq *) hashMustFindVal(dnaHash, header->target_name[c->tid]);
      refDna = refSeq->dna;
      lastTID = c->tid;
    }

    //Determine if we should use forward or reverse PSSMP
    if(c->flag & BAM_FREVERSE){
      currPSSMP = (int ***)(rpssmp->sm);
    }else{
      currPSSMP = (int ***)(fpssmp->sm);
    }

    packedQSeq = bam1_seq(b);

    int i,k,l,op,end;
    unsigned int tpos = c->pos;
    k = 0;
    op = cigar[k] & BAM_CIGAR_MASK;
    l = cigar[k] >> BAM_CIGAR_SHIFT;
    cigar = bam1_cigar(b);
    i = 0;
    int seqlen = c->l_qseq;

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
        alnScore -= (GOP + (GEP*l));
        i += l;
        break;

      case BAM_CDEL: //deletion relative to the reference
        alnScore -= (GOP + (GEP*l));
        tpos += l;
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
          char qbase = bam_nt16_rev_table[bam1_seqi(packedQSeq, i)];
          char rbase = refDna[tpos];
          //PSSM[position][reference base][query base]
          alnScore += currPSSMP[find_sm_depth(i,seqlen)][base2inx(rbase)][base2inx(qbase)];
        }
        break;

      default:
        break;
      }//done switching over alignment

    }//done loop through cigar alignment things

    //now we need to add the alignment score onto the bam aux data
    // AS:i:[alnScore]
    //see bam_aux_append usage in ../../thirdparty/samtools/bam_md.c
    bam_aux_append(b,"AS",'i',bamAuxAppendIntSize,(uint8_t*)&alnScore);
    //8*4=32=sizeof int, that's why heng li passes 4 as the len of the
    //aux data (after the tags). I tried to clarify the reasoning behind this
    //with the following:
    //int bamAuxAppendIntSize = sizeof(int)/sizeof(uint8_t);
    samwrite(bamFileOut,b);

  }//done with this bam item

}



int main(int argc, char *argv[])
/* Process command line. */
{

  optionInit(&argc, argv, options);
  boolean help = optionExists("help");
  char *refLst = optionVal("fai",NULL);
  if(help || argc != 5) usage();

  //read in our scoring matrices
  PSSMP fpssmp = read_pssm(argv[1]);
  PSSMP rpssmp = revcom_submat(fpssmp);


  samfile_t *bamFileIn = 0;
  samfile_t *bamFileOut = 0;
  char *inf = argv[3];
  int inlen = strlen(inf);
  char *outf = argv[4];
  int outlen = strlen(outf);

  if(inlen < 3 || toupper(inf[inlen-3]) == 'S'){
    bamFileIn = samopen(inf,"r",refLst);
  }else{
    bamFileIn = samopen(inf,"rb",refLst);
  }
  if(outlen < 3 || toupper(outf[outlen-3]) == 'S'){
    samfile_t *bamFileOut = samopen(outf,"w",refLst);
  }else{
    samfile_t *bamFileOut = samopen(outf,"wb",refLst);
  }
  //load our reference genome
  struct dnaSeq *refList = dnaLoadAll(argv[2]);
  struct hash *refListHash = dnaSeqHash(refList);

  //do the things here
  addAlignScoreToBamAux(bamFileIn, bamFileOut, fpssmp, rpssmp, refListHash);


  samclose(bamFileIn);
  samclose(bamFileOut);

  return 0;
} //end main()


