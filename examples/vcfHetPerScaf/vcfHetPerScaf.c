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
#include "uthash.h"

#define NUMCOLSVCF 20
#define VCFSAMPLESTART 9
#define NUMSAMPLES ( NUMCOLSVCF - VCFSAMPLESTART )

/*
Global variables
 */
bool verboseOut = false;
int sampleCount = 0;
struct sizes_hash {
    char name[25];                /* key (structure POINTS TO string */
    unsigned length;
    unsigned numHets[NUMSAMPLES];
    UT_hash_handle hh;         /* makes this structure hashable */
};


struct sizes_hash *chrSizes = NULL;

void usage()
/* Explain usage and exit. */
{
  errAbort(
      "vcfHetPerScaf -- Outputs the heterozygosity per sample per chromosome.\n"
      "usage:\n"
      "\vcfHetPerScaf [required options] [options] \n"
      "\n**required** options:\n"
      "\t-vcf=FILE\tFile name holding the vcf file to parse.\n"
      "\t-sizes=FILE\tFile name holding the chromosome sizes for calculating heterozygosity.\n"
      "\noptions:\n"
      "\t-verbose\tOutput verbose debug messages to stderr.\n\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"vcf",OPTION_STRING},
    {"sizes",OPTION_STRING},
    {"numSamples",OPTION_INT},
    {"verbose",OPTION_BOOLEAN},
    {NULL, 0}
}; //end options()


int fillSizes(char *sizes)
/**
 * Fills the hash table of chrom sizes
 *
 *
 *
 */
{
  struct lineFile *lf;
  lf = lineFileOpen(sizes,TRUE);
  char *line;
  struct sizes_hash *s;
  while(lineFileNextReal(lf,&line)){
    char *split[2];
    chopByWhite(line,split,2);
    s = malloc(sizeof(struct sizes_hash));
    strcpy(s->name, split[0]);
    s->length = atoi(split[1]);
    int i;
    for(i=0;i<NUMSAMPLES;i++)
      s->numHets[i] = 0;
    HASH_ADD_STR(chrSizes,name,s);
  }
  return 0;
}

int hetPerSeq(char *vcf, int numSamples)
/**
 *
 *
 * ##fileformat=VCFv4.1
##samtoolsVersion=0.1.18 (r982:295)
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">
##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configuration in the trio">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.">
##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.">
##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ../bwa.PE.Braun/amiss_v0.1c27.pe.braun.bam      ../bwa.PE.180.Feb2011/amiss_v0.1c27.pe.180.feb2011.bam
scaffold1       10504   .       C       G       169     .       DP=13;VDB=0.0416;AF1=1;AC1=4;DP4=0,0,7,4;MQ=60;FQ=-44.9 GT:PL:GQ        1/1:123,18,0:36 1/1:82,15,0:33
scaffold1       10509   .       C       T       139     .       DP=14;VDB=0.0349;AF1=1;AC1=4;DP4=0,0,7,3;MQ=60;FQ=-40.5 GT:PL:GQ        1/1:123,21,0:35 1/1:52,9,0:23
scaffold1       11211   .       G       A       133     .       DP=14;VDB=0.0278;AF1=1;AC1=4;DP4=0,0,11,1;MQ=44;FQ=-43.5        GT:PL:GQ        1/1:71,12,0:29  1/1:98,24,0:41
scaffold1       11226   .       T       G       163     .       DP=10;VDB=0.0272;AF1=1;AC1=4;DP4=0,0,6,2;MQ=39;FQ=-39.7 GT:PL:GQ        1/1:90,9,0:22   1/1:109,15,0:28
scaffold1       11233   .       A       G       125     .       DP=8;VDB=0.0183;AF1=1;AC1=4;DP4=0,0,6,2;MQ=34;FQ=-39.7  GT:PL:GQ        1/1:69,9,0:22   1/1:92,15,0:28
scaffold1       15972   .       C       T       67.2    .       DP=26;VDB=0.0384;AF1=0.2606;AC1=1;DP4=9,8,5,4;MQ=60;FQ=68.5;PV4=1,5.4e-08,1,1   GT:PL:GQ        0/1:100,0,255:98        0/0:0,12,118:14
scaffold1       18094   .       C       G       51      .       DP=35;VDB=0.0252;AF1=0.25;AC1=1;DP4=14,17,2,1;MQ=59;FQ=52.2;PV4=0.59,0.089,1,1  GT:PL:GQ        0/0:0,81,255:83 0/1:84,0,109:82
scaffold1       22776   .       T       C       158     .       DP=17;VDB=0.0284;AF1=0.3335;AC1=1;DP4=4,5,4,4;MQ=60;FQ=160;PV4=1,0.012,1,0.039  GT:PL:GQ        0/1:190,0,203:99        0/0:0,3,40:5
scaffold1       25400   .       T       C       47      .       DP=26;VDB=0.0208;AF1=0.25;AC1=1;DP4=14,9,0,3;MQ=60;FQ=48.2;PV4=0.085,1,1,0.21   GT:PL:GQ        0/0:0,57,255:59 0/1:80,0,117:78
scaffold2       12966   .       G       C       25.2    .       DP=5;VDB=0.0095;AF1=1;AC1=4;DP4=0,0,2,1;MQ=52;FQ=-32.3  GT:PL:GQ        1/1:37,6,0:12   1/1:22,3,0:9
scaffold2       12973   .       G       T,C     40      .       DP=6;VDB=0.0041;AF1=0.8893;AC1=3;DP4=0,1,2,2;MQ=55;FQ=-26.3;PV4=1,1,0.34,0.26   GT:PL:GQ        0/1:36,2,5,30,0,34:4    1/1:39,3,0,39,3,39:6
 */
{
  struct lineFile *lf;
  lf = lineFileOpen(vcf,TRUE);
  char *line;
  struct sizes_hash *s;
  //int numLines = 9+numSamples; //there are 9 + the number of samples columns in the vcf file
  char lastSeqId[200];
  strcpy(lastSeqId,"23121NoTASeqUeNCE22");


  while(lineFileNextReal(lf,&line)){
    // if line starts with "#" skip it
    char *split[NUMCOLSVCF];
    if(line[0] == '#')
      continue;
    int numRecords = chopByWhite(line,split,NUMCOLSVCF);
    if(sampleCount == 0){
      sampleCount = numRecords-VCFSAMPLESTART;
    }
    char *seqId = split[0];
    if(strcmp(seqId,lastSeqId)!= 0){
      strcpy(lastSeqId,seqId);
      HASH_FIND_STR(chrSizes,seqId,s);
    }
    // now check for het sites and increment count
    int i;
    for(i=VCFSAMPLESTART; i < numRecords; i++){
      char sampleStr[200];
      strcpy(sampleStr,split[i]);
      //samples look like: 1/1:37,6,0:12
#define PARTSSIZE 10
      char *parts[PARTSSIZE];
      int numSections = chopString(sampleStr,":",parts,PARTSSIZE);
      int j = 0;
      char *allele;
      for(j = 0;j<numSections;j++){
        if(strlen(parts[j]) == 3 && parts[j][1] == '/'){
          allele = parts[j];
          break;
        }
      }

      if((allele[0] == '0' && allele[2] != '0')||
          (allele[0] == '1' && allele[2] != '1')||
          (allele[0] == '2' && allele[2] != '2')){
        //The alleles don't match
        s->numHets[i-VCFSAMPLESTART]++;

      }

    }

  }
  return 0;
}

int printHeterozygosity(FILE *out){
  struct sizes_hash *s,*tmp;
  HASH_ITER(hh, chrSizes, s, tmp) {
//      HASH_DEL(chrSizes,s);  /* delete; users advances to next */
//      free(s);            /* optional- if you want to free  */
    int i;
    fprintf(out,"%s",s->name);
    for(i=0;i<sampleCount;i++){
      fprintf(out,"\t%lf",(double)s->numHets[i]/(double)s->length);
    }
    fprintf(out,"\n");


  }
  return 0;
}


int main(int argc, char *argv[])
/* Process command line. */
{
  char *vcf = NULL;
  char *sizes = NULL;
  int numSamples = 0;
  optionInit(&argc, argv, options);

  verboseOut = optionExists("verbose");
  vcf = optionVal("vcf",NULL);
  sizes = optionVal("sizes",NULL);
  numSamples = optionInt("numSamples",0);

  if (vcf == NULL || sizes == NULL){
//    if(numSamples == 0)
//       fprintf(stderr, "Error: must provide the number of samples in the vcf file >= 1\n");
    if(vcf == NULL)
      fprintf(stderr,"Error: Must provide vcf file\n");
    if(sizes == NULL)
      fprintf(stderr,"Error: Must provide sequence size file\n");
    usage();
  }
  fillSizes(sizes);
  hetPerSeq(vcf,numSamples);
  printHeterozygosity(stdout);
  return 0;
} //end main()


