/**
 * dukeUniquenessFromFasta
 *  Author: John St. John
 *  Date: 4/26/2012
 *  
 *  calculates the mean synteny decay in different range bins
 *
 *
 */

//Kent source C imports
extern "C" {

#include "common.h"
#include "options.h"
#include "fa.h"

}

#include <string>
#include <vector>
#include <tr1/unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace std::tr1;


/*
Global variables
 */

#define DEF_K (35)
#define DEF_HASH_INIT (2000000000)


void usage()
/* Explain usage and exit. */
{
  errAbort(
      (char*)"dukeUniquenessFromFasta -- computes uniqueness per base: 1=1 occurrence, 0.5=2, 0.33=3, 0.25=4, 0=5 or more (or containing sequence ambiguities)\n"
      "usage:\n"
      "\tdukeUniquenessFromFasta [options] input.fa output.wig\n"
      "Options:\n"
      "\t-help\tPrints this message.\n"
      "\t-K=NUM\tK-mer to compute mapability for (default 35)\n"
      "\t-initHashSize=NUM\tInitial size of hash for maybe a slight performance boost, larger numbers consume more memory (default 2000000000)"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {(char*)"help",OPTION_STRING},
    {(char*)"K",OPTION_INT},
    {NULL, 0}
}; //end options()

char op_complement_upper(char c){
  switch(c){
  case 'A':
    return('T');
  case 'T':
    return('A');
  case 'G':
    return('C');
  case 'C':
    return('G');
  default:
    return('N');
  }
}


string getMinKmerFromSeq(char *seq, int K){
  string kmer(seq,K);
  transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper); //uppercase
  //reverse complement
  string rkmer(kmer.rbegin(), kmer.rend());
  transform(rkmer.begin(),rkmer.end(),rkmer.begin(),op_complement_upper);
  return(min(kmer,rkmer));
}


void fillKmerFreqHash(char *fastaFile, const int K, unordered_map<string, unsigned int> &mapToBuild){
  struct lineFile *lf = lineFileOpen(fastaFile,TRUE);
  DNA *seq;
  int seqLen;
  char *seqName;
  while(faMixedSpeedReadNext(lf, &seq, &seqLen, &seqName)){
    //process this fasta entry, store all K-mers without 'N's in the map
    for(int i=0; i<seqLen; i++){
      string kmer = getMinKmerFromSeq(seq+i,K);

      size_t found = kmer.find('N');
      if(found != string::npos){
        i+=int(found);
      }else{
        //no N's
        mapToBuild[kmer]++;
      }
    }
  }
  faFreeFastBuf();
  lineFileClose(&lf);
}

void printDukeUniquenessWiggle(char * fastaFile, char * outFile, const int K, unordered_map<string, unsigned int> &mapToBuild){
  struct lineFile *lf = lineFileOpen(fastaFile,TRUE);
  DNA *seq;
  int seqLen;
  char *seqName;
  ofstream wigFile;
  wigFile.open(outFile);
  while(faMixedSpeedReadNext(lf, &seq, &seqLen, &seqName)){
    //process this fasta entry, store all K-mers without 'N's in the map
    wigFile << "fixedStep  chrom=" << seqName << "  start=0  step=1" << endl;
    for(int i=0; i<seqLen; i++){
      string kmer = getMinKmerFromSeq(seq+i,K);
      size_t found = kmer.find('N');
      if(found != string::npos){
        for(int j=i; j< (i+int(found)); j++){
          /*uniqueness=1=1 occurence,0.5=2,0.33=3,0.25=4,0=5 or more (or containing sequence ambiguities)*/
          wigFile << 0 << endl;
        }
        i+=int(found);
      }else{
        //no N's
        unsigned int count = mapToBuild[kmer];
        /*uniqueness=1=1 occurence,0.5=2,0.33=3,0.25=4,0=5 or more (or containing sequence ambiguities)*/
        switch(count){
        case 1:
          wigFile << 1 << endl;
          break;
        case 2:
          wigFile << 0.5 << endl;
          break;
        case 3:
          wigFile << 0.33 << endl;
          break;
        case 4:
          wigFile << 0.25 << endl;
          break;
        default:
          wigFile << 0 << endl;
          break;
        }
      }
    }
  }
  wigFile.close();
  faFreeFastBuf();
  lineFileClose(&lf);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if(optionExists((char*)"help") || argc != 2){
    usage();
  }
  int K = optionInt((char*)"K",DEF_K);
  int iHash = optionInt((char*)"initHashSize",DEF_HASH_INIT);

  //step 1: load up the hash
  unordered_map<string, unsigned int> kmerFreqHash(iHash);
  fillKmerFreqHash(argv[1],K,kmerFreqHash);

  //step 2: pass over assembly again and print the uniqueness
  printDukeUniquenessWiggle(argv[1],argv[2],K,kmerFreqHash);

  return(0);
} //end main()


