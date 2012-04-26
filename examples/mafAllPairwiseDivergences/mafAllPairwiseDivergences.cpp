/**
 * mafAllPairwiseDivergences
 *  Author: John St. John
 *  Date: 4/26/2012
 *  
 *  Program to calculate pairwise divergence from 1-1 alignment of
 *  MAF blocks.
 *
 *
 */

//Kent source C imports
extern "C" {

#include "common.h"
#include "options.h"
#include "maf.h"

}

#include <map>
#include <string>
#include <set>
#include <vector>
#include <sstream>
#include <iostream>

//#define NDEBUG
#import <assert.h>

using namespace std;


/*
Global variables
 */

class NucDivergence {
private:
  unsigned int data[4][4];
  int base2inx(const char base){
    switch(base){
    case 'A':
    case 'a':
      return(0);
    case 'C':
    case 'c':
      return(1);
    case 'G':
    case 'g':
      return(2);
    case 'T':
    case 't':
      return(3);
    default:
      return(4);
    }
  }

  bool isValid(const char base){
    switch(base){
        case 'A':
        case 'a':
        case 'C':
        case 'c':
        case 'G':
        case 'g':
        case 'T':
        case 't':
          return(true);
        default:
          return(false);
        }
  }
public:

  NucDivergence(){
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
        data[i][j] = 0;
  }

  void incVal(const char c1, const char c2){
    if(isValid(c1) && isValid(c2))
      data[base2inx(c1)][base2inx(c2)]++;
  }
  unsigned int getCount(const char c1, const char c2){
    return(data[base2inx(c1)][base2inx(c2)]);
  }
  unsigned int getTs(){
    return(data[0][2]+data[2][0]+data[1][3]+data[3][1]);
  }
  unsigned int getTv(){
    return(data[0][3]+data[3][0]+data[1][0]+data[0][1]+data[1][2]+data[2][1]+data[3][2]+data[2][3]);
  }
  double getTsTv(){
    return(double(getTs())/double(getTv()));
  }
  double getPID(){
    unsigned int total = 0;
    unsigned int ident = 0;
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++){
        total += data[i][j];
        if(i == j)
          ident += data[i][j];
      }
    return(double(ident)/double(total));
  }
};

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return(elems);
}


vector<string> split(const string &s, char delim) {
    vector<string> elems;
    return(split(s, delim, elems));
}

#define DEF_MIN_LEN (30)
typedef map<string,NucDivergence> stringToNucDivergence;
typedef map<string, stringToNucDivergence > stringPiarToNucDivergence;
stringPiarToNucDivergence pairwiseDivergence;


void usage()
/* Explain usage and exit. */
{
  errAbort(
      (char*)"mafAllPairwiseDivergences -- Calculates all pairiwise nucleotide difference spectra across each 1-1 alignment block in the list of maf files.\n"
      "usage:\n"
      "\tmafAllPairwiseDivergences [options] file1.maf[.gz] ... \n"
      "Options:\n"
      "\t-help\tPrints this message.\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {(char*)"help",OPTION_STRING},
    {NULL, 0}
}; //end options()

/**
 * Main function, takes filenames for paired qseq reads
 * and outputs three files.
 */
int iterateOverAlignmentBlocks(char *fileName){
  struct mafFile * mFile = mafOpen(fileName);
  struct mafAli * mAli;

  //loop over alignment blocks
  while((mAli = mafNext(mFile)) != NULL){
    struct mafComp *first = mAli->components;

    //First find and store set of duplicates in this block
    set<string> seen;
    set<string> dups;
    for(struct mafComp *item = first; item != NULL; item = item->next){
      string tmp(item->src);
      string tname = split(tmp,'.')[0];
      if(seen.count(tname)){
        //seen this item
        dups.insert(tname);
      }else{
        seen.insert(tname);
      }
    }
    for(struct mafComp *item1 = first; item1->next != NULL; item1 = item1->next){
      //stop one before the end
      string tmp1(item1->src);
      string name1 = split(tmp1,'.')[0];
      if(dups.count(name1))
        continue;
      for(struct mafComp *item2 = item1->next; item2 != NULL; item2 = item2->next){
        string tmp2(item2->src);
        string name2 = split(tmp2,'.')[0];
        if(dups.count(name2))
          continue;

        assert(name1 != name2);
        int seqlen = strlen(item1->text);
        assert(seqlen == strlen(item2->text));


        //make the first item of the index
        //the one that is alphanumerically less
        //for indexing consistency
        string namePos0;
        string namePos1;
        char *seq0;
        char *seq1;
        if(name1 < name2){
          namePos0 = name1;
          seq0 = item1->text;
          namePos1 = name2;
          seq1 = item2->text;
        }else{
          namePos0 = name2;
          seq0 = item2->text;
          namePos1 = name1;
          seq1 = item1->text;
        }

        assert(namePos0 < namePos1);

        for( int i = 0; i < seqlen; i++ ){
          char c0,c1;
          c0 = seq0[i];
          c1 = seq1[i];
          //increments only if both chars are valid
          pairwiseDivergence[namePos0][namePos1].incVal(c0,c1);
        } //end loop over pairwise block
      } //end loop over item2
    } //end loop over item1
  }//end loop over alignment blocks

  mafFileFree(&mFile);
  return(0);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if(optionExists((char*)"help") || argc <= 1){
    usage();
  }

  for(int i = 1; i<argc; i++){
    iterateOverAlignmentBlocks(argv[i]);
  }

  //print the pairwise divergence
  //stringPiarToNucDivergence pairwiseDivergence;
  stringPiarToNucDivergence::iterator outer;
  stringToNucDivergence::iterator inner;
  char nucLst[4] = {'A','C','G','T'};
  cout << "#sample1\tsample2";
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      char ci = nucLst[i];
      char cj = nucLst[j];
      cout << '\t' << ci << cj;
    }
  }
  cout << "\tts\ttv\tts/tv\tpid"  << endl;

  for(outer = pairwiseDivergence.begin(); outer != pairwiseDivergence.end(); outer++){
    for(inner = outer->second.begin(); inner != outer->second.end(); inner++){
      string ostr = outer->first;
      string istr = inner->first;
      NucDivergence nd = inner->second;
      cout << ostr << '\t' << istr;
      for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++){
          char ci = nucLst[i];
          char cj = nucLst[j];
          cout << '\t' << nd.getCount(ci,cj);
        }//end loop over things
      cout << '\t' << nd.getTs();
      cout << '\t' << nd.getTv();
      cout << '\t' << nd.getTsTv();
      cout << '\t' << nd.getPID();
      cout << endl;
    }
  }



  return(0);
} //end main()


