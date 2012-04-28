/**
 * mafPairwiseSyntenyDecay
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
#include "maf.h"

}

#include <map>
#include <string>
#include <set>
#include <vector>
#include <sstream>
#include <iostream>

//#define NDEBUG
#include <assert.h>

using namespace std;


/*
Global variables
 */

class PairAlnInfo {
public:
  string oname;
  int sstart;
  int send;
  int ostart;
  int oend;
  char strand;
  PairAlnInfo( const string _oname,
      const int _sstart, const int _send,
      const int _ostart, const int _oend,
      const char _strand ):
        oname(_oname),
        sstart(_sstart),
        send(_send),
        ostart(_ostart),
        oend(_oend),
        strand(_strand){}
  PairAlnInfo( const char * const _oname,
      const int _sstart, const int _send,
      const int _ostart, const int _oend,
      const char _strand ):
        oname(_oname),
        sstart(_sstart),
        send(_send),
        ostart(_ostart),
        oend(_oend),
        strand(_strand){}
  PairAlnInfo():
        oname(""),
        sstart(0),
        send(0),
        ostart(0),
        oend(0),
        strand('+'){}

  bool isPlus(){
    return(strand == '+');
  }
  bool isRC(){
    return(strand == '-');
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

#define DEF_MIN_LEN (200)
#define DEF_MIN_SCORE (200)

typedef map<int,PairAlnInfo> PairAlnInfoByPos;
typedef map<string, PairAlnInfoByPos > ChromToPairAlnInfoByPos;
ChromToPairAlnInfoByPos pairAlnInfoByPosByChrom;


void usage()
/* Explain usage and exit. */
{
  errAbort(
      (char*)"mafPairwiseSyntenyDecay -- Calculates pairwise syntenic decay from maf alignment containing at least the two specified species.\n"
      "usage:\n"
      "\tmafPairwiseSyntenyDecay [options] [*required options] file1.maf[.gz] ... \n"
      "Options:\n"
      "\t-help\tPrints this message.\n"
      "\t-minScore=NUM\tMinimum MAF alignment score to consider (default 200)\n"
      "\t-minAlnLen=NUM\tMinimum MAF alignment block length to consider (default 200)\n"
      "\t-speciesMain=NAME\t*Name of the main species (exactly as it appears before the '.') in the maf file (REQUIRED)\n"
      "\t-speciesOther=NAME\t*Name of the other species (exactly as it appears before the '.') in the maf file (REQUIRED)\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {(char*)"help",OPTION_STRING},
    {(char*)"minScore",OPTION_INT},
    {(char*)"minAlnLen",OPTION_INT},
    {(char*)"speciesMain",OPTION_STRING},
    {(char*)"speciesOther",OPTION_STRING},
    {NULL, 0}
}; //end options()

/**
 * Main function, takes filenames for paired qseq reads
 * and outputs three files.
 */
int iterateOverAlignmentBlocksAndStorePairInfo(char *fileName, const int minScore, const int minAlnLen, const string speciesMain, const string speciesOther){
  struct mafFile * mFile = mafOpen(fileName);
  struct mafAli * mAli;

  //loop over alignment blocks
  while((mAli = mafNext(mFile)) != NULL){
    struct mafComp *first = mAli->components;
    int seqlen = mAli->textSize;
    //First find and store set of duplicates in this block
    set<string> seen;
    set<string> dups;
    if(mAli->score < minScore || seqlen < minAlnLen)
      goto CLEANUPALNBLOCK;

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
      vector<string> nameSplit1 = split(tmp1,'.');
      string name1 = nameSplit1[0];
      if(dups.count(name1) || (name1 != speciesMain && name1 != speciesOther))
        goto CLEANUPALNBLOCK;
      for(struct mafComp *item2 = item1->next; item2 != NULL; item2 = item2->next){
        string tmp2(item2->src);
        vector<string> nameSplit2 = split(tmp2,'.');
        string name2 = nameSplit2[0];
        if(dups.count(name2) || (name2 != speciesMain && name2 != speciesOther))
          goto CLEANUPALNBLOCK;

        string chr1 = nameSplit1[1];
        string chr2 = nameSplit2[1];
        char strand;
        if(item1->strand == item2->strand)
          strand = '+';
        else
          strand = '-';

        int start1,end1,start2,end2;

        if(item1->strand == '+'){
          start1 = item1->start;
          end1 = start1 + item1->size;
        }else{
          end1 = item1->start;
          start1 = end1 - item1->size;
        }

        if(item2->strand == '+'){
          start2 = item2->start;
          end2 = start2+ item2->size;
        }else{
          end2 = item2->start;
          start2 = end2 - item2->size;
        }

        if(name1 == speciesMain){
          PairAlnInfo aln(chr2,start1,end1,start2,end2,strand);
          pairAlnInfoByPosByChrom[chr1][start1] = aln;
        }else{
          PairAlnInfo aln(chr1,start2,end2,start1,end1,strand);
          pairAlnInfoByPosByChrom[chr2][start2] = aln;
        }

      } //end loop over item2
    } //end loop over item1
CLEANUPALNBLOCK:
    mafAliFree(&mAli);
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
  int minAlnScore = optionInt((char*)"minScore",DEF_MIN_SCORE);
  int minAlnLen = optionInt((char*)"minAlnLen",DEF_MIN_LEN);

  string speciesMain = optionVal((char*)"speciesMain",NULL);
  string speciesOther = optionVal((char*)"speciesOther",NULL);

  if(speciesMain.empty() || speciesOther.empty())
    usage();

  //load the relevant alignment info from the maf(s)
  for(int i = 1; i<argc; i++){
    iterateOverAlignmentBlocksAndStorePairInfo(argv[i], minAlnScore, minAlnLen, speciesMain, speciesOther);
  }

  const int blockSize = 1000;
  const int blockCount = 100;

  int totalWindows[blockCount] = {0};
  int containBreak[blockCount] = {0};

  //we want the fraction of windows of each size that contain a break
  //

  ChromToPairAlnInfoByPos::iterator mainChromItter;
  for(mainChromItter = pairAlnInfoByPosByChrom.begin();
      mainChromItter != pairAlnInfoByPosByChrom.end();
      mainChromItter++){
    //process the alignments shared by this chromosome
    //note that map stores them sorted by begin position
    PairAlnInfoByPos::iterator posIter;
    for(posIter = mainChromItter->second.begin();
        posIter != mainChromItter->second.end();
        posIter++){
      //first check for trivial window (ie our block)
      PairAlnInfo pi1 = posIter->second;
      assert(pi1.send > pi1.sstart);
      assert(pi1.sstart == posIter->first);
      int numBucketsThisWindow = (pi1.send - pi1.sstart) % blockSize;
      int lastContigPos = pi1.send;
      for(int i = 0; i < numBucketsThisWindow && i < blockCount; i++)
        totalWindows[i]++;

      PairAlnInfoByPos::iterator posIterNext = posIter;
      posIterNext++; //move it up one position
      for( /* pre-initialized */; posIterNext != mainChromItter->second.end();
          posIterNext++){
        if((posIterNext->first - posIter->first) >= (blockSize * blockCount)){
          //posIter = posIterNext;
          break;
        }

        PairAlnInfo pi2 = posIterNext->second;

        assert(pi2.sstart == posIterNext->first);
        assert(pi2.send > pi2.sstart);
        assert(pi2.sstart > pi1.sstart);

        if(pi2.oname == pi1.oname){
          int moreToInc = (pi2.send - pi1.sstart) % blockSize;
          lastContigPos = pi2.send;
          for(int i = numBucketsThisWindow; i < moreToInc && i < blockCount; i++)
            totalWindows[i]++;
          numBucketsThisWindow = moreToInc; //so we don't double count
        }else{
          //from the last contiguous position until the start of this
          // block, we have a break somewhere
          int numDiscontigBuckets = (pi2.send - pi1.sstart) % blockSize;
          for(int i = numBucketsThisWindow; i < numDiscontigBuckets && i < blockSize; i++){
            containBreak[i]++;
            totalWindows[i]++;
          }
          numBucketsThisWindow = numDiscontigBuckets;
        }
      }
    }
  }



  cout << "#WindowSize\tNumContainBreak\tNumTotal\t1-(NumContainBreak/NumTotal)" << endl;
  for(int i = 0; i < blockCount; i++){
    cout << (i+1)*blockSize << '\t';
    cout << containBreak[i] << '\t';
    cout << totalWindows[i] << '\t';
    cout << 1.0 - (double(containBreak[i])/double(totalWindows[i])) << endl;
  }


  return(0);
} //end main()


