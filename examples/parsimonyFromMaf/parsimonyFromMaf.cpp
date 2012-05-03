/**
 * parsimonyFromMaf
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
#include <algorithm>
#include <vector>
#include <sstream>
#include <iostream>

//#define NDEBUG
#include <assert.h>

using namespace std;


/*
Global variables
 */


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

#define DEF_MIN_LEN (150)
#define DEF_MIN_SCORE (10000)



void usage()
/* Explain usage and exit. */
{
  errAbort(
      (char*)"parsimonyFromMaf -- claculates branch lengths of requested species for use in building parsimony trees.\n"
      "usage:\n"
      "\tparsimonyFromMaf [options] file1.maf[.gz] ... \n"
      "Required Options:\n"
      "\t-species=LST\tWhere LST is a comma separated list of genome names (between 2 and 32), like 'allMis1,croPor0'\n"
      "Options:\n"
      "\t-help\tPrints this message.\n"
      "\t-minScore=NUM\tMinimum MAF alignment score to consider (default 10000)\n"
      "\t-minAlnLen=NUM\tMinimum MAF alignment block length to consider (default 150)\n"

  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {(char*)"help",OPTION_STRING},
    {(char*)"species",OPTION_STRING},
    {(char*)"minScore",OPTION_INT},
    {(char*)"minAlnLen",OPTION_INT},
    {NULL, 0}
}; //end options()

/**
 * Function to iterate over a block that contains all of the species we are interested in.
 *
 *
 */
void incrementBranchCountsFromFilteredAlnBlock(
    map<unsigned int, string> &filteredSpeciesBitFlagToDNA,
    map<unsigned int,unsigned int> &branchMutationCounts){

  char DNAc[] = {'A','C','G','T'};
  const set<char> DNA(DNAc, DNAc + 4);

  vector<unsigned int> iteratorOrderedSpeciesBitFlags;



  //find columns that are all in [a,c,g,t]
  bool first = true;
  set<int> good_sites;
  for(map<unsigned int, string>::iterator it = filteredSpeciesBitFlagToDNA.begin();
      it != filteredSpeciesBitFlagToDNA.end(); ++it){
    string seq = it->second;
    iteratorOrderedSpeciesBitFlags.push_back(it->first);
    for(int i = 0; i < seq.size();i++){
      if(DNA.count(seq[i])){
        if(first){
          good_sites.insert(i);
        }
      }else{
        if(!first){
          good_sites.erase(i);
        }
      }
    }//end loop over sites
    if(first)
      first = false;
  }


  //loop over the good alignment columns
  //increasing the counts

  for(set<int>::iterator it = good_sites.begin();
      it != good_sites.end();++it){
    int site = *it;
    vector<char> alnColumn;
    set<char> seen;
    for(map<unsigned int, string>::iterator sit = filteredSpeciesBitFlagToDNA.begin();
        sit != filteredSpeciesBitFlagToDNA.end(); ++sit){
      alnColumn.push_back(sit->second[site]);
      seen.insert(sit->second[site]);
    }
    if(seen.size() > 2)
      continue; //interested in bialelic sites
    if(seen.size() == 1){
      branchMutationCounts[0]++;
    }else{
      //bialelic site
      unsigned int alleleAc = 0;
      unsigned int alleleBc = 0;
      set<char>::iterator tmpit = seen.begin();
      char alleleA = *tmpit;

      for(int i = 0; i < alnColumn.size(); i++){
        if(alnColumn[i] == alleleA){
          alleleAc |= iteratorOrderedSpeciesBitFlags[i];
        }else{
          alleleBc |= iteratorOrderedSpeciesBitFlags[i];
        }
      }//end for loop over column
      branchMutationCounts[min(alleleAc,alleleBc)]++;

    }


  }//end loop over columns


}




/**
 * Main maf iteration function
 */
int iterateOverAlignmentBlocks(char *fileName, int minScore, int minAlnLen,
    map<string, unsigned int> speciesToBitFlag,
    set<string> species,
    map<unsigned int,unsigned int> &branchMutationCounts){

  struct mafFile * mFile = mafOpen(fileName);
  struct mafAli * mAli;
  //loop over alignment blocks
  while((mAli = mafNext(mFile)) != NULL){
    struct mafComp *first = mAli->components;
    int seqlen = mAli->textSize;
    //First find and store set of duplicates in this block
    set<string> seen;
    set<string> dups;
    set<string> unique;
    map<unsigned int, string> filteredSpeciesBitFlagToString;
    if(mAli->score < minScore || seqlen < minAlnLen){
      //free here and pre-maturely end
      mafAliFree(&mAli);
      continue;
    }

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
    //put the set difference into unique
    set_difference(seen.begin(),seen.end(),dups.begin(),dups.end(),inserter(unique,unique.end()));

    //check that unique carries everything in species
    set<string> in_both;
    set_intersection(unique.begin(),unique.end(),species.begin(),species.end(),inserter(in_both,in_both.end()));
    if(species.size() != in_both.size()){
      continue;
    }


    for(struct mafComp *item = first; item->next != NULL; item = item->next){
      //stop one before the end
      string tmp1(item->src);
      string name1 = split(tmp1,'.')[0];
      if(species.count(name1)){
        string dna = string(item->text);
        transform(dna.begin(), dna.end(), dna.begin(), ::toupper); //make it upper case
        filteredSpeciesBitFlagToString[speciesToBitFlag[name1]] = dna;
      }
    } //end loop over items

    incrementBranchCountsFromFilteredAlnBlock(filteredSpeciesBitFlagToString, branchMutationCounts);



    mafAliFree(&mAli);
  }//end loop over alignment blocks

  mafFileFree(&mFile);
  return(0);
}

void writeTreeBranches(map<string, unsigned int>speciesToBitFlag,
    map<unsigned int, unsigned int> branchMutationCounts){
  for(map<unsigned int, unsigned int>::iterator it = branchMutationCounts.begin();
      it != branchMutationCounts.end(); ++it){
    //write out one of the groups of species that have this mutation
    unsigned int speciesBitSet = it->first;
    unsigned int count = it->second;
    cout << '(';
    bool first = true;
    for(map<string,unsigned int>::iterator bsit = speciesToBitFlag.begin();
        bsit != speciesToBitFlag.end(); ++bsit){
      if(bsit->second & speciesBitSet){
        if(!first)
          cout << ',';
        cout << bsit->first;
        if(first)
          first = false;
      }
    }
    cout << "):" << count << endl;
  }
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
  string speciesS = string(optionVal((char*)"species",NULL));

  vector<string> speciesV = split(speciesS,',');
  set<string> species;
  map<string, unsigned int>speciesToBitFlag;

  for(vector<string>::iterator it = speciesV.begin();
      it != speciesV.end(); ++it){
    if((*it).size()){
      species.insert(*it);
    }
  }

  if(species.size() < 2 || species.size() > 32){
    cerr << "Must request between 2 and 32 species" << endl;
    usage();
  }

  unsigned int i = 0;
  for(set<string>::iterator it = species.begin();
      it != species.end(); ++it){
    unsigned int id = 1 << i;
    speciesToBitFlag[*it] = id;
    ++i;
  }


  map<unsigned int, unsigned int> branchMutationCounts;

  for(int i = 1; i<argc; i++){
    iterateOverAlignmentBlocks(argv[i], minAlnScore, minAlnLen,
        speciesToBitFlag, species, branchMutationCounts);
  }

  writeTreeBranches(speciesToBitFlag, branchMutationCounts);

  return(0);
} //end main()


