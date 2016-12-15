#include <cmath>
#include "genomePosition.h"
#include "pFile.h"
#include "Error.h"

#define DEFAULT_NUM_CHRS 25

genomePosition genomePos;

genomePosition::genomePosition() {
  // hg19 by default
  const char* _chrs[DEFAULT_NUM_CHRS] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"};
  uint64_t _szs[DEFAULT_NUM_CHRS] = {
    249250621
    ,243199373
    ,198022430
    ,191154276
    ,180915260
    ,171115067
    ,159138663
    ,146364022
    ,141213431
    ,135534747
    ,135006516
    ,133851895
    ,115169878
    ,107349540
    ,102531392
    ,90354753 
    ,81195210 
    ,78077248 
    ,59128983 
    ,63025520 
    ,48129895 
    ,51304566 
    ,155270560
    ,59373566 
    ,16569
  };

  uint64_t cum = 0;
  for(int i=0; i < DEFAULT_NUM_CHRS; ++i) {
    mChrs[_chrs[i]] = (int)chrs.size();
    chrs.push_back(_chrs[i]);
    szChrs.push_back(_szs[i]);
    cum += _szs[i];
    cumChrs.push_back(cum);
  }
}

void genomePosition::loadFastaIndex(const char* fai) {
  pFile pf(fai);
  const char* line = NULL;
  std::vector<std::string> tokens;
  chrs.clear();
  szChrs.clear();
  cumChrs.clear();
  mChrs.clear();
  uint64_t sz = 0, cum = 0;
  while( ( line = pf.getLine() ) != NULL ) {
    pFile::tokenizeLine(line," \t\r\n", tokens);
    if ( tokens.size() < 2 ) error("Cannot parse fai file %s - observed line is '%s'",fai,line);
    sz = atoi(tokens[1].c_str());
    mChrs[tokens[0]] = (int)chrs.size();
    chrs.push_back(tokens[0]);
    szChrs.push_back(sz);
    cum += sz;
  }
}

uint64_t genomePosition::fromPos(const std::string& chr, int pos) {
  std::map<std::string,int>::iterator it = mChrs.find(chr);
  if ( it == mChrs.end() ) {
    warning("genomePosition::fromPos() Cannot find chromosome %s",chr.c_str());
    return 0;
  }
  else {
    return ( cumChrs[it->second] - szChrs[it->second] + (uint64_t)pos );
  }
}

uint64_t genomePosition::fromPos(const char* chr, int pos) {
  std::string schr(chr);
  return fromPos(schr,pos);
}

std::pair<std::string,int> genomePosition::toPos(uint64_t upos) {
  for(int i=0; i < (int)chrs.size(); ++i) {
    if ( upos < cumChrs[i] ) {
      return std::pair<std::string,int> (chrs[i],(int)(upos-(cumChrs[i]-szChrs[i])));
    }
  }
  error("UPOS %llu exceeds the size of the genome",upos);
  return std::pair<std::string,int> ("NA",0);
}

