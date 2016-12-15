#ifndef __GENOME_POSITION_H__
#define __GENOME_POSITION_H__

#include <utility>
#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <stdint.h>
#include "Error.h"

class genomePosition 
{
 protected:
  std::map<std::string,int> mChrs;
  std::vector<std::string> chrs;
  std::vector<uint64_t> szChrs;
  std::vector<uint64_t> cumChrs;

public:
  genomePosition();
  void loadFastaIndex(const char* fai);
  uint64_t fromPos(const char* chr, int pos);
  uint64_t fromPos(const std::string& chr, int pos);
  std::pair<std::string,int> toPos(uint64_t upos);
  uint64_t genomeSize() { return cumChrs.back(); }
};

extern genomePosition genomePos;

class genomePosLocus {
 public:
  uint64_t ubeg1;
  uint64_t uend0;

 genomePosLocus(uint64_t b, uint64_t e) : ubeg1(b), uend0(e) {}

  genomePosLocus(const char* c, int b, int e) {
    ubeg1 = genomePos.fromPos(c,b);
    uend0 = genomePos.fromPos(c,e);
  }

  genomePosLocus(const char* region) {
    const char* pcolon = strchr(region,':');
    const char* pminus = strchr(pcolon+1,'-');
    if ( ( pcolon == NULL ) || ( pminus == NULL ) ) 
      error("Cannot parse %s in genomePosLocus::genomePosLocus()",region);
    std::string chrom = std::string(region,0,pcolon-region);
    ubeg1 = genomePos.fromPos(chrom,atoi(pcolon+1));
    uend0 = genomePos.fromPos(chrom,atoi(pminus+1));
  }

  bool operator< (const genomePosLocus& l) const {
    if ( ubeg1 == l.ubeg1 ) {
      return ( uend0 < l.uend0 );
    }
    else {
      return ( ubeg1 < l.ubeg1 );
    }
  }

  uint64_t length() const {
    //error("%lld\t%lld\t%lld",uend0,ubeg1,uend0-ubeg1+1);
    if ( uend0 < ubeg1 ) {
      error("genomePosLocus::length() - %llu\t%llu",ubeg1,uend0);
    }
    return (uend0-ubeg1+1); 
  }

  bool overlaps (const genomePosLocus& l) const {
    if ( ( ubeg1 <= l.uend0 )  && ( l.ubeg1 <= uend0 ) ) {
      return true;
    }
    else {
      return false;
    }
  }

  bool merge (const genomePosLocus& l) {
    if ( ( ubeg1-1 <= l.uend0 )  && ( l.ubeg1-1 <= uend0 ) ) {
      if ( l.ubeg1 < ubeg1 ) ubeg1 = l.ubeg1;
      if ( l.uend0 > uend0 ) uend0 = l.uend0;
      return true;
    }
    else {
      return false;
    }
  }

  bool contains0(uint64_t upos0) const { 
    return contains1(upos0+1);
  }

  bool contains0(const char* chr, int pos0) const { 
    return contains1(chr,pos0+1); 
  }

  bool contains1(const char* chr, int pos1) const {
    return contains1(genomePos.fromPos(chr,pos1));
  }

  bool contains1(uint64_t upos1) const {
    return ( ( upos1 >= ubeg1 ) && ( upos1 <= uend0 ) );
  }
};

class genomePosLoci {
 public:
  std::set<genomePosLocus> loci;
  std::set<genomePosLocus>::iterator it;
  bool overlapResolved;

  genomePosLoci() : overlapResolved(false) {}

  void rewind() { it = loci.begin(); }
  bool next() { ++it; return ( it != loci.end() ); }
  bool isend() { return ( it == loci.end() ); }
  const genomePosLocus& currentLocus() { return (*it); }
  bool empty() { return loci.empty(); }

  int numLocus() const { return (int)loci.size(); }

  bool add(const char* chr, int beg1, int end0) {
    overlapResolved = false;
    if ( beg1 > end0 ) error("genomePosLoci::add() Negative interval %s:%d-%d",chr,beg1,end0);
    return loci.insert(genomePosLocus(chr,beg1,end0)).second;
  }

  bool add(const char* region) {
    overlapResolved = false;
    return loci.insert(genomePosLocus(region)).second;
  }
  
  bool add(uint64_t ubeg1, uint64_t uend0) {
    overlapResolved = false;
    if ( ubeg1 > uend0 ) error("genomePosLoci::add() Negative interval %llu-%llu",ubeg1,uend0);
    return loci.insert(genomePosLocus(ubeg1,uend0)).second;    
  }

  int resolveOverlaps() {
    if ( !overlapResolved ) {
      std::set<genomePosLocus>::iterator it;
      std::set<genomePosLocus>::iterator prev;
      int numMerged = 0;
      for(it = loci.begin(); it != loci.end(); ++it) {
	if ( it != loci.begin() ) {
	  if ( prev->overlaps(*it) ) {
	    // if overlaps, erase both and insert merged one
	    genomePosLocus locus = *prev;
	    locus.merge(*it);
	    loci.erase(it);
	    loci.erase(prev);
	    prev = it = loci.insert(locus).first;
	    ++numMerged;
	  }
	  else {
	    prev = it;
	  }
	}
	else {
	  prev = it;
	}
      }
      overlapResolved = true;
      return numMerged;
    }
    else {
      return 0;
    }
    return 0;
  }

  uint64_t totalLength() const {
    uint64_t sz = 0;
    std::set<genomePosLocus>::iterator it;
    for(it = loci.begin(); it != loci.end(); ++it) {
      sz += it->length();
      //notice("%llu\t%llu\t%llu\t%llu",it->ubeg1,it->uend0,it->length(),sz);
    }
    return sz;
  }

  bool moveTo(uint64_t upos1) {
    genomePosLocus locus(upos1,upos1);
    it = loci.lower_bound(locus);
    if ( it == loci.begin() ) { // do nothing
      return (it->contains1(upos1));
    }
    else if ( it == loci.end() ) {
      std::set<genomePosLocus>::iterator i = it;
      --i;
      if ( i->contains1(upos1) ) { it = i; return true; }
      else { return false; }
    }
    else {
      if ( it->contains1(upos1) ) return true;
      else {
	std::set<genomePosLocus>::iterator i = it;
	--i;
	if ( i->contains1(upos1) ) { it = i; return true; }
	else { return false; }
      }
    }
  }

  bool moveTo(const char* chr, int pos1) {
    return moveTo(genomePos.fromPos(chr,pos1));
  }

  bool contains1(uint64_t upos1) {
    genomePosLocus locus(upos1,upos1);
    std::set<genomePosLocus>::iterator i = loci.lower_bound(locus);
    --i;
    if ( i != loci.end() ) {
      return ( ( it->ubeg1 <= upos1 ) && ( it->uend0 >= upos1 ) );
    }
    else {
      return false;
    }
  }

  bool contains1(const char* chr, int pos1) {
    return contains1(genomePos.fromPos(chr,pos1));
  }
};

#endif


