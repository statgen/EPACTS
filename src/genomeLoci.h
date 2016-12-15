#include <vector>
#include <set>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "Error.h"

class genomeLocus {
 public:
  std::string chrom;
  int beg1; // includes 1-based, excludes 0-based
  int end0; // excludes 1-based, includes 0-based
  char buf[255];

  genomeLocus(const char* c, int b, int e) : chrom(c), beg1(b), end0(e) {
    sprintf(buf,"%s:%d-%d",c,b,e);
  }

  genomeLocus(const char* region) {
    strcpy(buf,region);
    const char* pcolon = strchr(region,':');
    const char* pminus = strchr(pcolon+1,'-');
    if ( ( pcolon == NULL ) || ( pminus == NULL ) ) 
      error("Cannot parse %s in genomeLocus::genomeLocus()");
    chrom = std::string(region,0,pcolon-region);
    beg1 = atoi(pcolon+1);
    end0 = atoi(pminus+1);
  }

  const char* toString() const { 
    return buf;
  }

  bool operator< (const genomeLocus& l) const {
    if ( chrom == l.chrom ) {
      if ( beg1 == l.beg1 ) {
	return ( end0 < l.end0 );
      }
      else {
	return ( beg1 < l.beg1 );
      }
    }
    else {
      int n1 = atoi(chrom.c_str());
      int n2 = atoi(l.chrom.c_str());
      if ( ( n1 == 0 ) && ( n2 == 0 ) ) {
	return chrom < l.chrom;
      }
      else if ( ( n1 > 0 ) && ( n2 > 0 ) ) {
	return n1 < n2;
      }
      else { // treat n1 == 0 as infinite
	return ( n1 > 0 ) ? true : false;
      }
    }
  }

  unsigned long length() const { return end0-beg1+1; }

  bool overlaps (const genomeLocus& l) const {
    if ( chrom == l.chrom ) {
      if ( ( beg1 <= l.end0 )  && ( l.beg1 <= end0 ) ) {
	return true;
      }
      else {
	return false;
      }
    }
    else {
      return false;
    }
  }

  bool merge (const genomeLocus& l) {
    if ( chrom == l.chrom ) {
      if ( ( beg1-1 <= l.end0 )  && ( l.beg1-1 <= end0 ) ) {
	if ( l.beg1 < beg1 ) beg1 = l.beg1;
	if ( l.end0 > end0 ) end0 = l.end0;
	return true;
      }
      else {
	return false;
      }
    }
    else {
      return false;
    }
  }

  bool contains0(const char* chr, int pos0) const { return contains1(chr,pos0+1); }

  bool contains1(const char* chr, int pos1) const {
    if ( chrom == chr ) {
      return ( ( pos1 >= beg1 ) && ( pos1 <= end0 ) );
    }
    else {
      return false;
    }
  }
};

class genomeLoci {
 public:
  std::set<genomeLocus> loci;
  std::set<genomeLocus>::iterator it;
  bool overlapResolved;

  genomeLoci() : overlapResolved(false) {}

  void rewind() { it = loci.begin(); }
  bool next() { ++it; return ( it != loci.end() ); }
  bool isend() { return ( it == loci.end() ); }
  const genomeLocus& currentLocus() { return (*it); }
  bool empty() { return loci.empty(); }

  int numLocus() const { return (int)loci.size(); }

  bool add(const char* chr, int beg1, int end0) {
    overlapResolved = false;
    return loci.insert(genomeLocus(chr,beg1,end0)).second;
  }

  bool add(const char* region) {
    overlapResolved = false;
    return loci.insert(genomeLocus(region)).second;
  }

  int resolveOverlaps() {
    if ( !overlapResolved ) {
      std::set<genomeLocus>::iterator it;
      std::set<genomeLocus>::iterator prev;
      int numMerged = 0;
      for(it = loci.begin(); it != loci.end(); ++it) {
	if ( it != loci.begin() ) {
	  if ( prev->overlaps(*it) ) {
	    // if overlaps, erase both and insert merged one
	    genomeLocus locus = *prev;
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

  unsigned long totalLength() const {
    //resolveOverlaps();
    unsigned long sz = 0;
    std::set<genomeLocus>::iterator it;
    for(it = loci.begin(); it != loci.end(); ++it) {
      sz += it->length();
    }
    return sz;
  }

  bool moveTo(const char* chr, int pos1) {
    genomeLocus locus(chr, pos1, pos1);
    it = loci.lower_bound(locus);
    if ( it == loci.begin() ) { // do nothing
      return (it->contains1(chr,pos1));
    }
    else if ( it == loci.end() ) {
      std::set<genomeLocus>::iterator i = it;
      --i;
      if ( i->contains1(chr,pos1) ) { it = i; return true; }
      else { return false; }
    }
    else {
      if ( it->contains1(chr,pos1) ) return true;
      else {
	std::set<genomeLocus>::iterator i = it;
	--i;
	if ( i->contains1(chr,pos1) ) { it = i; return true; }
	else { return false; }
      }
    }
  }

  bool contains1(const char* chr, int pos1) {
    genomeLocus locus(chr, pos1, pos1);
    std::set<genomeLocus>::iterator i = loci.lower_bound(locus);
    --i;
    //std::set<genomeLocus>::iterator it = loci.lower_bound(locus);
    if ( i != loci.end() ) {
      //if ( ( i->chrom == chr ) && ( i->beg1 <= pos1 ) && ( i->end0 >= pos1 ) ) {
      //notice("%s %s %d %d %d %d",i->chrom.c_str(), chr, i->beg1, i->end0, pos1, ( i->chrom == chr ) && ( i->beg1 <= pos1 ) && ( i->end0 >= pos1 ) );
      //}

      return ( ( it->chrom == chr ) && ( it->beg1 <= pos1 ) && ( it->end0 >= pos1 ) );
    }
    else {
      return false;
    }
  }
};
