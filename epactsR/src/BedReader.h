#ifndef _BEDREADER_H_
#define _BEDREADER_H_

#include "IO.h"
#include "TypeConversion.h"
#include "StringUtil.h"
#include <map>
#include <string>
#include <vector>
#include <algorithm>

struct Region{
  int beg; // inclusive
  int end; // exclusive
  std::string label;
};
struct RegionIndex{
  int beg;
  int end;
  std::vector<int> overlap; // store index of regions that overlaps
};

bool RegionComparator(const Region& a, const Region& b) {
  if (a.beg == b.beg)
    return (a.end < b.end);
  return (a.beg < b.beg);
};
bool RegionIndexComparator(const RegionIndex& a, const RegionIndex& b) {
  if (a.beg == b.beg)
    return (a.end < b.end);
  return (a.beg < b.beg);
};
bool inRegion(const int pos, const Region& r) {
  if ( r.beg <= pos && pos < r.end)
    return true;
  return false;
};

class BedReader {
 public:
  int open(const char* fn){
    LineReader lr(fn);
    std::string line;
    std::vector<std::string> fd;

    int totalRegion = 0;
    Region region;
    int& beg = region.beg;
    int& end = region.end;
    while (lr.readLine(&line)) {
      stringNaturalTokenize(line, "\t ", &fd);
      if (fd[0][0] == '#' || fd.size() < 3) // at least: chrom beg end [optional_label]
        continue;
      std::string& chrom = fd[0];
      if (!str2int(fd[1].c_str(), &beg)){
        fprintf(stderr, "Cannot convert [ %s ]\n", fd[1].c_str());
        continue;
      }
      if (!str2int(fd[2].c_str(), &end)){
        fprintf(stderr, "Cannot convert [ %s ]\n", fd[2].c_str());
        continue;
      }
      if (fd.size() >= 4) // default label is empty.
        region.label = fd[3];
      
      this->data[chrom].push_back(region);
      totalRegion ++;
    };

    this->createIndex();
    return totalRegion;
  };

  /**
   * @return >=1 if find; or 0 if not find.
   * @param ret_label: store all non-empty labels
   * NOTE: Here the BED boundary is defined as left-close, right-open.
   */
  int find(const char* chrom, int pos,
           std::vector<std::string>* ret_label) const {
    int nFound = 0;
    ret_label->clear();
    
    ConstIndexIter iterChrom = this->index.find(chrom);
    if (iterChrom == this->index.end()) {
      return 0;
    }
    RegionIndex temp;
    temp.beg = temp.end = pos;
    ConstRegionIndexIter iter = std::lower_bound(iterChrom->second.begin(), iterChrom->second.end(), temp, RegionIndexComparator);
    
    // iter always points to the region on the right of the region
    // e.g. find 100, but points to (200, 300) region
    // so we need to rewind to the left (subtract 1)
    if (iter == iterChrom->second.begin())
      return 0;

    iter--;
    // int idx;
    // if (iter == iterChrom->second.end()) {
    //   idx = 0;
    // } else {
    //   idx = iter - iterChrom->second.begin() - 1;
    // };
    // printf("check %d\n", iter - iterChrom->second.begin());
    const std::vector<int>& index = iter->overlap;
    ConstDataIter iterRegion = this->data.find(chrom);
    for (size_t i = 0; i < index.size(); ++i ) {
      const Region& r = (iterRegion->second)[index[i]];
      if (inRegion(pos, r)) {
        ++ nFound;
        // printf("i = %zu, pos = %d, beg = %d, end = %d\n", i, pos, r.beg, r.end);
        if (!r.label.empty()) 
          ret_label->push_back ( r.label);
      }
    }
    return nFound;
  };
  void dump() {
    puts("--------------------");
    dumpData();
    puts("--------------------");    
    dumpIndex();
  };
  void dumpData(){
    for (ConstDataIter iter = this->data.begin();
         iter != this->data.end();
         ++iter) {
      // int idx = 0;
      for (ConstRegionIter iter2 = iter->second.begin();
           iter2 != iter->second.end();
           ++iter2) {
        // fprintf(stderr, " [ %d ] %s:%d-%d -> %s \n", idx++, iter->first.c_str(), iter2->beg, iter2->end, iter2->label.c_str());   
      };
    };
  };
  void dumpIndex(){
    for (ConstIndexIter iter = this->index.begin();
         iter != this->index.end();
         ++iter) {
      int idx = 0;
      for (ConstRegionIndexIter iter2 = iter->second.begin();
           iter2 != iter->second.end();
           ++iter2) {
        std::string overlap;
        for (size_t i = 0; i < iter2->overlap.size(); ++i) {
          overlap += toStr(iter2->overlap[i]);
          overlap += ',';
        }
        fprintf(stderr, " [ %d ] %s:%d-%d -> %s \n", idx++, iter->first.c_str(), iter2->beg, iter2->end, overlap.c_str());   
      };
    };
  };

 private:
  // sort each region by region start position
  // and create index
  void createIndex() {
    DataIter iter = this->data.begin();
    for (; iter != this->data.end(); ++iter) {
      // sort
      std::vector<Region>& region = iter->second;
      std::sort(region.begin(), region.end(), RegionComparator);

      // create index
      size_t n = region.size();
      Region r = region[0];
      RegionIndex ri;
      ri.beg = r.beg;
      ri.end = r.end;
      ri.overlap.push_back(0);
      for (size_t i = 1; i < n; ++i) {
        // printf("%zu\t", i);
        if (inRegion(region[i].beg, r)) {// merge region
          r.end = std::max(r.end, region[i].end);
          ri.overlap.push_back(i);
          // printf("pushed in %zu\t", i);          
        } else { // 
          this->index[iter->first].push_back(ri);          
          r = region[i];
          ri.beg = r.beg;
          ri.end = r.end;
          ri.overlap.clear();
          ri.overlap.push_back(i);
          // printf("skip %zu\t", i);          
        };
        // printf("\n");
      }
      if (ri.overlap.size()){
        this->index[iter->first].push_back(ri);
      }
    };
    /* dumpData(); */
    /* dumpIndex(); */
  };
  static const char SEPARATOR = ',';
  std::map<std::string, std::vector< Region > > data;
  std::map<std::string, std::vector< RegionIndex> > index;

  typedef std::map<std::string, std::vector< Region > >::const_iterator ConstDataIter;
  typedef std::map<std::string, std::vector< Region > >::iterator DataIter;  
  typedef std::map<std::string, std::vector< RegionIndex> >::const_iterator ConstIndexIter;

  typedef std::vector< Region >::const_iterator ConstRegionIter;
  typedef std::vector< RegionIndex >::const_iterator ConstRegionIndexIter;
}; // BedReader


#endif /* _BEDREADER_H_ */
