#ifndef _TABIXREADER_H_
#define _TABIXREADER_H_

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include "StringUtil.h"
#include "TypeConversion.h"
#include "IO.h"

#include "tabix.h"

class TabixReader{
public:
  /**
   * @param argColChrom, argColPos, argColRef, argColAlt are all 1-based column index
   */
TabixReader(const char* fn, int argColChrom, int argColPos, int argColRef, int argColAlt):
  colChrom(argColChrom-1), colPos(argColPos-1), colRef(argColRef-1), colAlt(argColAlt-1),
      missingValue("NA"), fillMissing(true) {
    // open by tabix
    tabixHandle = ti_open(fn, 0);
    if (tabixHandle == 0) {
      fprintf(stderr, "Cannot use tabix to open %s\n", fn);
      abort();
    }
    if (ti_lazy_index_load(tabixHandle) < 0) {
      fprintf(stderr, "Cannot use tabix to read index file of %s\n", fn);
      abort();
    }

    // read header
    LineReader lr(fn);
    std::string line;
    std::string prevLine;
    while (lr.readLine(&line)) {
      if (!line.empty() && line[0] != '#') {
        break;
      } else{
        prevLine = line;
      }
    };
    std::vector<std::string> h;
    if (!prevLine.empty()) {
      stringTokenize(prevLine, "\t ", &h);
      for (size_t i = 0; i < h.size(); ++i) {
        header[h[i]] = i;
      }
    }

    // validity check
    if (argColChrom < 1 || argColPos < 1 || argColRef < 1 || argColAlt < 1) {
      fprintf(stderr, "Columns are 1-based, check validity for inputs [chrom=%d,pos=%d,ref=%d,alt=%d]\n", argColChrom, argColPos, argColRef, argColAlt);
      abort();
    };
    
  };

  ~TabixReader(){
    ti_close(tabixHandle);
  };
  /**
   * @param tag: tag to be outputted 
   * @param col: 1-based index
   */
  void addTag(const char* tag, const char* headerText) {
    if (this->header.empty()) {
      fprintf(stderr, "Input file does not have header part, so [ %s ] is invalid.\n", headerText);
      return;
    }
    if (this->header.find(headerText) == this->header.end()) {
      fprintf(stderr, "Input file does not have specified header, so [ %s ] is invalid.\n", headerText);
      return;
    };
    
    this->tag.push_back(tag);
    this->col.push_back(header[headerText]);
  };
  void addTag(const std::string& tag, const std::string& headerText) {
    addTag(tag.c_str(), headerText.c_str());
  }
  /**
   * @param tag: tag to be outputted 
   * @param col: 1-based index
   */
  void addTag(const char* tag, int col) {
    if (col < 1){
      fprintf(stderr, "Column are 1-based, so [ %d ] is invalid.\n", col);
      return;
    }
    this->tag.push_back(tag);
    this->col.push_back(col-1);
  };
  void addTag(const std::string& tag, int col) {
    addTag(tag.c_str(), col);
  }
  
  /**
   * @return 0 : success
   */
  int addAnnotation(const std::string& chrom,
                    const int pos,
                    const std::string& ref,
                    const std::string& alt) {
    val.clear();
    if (tag.empty())
      return 0;
    if (snprintf(range, sizeof(range), "%s:%d-%d", chrom.c_str(), pos, pos) >= (int) sizeof(range)) {
      fprintf(stderr, "Rnage buffer is too short for %s:%d\n", chrom.c_str(), pos);
      return -1;
    }
    const char* s;
    int tid, beg, end;
    int len;
    if (ti_parse_region(tabixHandle->idx, range, &tid, &beg, &end) != 0) {
      fprintf(stderr, "Tabix cannot parse range: %s\n", range);
      return -1;
    }

    tabixIter = ti_queryi(tabixHandle, tid, beg, end);
    while ((s = ti_read(tabixHandle, tabixIter, &len)) != 0) {
      // tokenize
      stringTokenize(s, "\t ", &fd);

      // check the first match
      if ( chrom == fd[colChrom] &&
           pos == toInt(fd[colPos]) &&
           ref == fd[colRef] &&
           alt == fd[colAlt]) {
        // add result
        for (size_t i = 0; i < this->col.size(); ++i){
          if (col[i] >= (int) fd.size()) {
            val.push_back("OutOfBoundIndex");
          } else{
            val.push_back(fd[col[i]]);
          }
        }
        ti_iter_destroy(tabixIter);
        return 0;
      }
    }
    ti_iter_destroy(tabixIter);

    // fill in NA if no match found
    if (fillMissing && val.size() == 0) {
      val.resize(tag.size());
      std::fill(val.begin(), val.end(), missingValue);
    }
    return 0;
  };
  const std::vector<std::string>& getTag() const {return this->tag;};  
  const std::vector<std::string>& getAnnotation() const {return this->val;};
private:
  tabix_t* tabixHandle;
  ti_iter_t tabixIter;

  int colChrom;
  int colPos;
  int colRef;
  int colAlt;

  std::vector<std::string> fd;
  char range[128];
  std::vector<std::string> tag;
  std::vector<int> col;
  std::vector<std::string> val;

  std::map<std::string, int> header;
  
  std::string missingValue;
  bool fillMissing;             //by default, use missingValue to fill tag value
};

#endif /* _TABIXREADER_H_ */
