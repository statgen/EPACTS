#ifndef _GENOMESCORE_H_
#define _GENOMESCORE_H_

#include <map>
#include <string>
#include <cstdio>
#include "IO.h"
#include "TypeConversion.h"

class GenomeScore {
  std::string dir;
  std::map<std::string,FILE*> fpmap;
  std::string curChrom;
public:
  static int convert(const char* chrom, const char* inFile, const char* outFile) {
    FILE* fp = fopen(outFile, "wb");
    LineReader lr (inFile);
    std::vector <std::string> fd;
    int lineNo = 0;
    float value;
    while (lr.readLineBySep(&fd, "\t ")) {
      value = toFloat(fd[1]);
      if (fwrite(&value, sizeof(float), 1, fp) == 0) {
        fprintf(stderr, "Cannot write base position %s:%d in %s\n", chrom, lineNo, outFile);
        return -1;
      } else {
        ++lineNo;
      }
      if (lineNo % 1000000l == 0) {
        fprintf(stderr, "\rFinished %s:%d ...", chrom, lineNo);
      }
    }
    return lineNo;
  }
  bool openChr(const char* chrom) {
    if ( fpmap.find(chrom) != fpmap.end() ) {
      // file pointer already exist. Do nothing
      return false;
    }
    else {
      std::string fname = dir + "/chr" + chrom + ".fbin";
      FILE* fp = fopen(fname.c_str(),"rb");
      if ( fp == NULL ) {
        fprintf(stderr, "Cannot open genomeScore file %s\n",fname.c_str());
        return false;
      }
      fpmap[chrom] = fp;
      return true;
    }
  }
public:
GenomeScore(const char* _dir) : dir(_dir) {}

  bool empty() { return dir.empty(); }

  void setDir(const char* _dir) {
    if ( !dir.empty() ) {
      for(std::map<std::string,FILE*>::iterator it = fpmap.begin();
          it != fpmap.end(); ++it) {
        fclose(it->second);
      }
      fpmap.clear();
    }
    dir = _dir;
  }

  ~GenomeScore() {
    for(std::map<std::string,FILE*>::iterator it = fpmap.begin();
        it != fpmap.end(); ++it) {
      fclose(it->second);
    }
    fpmap.clear();
  }
  /**
   * @param pos is 1-based index
   */
  float baseScore(const char* chrom, int pos) {
    openChr(chrom);
    FILE* fp = fpmap[chrom];
    if ( fseek(fp, (pos-1)*4, SEEK_SET) != 0 ) {
      fprintf(stderr, "Cannot access base position %s:%d in %s\n",chrom,pos,dir.c_str());
      return 0;
    }
    float ret;
    if ( fread(&ret,sizeof(float),1,fp) == 0 ) {
      fprintf(stderr, "Cannot read base position %s:%d in %s\n",chrom,pos,dir.c_str());
      return 0;
    }
    //error("baseScore(%s,%d)=%f",chrom,pos,ret);
    return ret;
  }

  float baseScore(const char* markerID) {
    char* pcolon = strchr((char*)markerID,':');
    if ( pcolon == NULL) {
      fprintf(stderr, "Cannot parse marker %s in %s\n",markerID,dir.c_str());
      return 0;
    }
    std::string chrom(markerID,pcolon-markerID);
    int pos = atoi(pcolon+1);
    return baseScore(chrom.c_str(),pos);
  }
};

#endif /* _GENOMESCORE_H_ */
