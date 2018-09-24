#include <map>
#include <string>
#include <cstdio>

class genomeScore {
  std::string dir;
  std::map<std::string,FILE*> fpmap;
  std::string curChrom;

  bool openChr(const char* chrom) {
    if ( fpmap.find(chrom) != fpmap.end() ) {
      // file pointer already exist. Do nothing
      return false;
    }
    else {
      std::string fname = dir + "/chr" + chrom + ".fbin";
      FILE* fp = fopen(fname.c_str(),"rb");
      if ( fp == NULL ) {
	error("Cannot open genomeScore file %s",fname.c_str());
	return false;
      }
      fpmap[chrom] = fp;
      return true;
    }
  }

 public:
   genomeScore() {}
   genomeScore(const char* _dir) : dir(_dir) {}

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

  ~genomeScore() {
    for(std::map<std::string,FILE*>::iterator it = fpmap.begin();
	it != fpmap.end(); ++it) {
      fclose(it->second);
    }
    fpmap.clear();
  }

  float baseScore(const char* chrom, int pos) {
    openChr(chrom);
    FILE* fp = fpmap[chrom];
    if ( fseek(fp, (pos-1)*4, SEEK_SET) != 0 ) {
      warning("Cannot access base position %s:%d in %s",chrom,pos,dir.c_str());
      return 0;
    }
    float ret;
    if ( fread(&ret,sizeof(float),1,fp) == 0 ) {
      warning("Cannot read base position %s:%d in %s",chrom,pos,dir.c_str());
      return 0;
    }
    //error("baseScore(%s,%d)=%f",chrom,pos,ret);
    return ret;
  }

  float baseScore(const char* markerID) {
    char* pcolon = strchr((char*)markerID,':');
    if ( pcolon == NULL) {
      warning("Cannot parse marker %s in %s",markerID,dir.c_str());
      return 0;
    }
    std::string chrom(markerID,pcolon-markerID);
    int pos = atoi(pcolon+1);
    return baseScore(chrom.c_str(),pos);
  }
};
