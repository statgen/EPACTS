#ifndef __TABIXED_FILE__H
#define __TABIXED_FILE__H

#include <zlib.h>
#include <string>
#include <cstdlib>
#include <cstdio>

#include <htslib/bgzf.h>
#include "tabix.h"

#include "Error.h"

// class to read tabixed
class pFile {
 protected:
  FILE *fp;
  gzFile gf;
  tabix_t *t;
  int tid;
  int beg;
  int end;
  int len;

  std::string fname;
  std::string chrom;
  std::string reg;
  bool head;
  ti_iter_t iter;
  int type;
  char* line;
  const ti_conf_t *idxconf;
  static const int MAX_LINE_SIZE = 10000000L;

 public:
  ~pFile() {
    close(); // TODO : calling close() generates error. Don't know why
  }

  static void tokenizeLine(const char* s, const char* delims, std::vector<std::string>& tokens) {
    const char* p = s;
    const char* c = p;
    int ndelims = strlen(delims);
    int i;
    tokens.clear();

    //fprintf(stderr,"s = '%s', strlen(s) = %d, delims = '%s'\n",s,(int)strlen(s), delims);
    while( *p != '\0' ) {
      for(i=0; i < ndelims; ++i) { if ( *p == delims[i] ) break; }
      if ( i != ndelims ) { // delimiter found
	if ( c < p )  { // unless delimiter is consencutive
	  //std::string s1(c,p-c);
	  tokens.push_back(std::string(c,p-c));
	}
	c = p+1;
      }
      ++p;
    }
    if ( c < p ) {
      tokens.push_back(std::string(c,p-c));
    }
  }

  // return 0 if not gzip, 1 if gz but not bgzf, 2 if bgzf, -1 if 10 bytes cannot be read
  static int fileType(const char *fn)
  {
    FILE* fp;
    //fprintf(stderr,"foo\n");
    uint8_t buf[10], magic[11]="\037\213\010\4\0\0\0\0\0\377";
    int n;

    if ( strcmp(fn,"-") == 0 ) { return 0; }

    if ((fp = fopen(fn, "rb")) == NULL)
    {
        return -1;
    }
    n = fread(buf, 1, 10, fp);

    //fprintf(stderr,"bar\n");

    fclose(fp);

    if ( n!=10 ) return -1;
    if ( !memcmp(magic, buf, 10) ) return 2;
    if ( !memcmp(magic, buf, 2) ) {
      return 1;
    }
    return 0;
  }

  void load(const char* filename, const char* region = NULL, bool printHeader = false) {
    head = printHeader;
    if ( region != NULL ) reg = region;
    open(filename);
    if ( reg.empty() ) {
      head = false; // head flag is valid only with specified region
      loadAll();
    }
    else if ( head ) {  // head flag is set with region specified
      loadIndex();
      loadAll();    // load the header first, and will load the region later
    }
    else {
      loadIndex();
      loadRegion();  // load the region right away
    }
  }

 pFile() : head(false), iter(NULL), type(-1), line(NULL), idxconf(NULL) {
  }

  // read specifying one region
 pFile(const char* filename, const char* region = NULL, bool printHeader = false) : head(printHeader), iter(NULL), type(-1), line(NULL), idxconf(NULL) {
    load(filename, region, printHeader);
    /*
    // open the file
    if ( region != NULL ) reg = region;
    open(filename);
    if ( reg.empty() ) {
      head = false; // head flag is valid only with specified region
      loadAll();
    }
    else if ( head ) {  // head flag is set with region specified
      loadIndex();
      loadAll();    // load the header first, and will load the region later
    }
    else {
      loadIndex();
      loadRegion();  // load the region right away
    }
    */
  }

  void updateRegion(const char* region, bool sepchr = false) {
    reg = region;
    if ( reg.empty() ) {
      loadAll();
    }
    else {
      if ( idxconf == NULL ) loadIndex();
      loadRegion(sepchr);
    }
  }

  int getLength() {
    return len;
  }

  size_t read(void* ptr, size_t count) {
    switch(type) {
    case 0:
      return fread(ptr, 1, count, fp);
    case 1:
      return gzread(gf, ptr, count);
    case 2: // bgzipped files
      return bgzf_read(t->fp, ptr, count);
    default:
      error("pFile::read() - unknown type %d\n",type);
    }
    return 0;
  }

  const char* peekLine() { return line; }

  const char* getLine() {
    //fprintf(stderr,"gerLine() called\n");

    switch(type) {
    case 0:
      if ( line == NULL ) line = new char[MAX_LINE_SIZE];
      if ( fgets(line, MAX_LINE_SIZE, fp) != NULL ) {
	//fputs(line,stderr);
	len = strlen(line); // TODO : convert to lazy evaluation
	if ( line[len-1] == '\n' ) {
	  line[len-1] = '\0';
	  --len;
	}
      }
      else {
	if ( line != NULL ) delete [] line;
	len = 0;
	line = NULL;
      }
      return line;
      /*
      size_t tn;
      //if ( line == NULL ) line = new char[MAX_LINE_SIZE];
      if ( (len = getline(&line, &tn, fp)) != -1 ) {
	if ( line[len-1] == '\n' ) {
	  line[len-1] = '\0';  // update carriage return to null character
	  --len;
	}
	return line;
      }
      else {
	//if ( line != NULL ) delete [] line;
	len = 0;
	return NULL;
      }
      */
    case 1:
      if ( line == NULL ) line = new char[MAX_LINE_SIZE];
      if ( gzgets(gf, line, MAX_LINE_SIZE) ) {
	len = strlen(line); // TODO : convert to lazy evaluation
	if ( line[len-1] == '\n' ) {
	  line[len-1] = '\0';
	  --len;
	}
      }
      else {
	if ( line != NULL ) delete [] line;
	len = 0;
	line = NULL;
      }
      return line;
    case 2:
      if ( iter == NULL ) return NULL;
      line = (char*)ti_read(t, iter, &len);
      if ( head ) { // if reached the end of the header
	if ( (int)(*line) != idxconf->meta_char ) {
	  //if ( iter != NULL ) ti_iter_destroy(iter); // close existing iterator
	  head = false;
	  loadRegion();
	  return getLine();
	}
      }
      else if ( line == NULL ) {
	//fprintf(stderr,"foo\n");
	ti_iter_destroy(iter);
	iter = NULL;
	return NULL;
      }
      return line;
    default:
      warning("Attempt to read empty file unknown file type.\n");
      return NULL;
    }
  }

  bool open(const char* filename, bool forcegz = false) {  // return true if gzipped, false otherwise
    fname = filename;
    type = fileType(filename);
    if ( forcegz ) type = 1;
    switch(type) {
    case 0:
      t = NULL;
      gf = NULL;
      if ( strcmp(filename,"-") == 0 ) {
	fp = stdin;
      }
      else {
	fp = fopen(filename,"r");
      }
      return (fp != NULL);
    case 1:
      t = NULL;
      gf = gzopen(filename,"rb");
      fp = NULL;
      return (gf != NULL);
    case 2:
      //notice("open() is called");
      if ( (t = ti_open(filename,0)) == 0 ) {
	warning("Cannot open %s with tabix..\n",filename);
	return false;
      }
      gf = NULL;
      fp = NULL;
      //notice("open() is successful");
      return true;
    default:
      warning("Cannot open %s. File is not accessible\n",filename);
      return false;
      break;
    }
    if ( !reg.empty() && type < 2 ) {
      error("File %s is not indexed, so cannot be acessed with specified region",filename);
    }
  }

  void close() {
    switch( type ) {
    case 0:
      if ( fp != NULL ) fclose(fp);
      fp = NULL;
      break;
    case 1:
      if ( gf != NULL ) gzclose(gf);
      if ( line != NULL ) delete [] line;
      gf = NULL;
      line = NULL;
      break;
    case 2:
      //notice("close() is called %d %d",iter,t->idx);
      //if ( iter != NULL ) ti_iter_destroy(iter);
      if ( t != NULL ) {
	//ti_index_destroy(t->idx);
	ti_close(t);
      }
      idxconf = NULL;
      t = NULL;
      iter = NULL;
      break;
    }
  }

  void loadAll() {
    if ( type == 2 ) {
      iter = ti_query(t,0,0,0);
    }
  }

  void loadIndex() {
    if (ti_lazy_index_load(t) < 0 ) {
      error("Failed to load the index file");
    }
    idxconf = ti_get_conf(t->idx);
  }

  void loadRegion(bool sepchr = false) {
    //notice("ti_parse_region( %s )",reg.c_str());
    if ( iter != NULL ) ti_iter_destroy(iter); // close existing iterator
    if ( ti_parse_region(t->idx, reg.c_str(), &tid, &beg, &end) != 0 ) {
      if ( sepchr ) {
	// changes all "chrAA." to "chrBB." from the files
	std::string newfname;
	int pos = 0;
	size_t ichr = 0;
	while ( (ichr = fname.find("chr",pos)) != std::string::npos ) {
	  size_t idot = fname.find_first_of("-_./",ichr);
	  std::string newchr = reg.substr(0,reg.find(':'));
	  if ( idot == std::string::npos )
	    error("Cannot find '.','_','-', or '/' after chr in the filename with --sepchr option");
	  if ( newchr.compare(0,3,"chr") == 0 )
	    newfname += (fname.substr(pos,ichr-pos) + newchr);	    
	  else
	    newfname += (fname.substr(pos,ichr-pos) + "chr" + newchr);
	  pos = idot;
	}
	newfname += fname.substr(pos);
	fname = newfname;

	notice("Changing the VCF file name to %s",fname.c_str());

	/*
	//notice("loadRegion(true) %s",reg.c_str());
	// assume that current filename is [prefix]chr[chr].[suffix]
	int ichr = fname.find("chr");
	int idot = fname.find('.',ichr);
	std::string newchr = reg.substr(0,reg.find(':'));
	std::string prefix = fname.substr(0,ichr);
	std::string suffix = fname.substr(idot);
	fname = prefix + "chr" + newchr + suffix;
	//notice("open(%s)",fname.c_str());
	*/

	if ( fileType(fname.c_str()) < 0 ) {
	  warning("Cannot parse region %s.. Returning empty",reg.c_str());
	  iter = NULL;
	}
	else {
	  close();
	  open(fname.c_str());
	  loadIndex();
	  loadRegion();
	}
      }
      else {
	warning("Cannot parse region %s.. Returning empty",reg.c_str());
	iter = NULL;
      }
    }
    else {
      //notice("ti_query(%x, %d, %d, %d)",t,tid,beg,end);
      iter = ti_queryi(t,tid,beg,end);
    }
  }
};

#endif // __TABIXED_FILE
