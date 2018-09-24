#ifndef __TABIXED_WFILE__H
#define __TABIXED_WFILE__H

#include <zlib.h>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstdarg>

#include <htslib/bgzf.h>
#include "Error.h"

// class to read tabixed 
class wFile {
 protected:
  FILE *fp;
  gzFile gf;
  BGZF* bgf;
  int type;
  static const int MAX_LINE_SIZE = 10000000L;
  char* buf;

 public:
  ~wFile() {
    if ( buf != NULL ) delete [] buf;
    close(); // TODO : calling close() generates error. Don't know why
  }

  // read specifying one region
 wFile(const char* filename, int mode = -1)  {
   if ( mode == -1 ) { // use bgzf if ends with .gz, otherwise, plain mode
     int n = strlen(filename);
     if ( strncmp(filename+n-3,".gz",3) == 0 ) {
       type = 2;
     }
     else {
       type = 0;
     }
   }
   else {
     type = mode;
   }

   switch(type) {
   case 0:
     if ( strcmp(filename,"-") == 0 ) {
       fp = stdout;
     }
     else {
       fp = fopen(filename,"wb"); 
     }
     gf = NULL;
     bgf = NULL;
     if ( fp == NULL ) {
       fprintf(stderr,"Cannot open file %s in plain mode\n",filename);
       abort();
     }
     buf = NULL;
     break;
   case 1:
     fp = NULL;
     gf = gzopen(filename,"wb");
     bgf = NULL;
     if ( gf == NULL ) {
       fprintf(stderr,"Cannot open file %s in gzip mode\n",filename);
       abort();
     }
     buf = new char[MAX_LINE_SIZE];
     break;
   case 2:
     fp = NULL;
     gf = NULL;
     bgf = bgzf_open(filename,"wb");
     if ( bgf == NULL ) {
       fprintf(stderr,"Cannot open file %s in bgzip mode\n",filename);
       abort();
     }
     buf = new char[MAX_LINE_SIZE];
     break;
   }

   //fprintf(stderr,"*** type = %d\n",type);
  }

  int write(void* ptr, size_t count) {
    switch(type) {
    case 0:
      return (int)fwrite(ptr, 1, count, fp);
    case 1: 
      return gzwrite(gf, ptr, count);
    case 2: // bgzipped files
      return bgzf_write(bgf, ptr, count);
    default:
      fprintf(stderr,"Invalid file type %d in wFile::write()\n",type);
      abort();
    }
  }

  int putc(char c) {
    //error("putc(%c)\n",c);
    switch(type) {
    case 0:
      return fputc(c,fp);
    case 1:
      return gzputc(gf,c);
    case 2:
      return bgzf_write(bgf, (void*)&c, 1);
    default:
      fprintf(stderr,"Invalid file type %d in wFile::putc()\n",type);
      abort();
    }
  }

  int puts(const char* s) {
    int n;

    //error("puts(%s)\n",s);

    switch(type) {
    case 0:
      return fputs(s,fp);
    case 1:
      return gzputs(gf,s);
    case 2:
      n = strlen(s);
      return bgzf_write(bgf, s, n);
    default:
      fprintf(stderr,"Invalid file type %d in wFile::puts()\n",type);
      abort();
    }
  }

  void vprintf(const char* msg, va_list ap) {
    //fputs("foo()\n",stderr);
    if ( type == 0 ) {
      vfprintf(fp, msg, ap);
      //::vprintf(msg, ap);
    }
    else {
      int n = vsprintf(buf, msg, ap);
      if ( type == 1 ) {
	gzwrite(gf, buf, n);
      }
      else if ( type == 2 ) {
	bgzf_write(bgf, buf, n);
      }
      else {
	fprintf(stderr,"Invalid file type %d in wFile::printf()\n",type);
	abort();
      }
    }
  }

  void printf(const char* msg, ...) {
    va_list ap;
    va_start(ap,msg);
    vprintf(msg,ap);
    va_end(ap);
  }

  void close() {
    switch( type ) {
    case 0:
      if ( fp != NULL ) { fclose(fp); fp = NULL; }
      break;
    case 1:
      if ( gf != NULL ) { gzclose(gf); gf = NULL; } 
      break;
    case 2:
      if ( bgf != NULL ) { bgzf_close(bgf); bgf = NULL; }
      break;
    }
  }
};

#endif // __TABIXED_WFILE__H
