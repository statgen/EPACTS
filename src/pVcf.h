#ifndef __TABIXED_PVCF_H
#define __TABIXED_PVCF_H

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <set>
#include "pFile.h"
#include "Error.h"
#include "boolParser.h"

class pVcf {
public:
  // member variables
  int nInds;           // number of individuals subselected
  int nMarkers;        // number of markers
  std::string key;     // key string
  boolParser parser;   // per-marker pattern matching parser
  bool passOnly;       // filtering option
  pFile tf;            // file handle
  std::vector<std::string> inds;     // individual IDs
  std::vector<std::string> markers;  // marker IDs
  std::vector<std::string> chroms;   // chromosome names
  std::vector<int> pos1s;            // marker positions
  std::vector<std::string> refs;     // reference alleles
  std::vector<std::string> alts;     // non-reference alleles
  std::vector<float> genos;          // GT : genotypes
  std::vector<uint8_t> PLs;          // GL or PL : genotype likelihoods
  std::vector<int> depths;           // GD or DP : genotype depth
  std::vector<int> numAlleles;       // AN 
  std::vector<float> sumAlleles;     // AC
  std::vector<float> sumsqAlleles;   // squared sum of AC
  std::vector<int> icols;            // individual indices to subset
  std::vector<float> AFs;            // allele frequency estimates

  std::set<std::string> markerSet;
  std::string groupID;

  bool glFlag;
  bool dpFlag;
  bool phredFlag;
  std::vector<double> p2e;

  static const float NAN_FLT;
  static const double NAN_DBL;

  static int parseMarkerID(const char* markerID, std::string& chrom, std::string& beg, std::string& end, std::string& name, std::string& anno) {
    const char* pp = markerID;
    const char* pn = markerID;
    anno.clear();
    int step = 0; // 0 : CHROM[:], 1 : BEG [-_], 2 : END [_], 3 : NAME [:]
    while( *pn != '\0' ) {
      switch(step) {
      case 0: // CHROM[:]
	if ( *pn == ':' ) { // END PARSING
	  chrom.assign( pp, pn - pp ); // copy chrom
	  step = 1;
	  ++pn;
	  pp = pn;
	}
	// otherwise, just skip else {}
	break;
      case 1: // BEG
	if ( *pn == '-' ) { // INTERVAL is given
	  beg.assign( pp, pn - pp );
	  step = 2;
	  ++pn;
	  pp = pn;
	}
	else if ( *pn == '_' ) { // BEG==END
	  beg.assign( pp, pn - pp );
	  end = bed;
	  step = 3;
	  ++pn;
	  pp = pn;
	}
	break;
      case 2: // END
	if ( *pn == '_' ) { // BEG==END
	  end.assign( pp, pn - pp );
	  end = bed;
	  step = 3;
	  ++pn;
	  pp = pn;
	}
	break;
      case 3: // NAME
	if ( pn == ':' ) {
	  name.assign( pp, pn - pp );
	  step = 4;
	  ++pn;
	  pp = pn;
	}
	break;
      }
    }

    if ( pn > pp ) {
      switch(step) {
      case 0:
	chrom.assign( pp, pn - pp );
	break;
      case 1:
	beg.assign( pp, pn - pp );
	end = bed;
	step = 2;
	break;
      case 2:
	end.assign( pp, pn - pp );
	break;
      case 3:
	name.assign( pp, pn - pp );
	break;
      case 4:
	anno.assign( pp, pn - pp );
	break;
      }
    }
    return step;
  }

  void load(const char* vcf, const char* region, const char* _key, const char* rule, bool pass, std::set<std::string>& idset) {
    nInds = 0;
    nMarkers = 0;
    key = _key;
    parser = rule;
    passOnly = pass;
    tf.load(vcf, region, true);

    if ( key == "GL" ) {
      glFlag = true;
      phredFlag = false;
      initPhred2Err();
    }
    else if ( key == "PL" ) {
      glFlag = phredFlag = true;
      initPhred2Err();
    }
    else {
      glFlag = phredFlag = false;
    }
    
    // read VCF header
    char* line;
    int nc;
    while ( (line = (char*)tf.getLine()) != NULL ) {
      nc = tf.getLength();
      if ( line[0] == '#' ) {
	if ( strncmp(line,"#CHROM",6) == 0 ) {
	  nInds = parseInds(line, idset);
	  break;
	}
	else if ( line[1] == '#' ) {
	  // parse meta line???
	}
      }
      else {
	error("Non-header line %s is observed before #CHROM in %s",line,vcf);
      }
    }
    
    if ( idset.size() != icols.size() ) {
      warning("Identified %d individuals from index file, and only %d overlaps with VCF",(int)idset.size(),(int)icols.size());
    }
  }

  void load(const char* vcf, const char* region, const char* _key, const char* rule, bool pass, std::vector<int>& subcols) {
    nInds = 0;
    nMarkers = 0;
    key = _key;
    parser = rule;
    passOnly = pass;
    tf.load(vcf, region, true);

    if ( key == "GL" ) {
      glFlag = true;
      phredFlag = false;
      initPhred2Err();
    }
    else if ( key == "PL" ) {
      glFlag = phredFlag = true;
      initPhred2Err();
    }
    else {
      glFlag = phredFlag = false;
    }
    
    // read VCF header
    char* line;
    int nc;
    while ( (line = (char*)tf.getLine()) != NULL ) {
      nc = tf.getLength();
      if ( line[0] == '#' ) {
	if ( strncmp(line,"#CHROM",6) == 0 ) {
	  nInds = parseInds(line, subcols);
	  break;
	}
	else if ( line[1] == '#' ) {
	  // parse meta line???
	}
      }
      else {
	error("Non-header line %s is observed before #CHROM in %s",line,vcf);
      }
    }
    
    if ( ( subcols.size() > 0 ) && ( (int)subcols.size() != nInds ) ) {
      warning("Identified %d individuals from index file, and nInds = %d",(int)subcols.size(),nInds);
    }
  }

 pVcf() : nInds(0), nMarkers(0), passOnly(false) {}

 pVcf(const char* vcf, const char* region, const char* _key, const char* rule, bool pass, std::vector<std::string> oids) :
  nInds(0), nMarkers(0), key(_key), parser(rule), passOnly(pass), tf(vcf, region, true) {
      if ( key == "GL" ) {
	glFlag = true;
	phredFlag = false;
	initPhred2Err();
      }
      else if ( key == "PL" ) {
	glFlag = phredFlag = true;
	initPhred2Err();
      }
      else {
	glFlag = phredFlag = false;
      }

      // read VCF header
      char* line;
      int nc;
      while ( (line = (char*)tf.getLine()) != NULL ) {
	nc = tf.getLength();
	if ( line[0] == '#' ) {
	  if ( strncmp(line,"#CHROM",6) == 0 ) {
	    nInds = parseInds(line, oids);
	    break;
	  }
	  else if ( line[1] == '#' ) {
	    // parse meta line???
	  }
	}
	else {
	  error("Non-header line %s is observed before #CHROM in %s",line,vcf);
	}
      }

      if ( oids.size() != icols.size() ) {
	warning("Identified %d individuals from index file, and %d overlaps with VCF",(int)oids.size(),(int)icols.size());
      }
    }

 pVcf(const char* vcf, const char* region, const char* _key, const char* rule, bool pass, std::set<std::string>& idset) :
  nInds(0), nMarkers(0), key(_key), parser(rule), passOnly(pass), tf(vcf, region, true) {
      if ( key == "GL" ) {
	glFlag = true;
	phredFlag = false;
	initPhred2Err();
      }
      else if ( key == "PL" ) {
	glFlag = phredFlag = true;
	initPhred2Err();
      }
      else {
	glFlag = phredFlag = false;
      }

      // read VCF header
      char* line;
      int nc;
      while ( (line = (char*)tf.getLine()) != NULL ) {
	nc = tf.getLength();
	if ( line[0] == '#' ) {
	  if ( strncmp(line,"#CHROM",6) == 0 ) {
	    nInds = parseInds(line, idset);
	    break;
	  }
	  else if ( line[1] == '#' ) {
	    // parse meta line???
	  }
	}
	else {
	  error("Non-header line %s is observed before #CHROM in %s",line,vcf);
	}
      }

      if ( idset.size() != icols.size() ) {
	warning("Identified %d individuals from index file, and %d overlaps with VCF",(int)idset.size(),(int)icols.size());
      }
  }

  // always parse the header
  pVcf(const char* vcf, const char* region, const char* _key, const char* rule, bool pass, const char* indf) :
  nInds(0), nMarkers(0), key(_key), parser(rule), passOnly(pass), tf(vcf, region, true) {
      if ( key == "GL" ) {
	glFlag = true;
	phredFlag = false;
	initPhred2Err();
      }
      else if ( key == "PL" ) {
	glFlag = phredFlag = true;
	initPhred2Err();
      }
      else {
	glFlag = phredFlag = false;
      }
      
      // read indf files if needed
      std::set<std::string> idset;
      char* line = NULL;
      int nc;
      if ( indf != NULL ) {
	pFile tind(indf);
	while( (line = (char*)tind.getLine()) != NULL ) {
	  if ( line[0] == '#' ) continue;
	  char* p = line;
	  while( ( *p != ' ' ) && ( *p & 0xf0 ) ) { ++p; }
	  std::string id(line,p-line);
	  if ( idset.find(id) != idset.end() ) {
	    error("ERROR: Duplicate individual ID %s in %s",id.c_str(),indf);
	  }
	  idset.insert(id);
	}
      }

      // read VCF header
      while ( (line = (char*)tf.getLine()) != NULL ) {
	nc = tf.getLength();
	if ( line[0] == '#' ) {
	  if ( strncmp(line,"#CHROM",6) == 0 ) {
	    nInds = parseInds(line, idset);
	    break;
	  }
	  else if ( line[1] == '#' ) {
	    // parse meta line???
	  }
	}
	else {
	  error("Non-header line %s is observed before #CHROM in %s",line,vcf);
	}
      }

      if ( ( idset.size() > 0 ) && ( idset.size() != icols.size() ) ) {
	warning("Identified %d individuals from index file, and %d overlaps with VCF",(int)idset.size(),(int)icols.size());
      }
  }

  pVcf(const char* filename, const char* region, const char* _key, const char* rule, bool pass, bool header, int* pIcols, int nIcols) :
  nInds(0), nMarkers(0), key(_key), parser(rule), passOnly(pass), tf(filename, region, header) {
    if ( key == "GL" ) {
      glFlag = true;
      phredFlag = false;
      initPhred2Err();
    }
    else if ( key == "PL" ) {
      glFlag = phredFlag = true;
      initPhred2Err();
    }
    else {
      glFlag = phredFlag = false;
    }

    // assume individual index is in increasing order
    if ( pIcols != NULL ) {
      //std::cout << nIcols << std::endl;
      for(int i=0; i < nIcols; ++i) {
	icols.push_back(pIcols[i]);  // assume 0-based index
	if ( ( i > 0 ) && ( icols[i] <= icols[i-1] ) ) {
	  error("Column index (ind file) must be ordered");
	}
      }
      nInds = nIcols;
    }
    //notice("pVcf::pVcf(..) ended");
  }

 pVcf(const char* filename, const char* region, const char* _key, bool header, int* pIcols, int nIcols, const char** pMarkers, int nMarkers) :
  nInds(0), nMarkers(0), key(_key), tf(filename, region, header) {    
    for(int i=0; i < nMarkers; ++i) {
      markerSet.insert(std::string(pMarkers[i]));
      //fprintf(stderr,"%s\n",pMarkers[i]);
    }

    if ( key == "GL" ) {
      glFlag = true;
      phredFlag = false;
      initPhred2Err();
    }
    else if ( key == "PL" ) {
      glFlag = phredFlag = true;
      initPhred2Err();
    }
    else {
      glFlag = phredFlag = false;
    }

    // assume individual index is in increasing order
    if ( pIcols != NULL ) {
      //std::cout << nIcols << std::endl;
      for(int i=0; i < nIcols; ++i) {
	icols.push_back(pIcols[i]);  // assume 0-based index
	if ( ( i > 0 ) && ( icols[i] <= icols[i-1] ) ) {
	  error("Column index (ind file) must be ordered");
	}
      }
      nInds = nIcols;
    }
  }

  void initPhred2Err() {
    for(int i=0; i < 256; ++i) {
      p2e.push_back(pow(10,-0.1*i));
    }
  }

  void clear() {
    markers.clear();
    chroms.clear();
    pos1s.clear();
    refs.clear();
    alts.clear();
    genos.clear();
    PLs.clear();
    numAlleles.clear();
    sumAlleles.clear();
    sumsqAlleles.clear();
    AFs.clear();
    markerSet.clear();
    nMarkers = 0;
  }

  float callRate(int m) {
    return ((float)numAlleles[m] / (float)nInds / 2.);
  }

  float alleleFreq(int m) {
    if ( glFlag ) {
      if ( (int)AFs.size() < m + 1 ) {
	AFs.resize(m+1,-1);
      }
      if ( AFs[m] < 0 ) {
	AFs[m] = (float)(emAF(m,1e-6).first);
      }
      return AFs[m];
    }
    else {
      return ( numAlleles[m] > 0 ) ? (sumAlleles[m] / (float)numAlleles[m]) : 0;
    }
  }

  int MAC(int m) {
    //if ( glFlag ) {
      //return (int)floor(MAF(m) * nInds + .5);
    //}
    //else {
    return (int)floor(((sumAlleles[m] < numAlleles[m]-sumAlleles[m]) ? sumAlleles[m] : numAlleles[m]-sumAlleles[m])+.5);
      //}
  }

  // works only for discrete genotypes
  // numAlleles[m] = r + h + a
  // sumAlleles[m] = h + 2a
  // sumsqAlleles[m] = h + 4a
  // a = ((h+4a)-(h+2a))/2
  // r = (r+h+a) - (h+2a) + ((h+4a)-(h+2a))/2
  //   = (r+h+a) - 1.5(h+2a) + 0.5(h+4a)
  int HOMMINC(int m) {
    return (int)floor( ((sumAlleles[m] < numAlleles[m]-sumAlleles[m]) ? (sumsqAlleles[m]-sumAlleles[m])/2. : ( numAlleles[m] - 1.5 * sumAlleles[m] + 0.5 * sumsqAlleles[m] )) + .5);
  }

  // works only 
  // numAlleles[m] = 2*(r + h + a)
  // sumAlleles[m] = h + 2a
  // sumsqAlleles[m] = h + 4a
  // a = ((h+4a)-(h+2a))/2
  // r = (r+h+a) - (h+2a) + ((h+4a)-(h+2a))/2
  //   = (r+h+a) - 1.5(h+2a) + 0.5(h+4a)
  // h = (r+h+a) - a - r
  void GENOCNT(int m, float* cnts) {
    cnts[0] = ( numAlleles[m]/2 - 1.5 * sumAlleles[m] + 0.5 * sumsqAlleles[m] );
    cnts[2] = (sumsqAlleles[m]-sumAlleles[m])/2.;
    cnts[1] = numAlleles[m]/2 - cnts[0] - cnts[2];
  }

  float MAF(int m) {
    float f = alleleFreq(m);
    return (f > .5 ? 1.-f : f);
  }

  float RSQ(int m) {
    int n = numAlleles[m]/2;
    if ( n > 1 ) {
      float s = sumAlleles[m];
      float varObs = (sumsqAlleles[m]-s*s/n)/(n-1);     // var(x)
      float varExp = 2. * s * ( n + n - s ) / ( n * n ); // 2pq
      return (varObs > varExp) ? 1.0 : varObs/varExp;
    }
    else {
      return 0;
    }
  }

  void print(FILE* fp) {
    float v;
    //fprintf(stderr,"%d %d %s\n",nInds,nMarkers,inds[nInds-1].c_str());
    if ( inds.size() > 0 ) {
      fprintf(fp, "#MARKER");
      for(int i=0; i < nInds; ++i) {
	if ( inds.empty() ) {
	  fprintf(fp, "\tID%d", i+1);
	}
	else {
	  fprintf(fp, "\t%s", inds[i].c_str());
	}
      }
      fprintf(fp,"\n");
    }
    for(int i=0; i < nMarkers; ++i) {
      fprintf(fp, "%s", markers[i].c_str());
      for(int j=0; j < nInds; ++j) {
	v = genos[i*nInds + j];
	if ( isnan(v) ) {
	  fprintf(fp,"\tNA");
	}
	else {
	  fprintf(fp, "\t%.4f", genos[i*nInds + j]);
	}
      }
      fprintf(fp,"\n");
    }
  }

  int parseMarkers(char* line, int startIdx = 9) {
    //notice("pVcf::parseMarkers() called");
    //, std::vector<std::string>& m, std::vector<float>& v, int startIdx = 9, const char* key = "GT", const char** infoSubstrs = NULL, int nSubstr = 0) {
    char* p;
    char* n;
    char* pch = line;
    char* nch = NULL;
    std::string chrom;
    std::string pos;
    std::string ref;
    std::string alt;
    std::string markerId;
    std::string markerKey;
    std::string s;
    int keyIdx = 0;
    int i, j, k;
    int lKey = (int)key.size();
    int AN = 0;
    float AC = 0, sqAC = 0;
    float f;

    //int* ls = (nSubstr > 0 ) ? new int[nSubstr] : NULL;
    //int* js = (nSubstr > 0 ) ? new int[nSubstr] : NULL;
    //for(int i=0; i < nSubstr; ++i) {
    //  ls[i] = strlen(infoSubstrs[i]);
      //fprintf(stderr,"%s\t%d\n",infoSubstrs[i],ls[i]);
    //}
    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i < startIdx ) {
	if ( i < 7 ) {
	  if ( nch == NULL ) s.assign( pch );
	  else s.assign( pch, nch - pch );

	  switch(i) {
	  case 0:
	    chrom = s; break;
	  case 1:
	    pos = s; break;
	  case 2:
	    markerId = s; break;
	  case 3:
	    ref = s; break;
	  case 4:
	    alt = s; 
	    if ( ! markerSet.empty() ) {
	      markerKey = chrom+":"+pos+"_"+ref+"/"+alt;
	      if ( markerSet.find(markerKey) == markerSet.end() ) {
		return -1;
	      }
	    }
	    break;
	  case 6:
	    if ( ( passOnly ) && ( s != "PASS" ) ) return -1;
	  }
	}
	else if ( i == 7 ) {
	  if ( nch == NULL ) { if ( !parser.parse( pch ) ) return -1; }
	  else { if ( !parser.parse( pch, nch - pch ) ) return -1; }
	}
	else if ( i == 8 ) { // parse FORMAT field
	  if ( pch[0] == key[0] && pch[1] == key[1] ) { // comparing two characters are not exactly right
	    keyIdx = 0;
	  }
	  else if ( nch == NULL ) {
	    error("VCF has FORMAT field but does not have any genotype");
	  }
	  else {
	    k = 0;
	    keyIdx = 0;
	    p = pch;
	    while( p < nch ) {
	      if ( *p == ':' ) {
		if ( k >= lKey ) {
		  break;
		}
		else {
		  ++keyIdx;
		  k = 0;
		}
	      }
	      else {
		if ( ( k == 2 ) || ( key[k] == *p ) ) {
		  //if ( key[k] == *p ) {
		  ++k;
		}
		else {
		  k = 0;
		}
	      }
	      ++p;
	    }
	    if ( ( p == nch ) && ( k != lKey ) ) {
	      warning("Cannot find %s in the FORMAT field at marker %s:%s_%s/%s .. Skipping",key.c_str(),chrom.c_str(),pos.c_str(),ref.c_str(),alt.c_str());
	      return -1;
	    }
	  }
	}
      }
      else {
	if ( icols.empty() || ( ( j < (int)icols.size() ) && ( icols[j] == i - startIdx ) ) ) {
	  p = pch;
	  n = NULL;
	  
	  // reach to the key index
	  if ( keyIdx > 0 ) {
	    for(int i=0; (i < keyIdx) && (p != NULL); ++i) {
	      n = strchr(p, ':');
	      p = (n == NULL) ? NULL : n+1;
	    }
	  }
	  
	  // if the field contains '/' or '|', add two values
	  // if the field is '.', return NA
	  // otherwise, convert the field into float and return
	  if ( ( p == NULL ) || ( p[0] == '.' ) ) { // missing
	    if ( glFlag ) {
	      PLs.push_back(0);
	      PLs.push_back(0);
	      PLs.push_back(0);
	    }
	    else {
	      genos.push_back(NAN_FLT);
	    }
	  }
	  else if ( ( p[1] == '/' ) || ( p[1] == '|' ) ) {
	    f = (float)((p[0] - '0') + ( p[2] - '0' )); // ignore phase info
	    genos.push_back(f);
	    AN += 2;
	    AC += f;
	    sqAC += f*f;
	  }
	  else {
	    if ( glFlag ) { // search for three commas
	      char* c1 = strchr(p,',');
	      if ( c1 == NULL ) error("Cannot parse %s field (currently suppports only biallelic autosomal variants",key.c_str());
	      char* c2 = strchr(c1+1,',');
	      if ( c2 == NULL ) error("Cannot parse %s field (currently suppports only biallelic autosomal variants",key.c_str());
	      if ( phredFlag ) {
		int pl = atoi(p);
		if ( pl > 255 ) pl = 255;
		PLs.push_back(pl);
		pl = atoi(c1+1);
		if ( pl > 255 ) pl = 255;
		PLs.push_back(pl);
		pl = atoi(c2+1);
		if ( pl > 255 ) pl = 255;
		PLs.push_back(pl);
	      }
	      else {
		float gl = strtof(p,NULL);
		if ( gl < -25.5 ) PLs.push_back(255);
		else PLs.push_back((int)(-10*gl+.5));

		gl = strtof(c1+1,NULL);
		if ( gl < -25.5 ) PLs.push_back(255);
		else PLs.push_back((int)(-10*gl+.5));

		gl = strtof(c2+1,NULL);
		if ( gl < -25.5 ) PLs.push_back(255);
		else PLs.push_back((int)(-10*gl+.5));
	      }
	    }
	    else {
	      f = strtof(p,NULL);
	      genos.push_back(f);
	      AN += 2;
	      AC += f;
	      sqAC += f*f;
	    }
	  }
	  ++j;
	}
	else {
	  //std::cout << "Skipping " << i << "\t" << j << std::endl;
	}
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }

    if ( (nInds == 0) && icols.empty() && ((int)genos.size() == j) ) nInds = j;

    if ( nInds != j ) {
      fprintf(stderr,"i=%d, j=%d, nInds=%d, icols.size()=%d, genos.back()=%d\n",i,j,nInds,(int)icols.size(),(int)genos.back());
      abort();
    }

    // if GL or PL flag is set, automatically put dosage as quantitative genotypes
    if ( glFlag ) {
      //notice("pVcf::parseMarkers() - glFlag is on");

      // set allele frequency
      int m = markers.size();  //notice("pVcf::parseMarkers() - m = %d",m);
      float af = (float)(emAF(m,1e-6).first); //notice("pVcf::parseMarkers() - af = %f",af);
      if ( (int)AFs.size() > m ) { AFs[m] = af; }
      else if ( (int)AFs.size() == m ) { AFs.push_back(af); }
      else { error("pVcf::parseMarkers() -- AFs.size() < m"); }

      double p0,p1,p2;
      int kos;
      for(int k=0; k < j; ++k) {
	// calculate genotype dosages under HWE
	// Pr(G|D) = Pr(D|G)Pr(G)
	kos = 3*(m*nInds+k);
	p0 = p2e[PLs[kos+0]] * (1.-AFs[m]) * (1.-AFs[m]);
	p1 = p2e[PLs[kos+1]] * 2. * AFs[m] * (1.-AFs[m]);
	p2 = p2e[PLs[kos+2]] * AFs[m] * AFs[m];
	if ( p0+p1+p2 > 0 ) {
	  genos.push_back((p1+2*p2)/(p0+p1+p2));
	}
	else {
	  genos.push_back(AFs[m]*2.);
	}
	AN += 2.;
	AC += genos.back();
	sqAC += genos.back()*genos.back();
      }
    }

    markers.push_back(chrom+":"+pos+"_"+ref+"/"+alt);
    chroms.push_back(chrom);
    pos1s.push_back(atoi(pos.c_str()));
    refs.push_back(ref);
    alts.push_back(alt);
    numAlleles.push_back(AN);
    sumAlleles.push_back(AC);
    sumsqAlleles.push_back(sqAC);

    //fprintf(stderr,"%s\n",markers.back().c_str());

    return j;
  }

  int parseInds(char* line, std::set<std::string>& idset, int startIdx = 9) {
    char* pch = line;
    char* nch = NULL;
    int i, j;

    icols.clear();
    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i >= startIdx ) {
	std::string id = (nch == NULL) ? std::string(pch) : std::string(pch,nch-pch);
	if ( idset.empty() || (idset.find(id) != idset.end()) ) {
	  icols.push_back(i - startIdx);
	  inds.push_back(id);
	  ++j;
	}
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }
    return j;
  }

  int parseInds(char* line, std::vector<std::string>& subids, int startIdx = 9) {
    char* pch = line;
    char* nch = NULL;
    int i, j;

    icols.clear();
    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i >= startIdx ) {
	std::string id = (nch == NULL) ? std::string(pch) : std::string(pch,nch-pch);
	if ( subids[(int)icols.size()] == id ) {
	  icols.push_back(i - startIdx);
	  inds.push_back(id);
	  ++j;
	}
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }
    if ( j != (int)subids.size() ) {
      error("ERROR in pVcf::parseInds() - not all subids matches");
    }
    return j;
  }

  int parseInds(char* line, std::vector<int>& subcols, int startIdx = 9) {
    //notice("pVcf::parseInds(char*, std::vector<int>&) called");

    char* pch = line;
    char* nch = NULL;
    int i, j;

    icols.clear();
    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i >= startIdx ) {
	if ( subcols.empty() || subcols[j] == i - startIdx ) {
	  icols.push_back(i - startIdx);
	  std::string id = (nch == NULL) ? std::string(pch) : std::string(pch,nch-pch);
	  inds.push_back(id);
	  ++j;
	}
	else { // Skip the individual
	  //std::cerr << "Skipping " << i << "\t" << j << std::endl;
	  //abort();
	}
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }
    return j;
  }

  // read markers until reach the end of file
  int readMarkerGroup(const char** pMarkers, int nMarkers, int startIndex = 0, bool sepchr = false) {
    std::vector<std::string> markerIDs;
    for(int i=0; i < nMarkers; ++i) {
      markerIDs.push_back(pMarkers[i]);
    }
    return readMarkerGroup(markerIDs,startIndex,sepchr);
  }

  int readMarkerGroup(std::vector<std::string>& markerIDs, 
		      int startIndex = 0, bool sepchr = false) 
  {
    int nc;
    char* line = NULL;

    std::string curChrom;
    int beg = 1e9, end = 0;

    std::vector<std::string> tokens;
    for(int i=startIndex; i < (int)markerIDs.size(); ++i) {
      pFile::tokenizeLine(markerIDs[i].c_str(),":_/",tokens);
      if ( tokens.size() != 4 ) {
	error("Cannot parse markerID %s",markerIDs[i].c_str());
      }
      std::string region(tokens[0]);
      region += ":";
      region += tokens[1];
      region += "-";
      region += tokens[1];

      if ( i == startIndex ) {
	curChrom = tokens[0];
	beg = end = atoi(tokens[1].c_str());
      }
      else {
	if ( curChrom != tokens[0] ) { 
	  curChrom = "multichrs";
	  beg = end = 0;
	  //error("Currently group across different chromosomes are not supported"); 
	}
	else {
	  int bp = atoi(tokens[1].c_str());
	  if ( beg > bp ) { beg = bp; }
	  if ( end < bp ) { end = bp; }
	}
      }

      tf.updateRegion(region.c_str(),sepchr);

      markerSet.clear();
      markerSet.insert(markerIDs[i]);
      while ( (line = (char*)tf.getLine()) != NULL ) {
	nc = tf.getLength();
	int cols2 = parseMarkers(line); 
	if ( cols2 >= 0 ) {
	  if ( nInds != cols2 ) {
	    error("pVcf::readMarkers() : Column size does not match : %d vs %d at marker %s", nInds, cols2, markers[markers.size()-1].c_str() );
	    abort();
	  }
	  ++nMarkers;
	}
      }
    }

    if ( startIndex > 0 ) {
      char tmp[1024];
      sprintf(tmp,"%s:%d-%d_%s",curChrom.c_str(),beg,end,markerIDs[0].c_str());
      //markerIDs[0] = tmp;
      groupID = tmp;
    }
    return nMarkers;
  }

  // read markers until reach the end of file
  int readMarkers(int m = 0, bool del = true) {
    int nc;
    char* line = NULL;
    
    if ( del ) clear();
    
    int mStart = nMarkers;

    //notice("pVcf::readMarkers(%d) called",m);

    while ( (line = (char*)tf.getLine()) != NULL ) {
      //notice("bar");
      nc = tf.getLength();
      //notice("nMarkers=%d, line[0] = %c, icols.size() = %d, nc=%d",nMarkers,line[0],(int)icols.size(), nc);
      if ( line[0] == '#' ) {
	if ( strncmp(line,"#CHROM",6) == 0 ) {
	  nInds = parseInds(line, icols);
	}
      }
      else {
	int cols2 = parseMarkers(line); 
	//notice("cols2 = %d",cols2);
	if ( cols2 >= 0 ) {
	  if ( nInds != cols2 ) {
	    error("pVcf::readMarkers() : Column size does not match : %d vs %d at marker %s", nInds, cols2, markers[markers.size()-1].c_str() );
	    abort();
	  }
	  ++nMarkers;
	  if ( ( m > 0 ) && ( m <= nMarkers - mStart ) ) break;
	}
      }
    }
    //notice("Returning %d",nMarkers-mStart);
    return ( nMarkers-mStart );
  };

  void close() {
    tf.close();
  }

  /*
  static FILE* openVCF(const char* filename, const char* region = NULL, const char* tabix = "/usr/cluster/bin/tabix" ) {
    // check if the file exists
    if ( exists(filename) == 0 ) {
      std::cerr << "Cannot open file " << filename << " for reading" << std::endl;
      return NULL;
    }
    
    std::string fn(filename);
    FILE* fp;
    if ( fn.substr(fn.size()-3) == ".gz" ) {
      std::string cmd = tabix;
      if ( region == NULL ) {
	cmd = "zcat ";
      }
      else if ( exists(tabix) == 0 ) {
	std::cerr << "Cannot find tabix binary " << tabix << ". Failed opening with region specified";
	return NULL;
      }
      else {
	cmd += " -h ";
      }
      cmd += fn;
      if ( region != NULL ) {
	cmd += " ";
	cmd += region;
      }
      //std::cout << cmd << std::endl;
      fp = popen(cmd.c_str(),"r");
      return fp;
    }
    else {
      fp = fopen(filename,"r");
      return fp;
    }
  }
  */

  /**
     Computes HWE allele frequencies using EM algorithm
     Input:
     1)Genotype Likelihoods (qscores)
     
     Output: 
     1) estimated HWE allele frequencies
     2) sample size
  */
  std::pair<double,double> emAF(int m, double eps) {
    // initial step : start from randomly assigned allele frequencies
    std::vector<int> indices;
    for(int i=0; i < nInds; ++i) {
      indices.push_back(i);
    }
    return emAF(m, indices,eps);
  }

  std::pair<double,double> emAF(int m, std::vector<int>& indices, double eps) {
    // initialization : pick an arbitrary AF from 
    //notice("pVcf::emAF(%d, %d, %lf), p2e.size() = %d",m,(int)indices.size(),eps,(int)p2e.size());

    std::vector<int>::iterator it;
    int i, c, n, r;
    double p = .5 + rand()/(RAND_MAX+1.)*.3; // random AF
    double q;
    double f0,f1,f2, fsum, sum;
    double llk = 0;
    std::vector<double> post;
    n = (int)indices.size();
    post.resize(n*3);

    for(r = 0; r < 100; ++r) {
      sum = 0;
      c = 0;
      q = 1.-p;
      for(it = indices.begin(); it != indices.end(); ++it) {
	i = 3*m*nInds + (*it)*3;
	f0 = q * q * p2e[PLs[i]];
	f1 = 2. * p * q * p2e[PLs[i+1]];
	f2 = p * p * p2e[PLs[i+2]];
	fsum = f0+f1+f2;
	post[c++] = f0/fsum;
	sum += (post[c++] = f1/fsum);
	sum += (2 * (post[c++] = f2/fsum));
      }
      p = sum / (2*n);
      if ( fabs(p + q - 1.) < eps ) break;
    }

    //notice("pVcf::emAF() - p = %lf",p);

    // Pr(Data|AF) = \sum_g Pr(Data|g)Pr(g|AF)
    q = 1.-p;
    for(it = indices.begin(); it != indices.end(); ++it) {
      i = 3*m*nInds + (*it)*3;
      f0 = q * q * p2e[PLs[i]];
      f1 = 2 * p * q * p2e[PLs[i+1]];
      f2 = p * p * p2e[PLs[i+2]];
      fsum = f0+f1+f2;
      llk += log(fsum);
    }
    //notice("pVcf::emAF() finished - returning (%lf, %lf)",p,llk);
    return std::pair<double,double>(p,llk);
  }

  double LRT(int m, std::vector<int>& g1, std::vector<int>& g2, double& af1, double& af2) {
    std::vector<int> g;
    std::pair<double,double> p;
    for(int i=0; i < (int)g1.size(); ++i) { g.push_back(g1[i]); }
    for(int i=0; i < (int)g2.size(); ++i) { g.push_back(g2[i]); }
    double llk0 = emAF(m, g, 1e-6).second;
    p = emAF(m, g1, 1e-6);
    af1 = p.first;
    double llk1 = p.second;
    p = emAF(m, g2, 1e-6);
    af2 = p.first;
    llk1 += p.second;

    //printf("%d\t%lf\t%lf\t%lf\n",m,llk0,llk1,llk1-llk0);
    
    return( 2*(llk1 - llk0) );
  }
};

const float pVcf::NAN_FLT = sqrtf(-1.);     // assign float  NAN value
const double pVcf::NAN_DBL = sqrt(-1.);  // assign double NAN value

#endif // __TABIXED_VCF_H
