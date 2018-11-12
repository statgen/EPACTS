#include "Parameters.h"
#include "genomeLoci.h"
#include "genomeScore.bak.h"
#include "Error.h"
#include "fVcf.h"
#include "wFile.h"
#include "hDist.h"
#include "genomePosition.h"
//#include "Rmath.h"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <getopt.h>

#include "IO.h"
#include "TypeConversion.h"
#include "GenomeSequence.h"
#include "cdsStat.h"
#include "ncStat.h"
#include <htslib/bgzf.h>

#define WINDOW_SIZE 65536

class pVPHArgs {
public:
  // VCF-related string arguments
  std::string vcf;
  std::string region;
  std::string rule;
  std::string field;
  std::string scoref;

  // Other input files
  std::string indf;
  std::string itvf;
  std::string bedf;
  std::string genef;
  std::string ref;

  std::string outf;
  int unit;
  bool verbose;
  bool ignoreFilter;
  bool ignoreMissing;
  bool includeMultiAllelic;

  int minAC;
  int minMAC;
  int maxAC;
  double minMAF;
  double maxMAF;
  double minCallRate;

  bool genoFlag;
  bool acFlag;
  bool anFlag;
  bool aldFlag;
  bool tstvFlag;

  bool sepchr;

  static int const DEFAULT_UNIT = 10000L;
  constexpr static double const DEFAULT_MIN = 1e-6;
  constexpr static double const DEFAULT_MAX_MAF = 1;

  pVPHArgs() :
    unit(DEFAULT_UNIT), verbose(false), ignoreFilter(false), ignoreMissing(false), includeMultiAllelic(false), minAC(0), minMAC(0), maxAC(INT_MAX), minMAF(DEFAULT_MIN), maxMAF(DEFAULT_MAX_MAF), minCallRate(DEFAULT_MIN), genoFlag(false), acFlag(false), anFlag(false), aldFlag(false), tstvFlag(false),sepchr(false)
  {}
};

class vcfHashKey {
public:
  std::string genos;
  int ac;
  int an;
  int ald;
  int tstv;

  vcfHashKey(fVcf& vcf, int i, bool genoFlag, bool acFlag, bool anFlag, bool aldFlag, bool tstvFlag, bool ignoreMissing)  {
    if ( genoFlag ) {
      float g;
      genos.resize(vcf.nInds);
      for(int j=0; j < vcf.nInds; ++j) {
        g = vcf.genos[i * vcf.nInds + j];
        if ( std::isnan(g) ) {
          genos[j] = ignoreMissing ? '0' : '.';
        }
        else {
          genos[j] = (int)(g+.5) + '0';
        }
      }
    }
    if ( acFlag ) { ac = vcf.sumAlleles[i]; }
    if ( anFlag ) { an = vcf.numAlleles[i]; }
    if ( aldFlag || tstvFlag ) {
      int i_us = vcf.markers[i].find_first_of('_');
      int i_sl = vcf.markers[i].find_first_of('/',i_us+1);
      int l = (int)vcf.markers[i].size();

      if ( aldFlag ) { ald = (l-i_sl-i_sl+i_us); }
      if ( tstvFlag ) {
        if ( ( i_sl-i_us-1 == 1 ) && ( l-i_sl-1 == 1 ) ) {
          tstv = determineTsTv(vcf.markers[i][i_us+1],vcf.markers[i][i_sl+1]);
        }
        else {
          tstv = -1;
        }
      }
    }
  }

  static int determineTsTv(char ref, char alt) {
    switch(ref) {
      case 'A':
        switch(alt) {
          case 'A': return -1;
          case 'C': return 0;
          case 'G': return 1;
          case 'T': return 0;
          default: return -1;
        }
      case 'C':
        switch(alt) {
          case 'A': return 0;
          case 'C': return -1;
          case 'G': return 0;
          case 'T': return 1;
          default: return -1;
        }
      case 'G':
        switch(alt) {
          case 'A': return 1;
          case 'C': return 0;
          case 'G': return -1;
          case 'T': return 0;
          default: return -1;
        }
      case 'T':
        switch(alt) {
          case 'A': return 0;
          case 'C': return 1;
          case 'G': return 0;
          case 'T': return -1;
          default: return -1;
        }
      default:
        return -1;
    }
  }

  bool operator< (const vcfHashKey& x) const {
    int c = genos.compare(x.genos);
    if ( c != 0 ) {
      return c < 0 ? true : false;
    }
    if ( ac != x.ac ) {
      return ac < x.ac ? true : false;
    }
    if ( an != x.an ) {
      return an < x.an ? true : false;
    }
    if ( ald != x.ald ) {
      return ald < x.ald ? true : false;
    }
    if ( tstv != x.tstv ) {
      return tstv < x.tstv ? true : false;
    }
    return false;
  }

  void print(wFile& wf, bool genoFlag, bool acFlag, bool anFlag, bool aldFlag, bool tstvFlag) const {
    bool tabFlag = false;
    if ( genoFlag ) { wf.printf("%s",genos.c_str()); tabFlag = true; }
    if ( tabFlag ) wf.printf("\t");
    if ( acFlag ) { wf.printf("%d",ac); tabFlag = true; }
    if ( tabFlag ) wf.printf("\t");
    if ( anFlag ) { wf.printf("%d",an); tabFlag = true; }
    if ( tabFlag ) wf.printf("\t");
    if ( aldFlag ) { wf.printf("%d",ald); tabFlag = true; }
    if ( tabFlag ) wf.printf("\t");
    if ( tstvFlag ) {
      switch(tstv) {
        case 0:
          wf.printf("TS");
          break;
        case 1:
          wf.printf("TV");
          break;
        default:
          wf.printf("OTHER");
          break;
      }
      tabFlag = true;
    }
  }
};

void hashBit(fVcf& vcf, unsigned char* bytes, int isnp)  {
  // we fill from the higher order
  float g;
  int nbytes = (vcf.nInds + 3) / 4;
  for(int i=0; i < vcf.nInds; ++i) {
    g = vcf.genos[isnp * vcf.nInds + i];
    unsigned char ng = 0x00;
    if ( ( !std::isnan(g) ) && ( g > 0 ) ) {  // missing genotypes are converted to HOMREF
      ng = (g == 1.) ? 0x01 : 0x02;
    }
    if ( i % 4 == 0 ) bytes[i] = 0;
    bytes[nbytes * isnp + i] |= (ng << (6 - (i % 4) * 2));
  }
}

int runSummary(int argc, char** argv) {
  pVPHArgs arg;
  ParameterList pl;
  double minScore = 0;
  bool fullAFS = false;
  bool nonsnps = false;

  BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("VCF Input Options")
      LONG_STRINGPARAMETER("vcf",&arg.vcf)
      LONG_STRINGPARAMETER("region",&arg.region)
      LONG_STRINGPARAMETER("indf",&arg.indf)
      LONG_STRINGPARAMETER("rule",&arg.rule)
      LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
      LONG_STRINGPARAMETER("scoref",&arg.scoref)
      LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
      LONG_DOUBLEPARAMETER("minScore",&minScore)
      LONG_INTPARAMETER("minAC",&arg.minAC)
      LONG_INTPARAMETER("maxAC",&arg.maxAC)
      LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)
      LONG_PARAMETER("full-afs",&fullAFS)
      LONG_PARAMETER("nonsnps",&nonsnps)

      LONG_PARAMETER_GROUP("Output Options")
      LONG_STRINGPARAMETER("out",&arg.outf)
      LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.vcf.empty() || arg.outf.empty()  ) {
    error("--vcf, --out are required parameters (--indf are also recommended)");
  }

  fVcf tvcf;
  tvcf.load(arg.vcf.c_str(), arg.region.c_str(), "GT", arg.rule.c_str(), !arg.ignoreFilter, arg.indf.empty() ? NULL : arg.indf.c_str());
  int n = tvcf.nInds;
  int i = 0, m = 0, j = 0;
  double af, maf;

  std::vector<int> cnts(n * (3 * 7 + 1), 0); // Ts/Tv/Oth, REF/HET/ALT/SNG/DBL
  std::vector<int> acnts;
  if ( fullAFS ) acnts.resize(n * (n+n+1) * 3, 0); // Ts/Tv/Oth, REF/HET/ALT/SNG/DBL

  //fprintf(stderr,"acnts.size() = %d\n",(int)acnts.size());

  genomeScore gScore;
  if ( !arg.scoref.empty() ) {
    gScore.setDir(arg.scoref.c_str());
  }

  for(int M=0; tvcf.readMarkers(arg.unit); ) {
    M += tvcf.nMarkers;
    fprintf(stderr,"Processing %d markers across %d individuals..\n", M, tvcf.nInds);
    for(i=0, m=0; i < tvcf.nMarkers; ++i) { // for each marker
      af = tvcf.alleleFreq(i);
      maf = af > 0.5 ? 1-af : af;
      if ( ( maf >= arg.minMAF ) &&
           ( tvcf.callRate(i) >= arg.minCallRate ) &&
           ( tvcf.sumAlleles[i] >= arg.minAC ) &&
           ( tvcf.sumAlleles[i] <= arg.maxAC ) &&
           ( arg.scoref.empty() || ( gScore.baseScore(tvcf.chroms[i].c_str(), tvcf.pos1s[i]) >= minScore ) )
        )
      { // if pass the criteria
        int vt = 0; // 0 OTHER 1 Ts 2 Tv
        int frq = 0;  // 0 COMMON 1 SING 2 DBL
        if ( ( tvcf.refs[i].length() == 1 ) && ( tvcf.alts[i].length() == 1 ) ) {
          switch(tvcf.refs[i][0]) {
            case 'A':
              if ( tvcf.alts[i][0] == 'G' ) { vt = 1; }
              else if ( tvcf.alts[i][0] == 'C' ) { vt = 2; }
              else if ( tvcf.alts[i][0] == 'T' ) { vt = 2; }
              break;
            case 'C':
              if ( tvcf.alts[i][0] == 'T' ) { vt = 1; }
              else if ( tvcf.alts[i][0] == 'A' ) { vt = 2; }
              else if ( tvcf.alts[i][0] == 'G' ) { vt = 2; }
              break;
            case 'G':
              if ( tvcf.alts[i][0] == 'A' ) { vt = 1; }
              else if ( tvcf.alts[i][0] == 'C' ) { vt = 2; }
              else if ( tvcf.alts[i][0] == 'T' ) { vt = 2; }
              break;
            case 'T':
              if ( tvcf.alts[i][0] == 'C' ) { vt = 1; }
              else if ( tvcf.alts[i][0] == 'A' ) { vt = 2; }
              else if ( tvcf.alts[i][0] == 'G' ) { vt = 2; }
              break;
          }
        }

        frq = (tvcf.sumAlleles[i] == 1) ? 1 : (tvcf.sumAlleles[i] == 2) ? 2 : 0;
        for(j=0; j < n; ++j) {
          float v = tvcf.genos[i*n + j];
          if ( !std::isnan(v) ) {
            int g = (v == 0 ? 0 : ( v == 1 ) ? 1 : 2 );
            ++cnts[g * 3 + vt + j * 22];
            if ( ( frq > 0 ) && ( g > 0 ) ) {
              ++cnts[6 + frq * 3 + vt + j * 22];
            }
            if ( g > 0 ) {
              if ( af >= 0.005 ) {
                ++cnts[18 + vt + j * 22];
              }
              else {
                ++cnts[15 + vt + j * 22];
              }
            }
            if ( fullAFS ) ++acnts[j*(n+n+1)*3 + (int)(tvcf.sumAlleles[i])*3 + g];
          }
          else {
            ++cnts[21 + j * 22];
          }
        }
        ++m;
      }
    }
  }

  wFile wf(arg.outf.c_str());
  wf.printf("#IND_ID\t#SNPS\t#SING\t#DBL\t#AF<.5%\tAF>.5%\t#REF\t#HET\t#ALT\tMISS\tREFTsTv\tHETTsTv\tALTTsTv\tSNGTsTv\tRARTsTv\n");
  for(j=0; j < n; ++j) {
    wf.printf("%s",tvcf.inds[j].c_str());
    wf.printf("\t%d",cnts[j*22+3]+cnts[j*22+4]+cnts[j*22+5]+cnts[j*22+6]+cnts[j*22+7]+cnts[j*22+8]);
    wf.printf("\t%d",cnts[j*22+9]+cnts[j*22+10]+cnts[j*22+11]);
    wf.printf("\t%d",cnts[j*22+12]+cnts[j*22+13]+cnts[j*22+14]);
    wf.printf("\t%d",cnts[j*22+15]+cnts[j*22+16]+cnts[j*22+17]);
    wf.printf("\t%d",cnts[j*22+18]+cnts[j*22+19]+cnts[j*22+20]);
    wf.printf("\t%d",cnts[j*22+0]+cnts[j*22+1]+cnts[j*22+2]);
    wf.printf("\t%d",cnts[j*22+3]+cnts[j*22+4]+cnts[j*22+5]);
    wf.printf("\t%d",cnts[j*22+6]+cnts[j*22+7]+cnts[j*22+8]);
    wf.printf("\t%d",cnts[j*22+21]);
    wf.printf("\t%.3lf",(double)(cnts[j*22+1]+.1)/(cnts[j*22+2]+.2));
    wf.printf("\t%.3lf",(double)(cnts[j*22+4]+.1)/(cnts[j*22+5]+.2));
    wf.printf("\t%.3lf",(double)(cnts[j*22+7]+.1)/(cnts[j*22+8]+.2));
    wf.printf("\t%.3lf",(double)(cnts[j*22+10]+.1)/(cnts[j*22+11]+.2));
    wf.printf("\t%.3lf",(double)(cnts[j*22+16]+.1)/(cnts[j*22+17]+.2));
    wf.printf("\n");
  }
  wf.close();

  if ( fullAFS ) {
    wFile acf((arg.outf+".indAFS").c_str());
    acf.printf("#IND_ID");
    for(j=0; j <= n+n; ++j) {
      acf.printf("\t%d",j);
    }
    acf.printf("\n");

    for(j=0; j < n; ++j) {
      acf.printf("%s",tvcf.inds[j].c_str());
      for(m=0; m <= n+n; ++m) {
        acf.printf("\t%d,%d,%d",acnts[j*(n+n+1)*3+m*3+0],acnts[j*(n+n+1)*3+m*3+1],acnts[j*(n+n+1)*3+m*3+2]);
      }
      acf.printf("\n");
    }
    acf.close();
  }

  return 0;
}

int schr2nchr(const char* schr) {
  int n = atoi(schr);
  if ( n > 0 ) { return n; }
  else if ( strcmp(schr,"X") == 0 ) { return 23; }
  else if ( strcmp(schr,"Y") == 0 ) { return 24; }
  else if ( strcmp(schr,"XY") == 0 ) { return 25; }
  else if ( ( strcmp(schr,"MT") == 0 ) || ( strcmp(schr,"M") == 0 ) || ( strcmp(schr,"mt") == 0 ) ) { return 26; }
  else if ( strncmp(schr,"chr",3) == 0 ) { return schr2nchr(schr+3); }
  else { error("Cannot convert chromosome %s into PLINK format",schr); return 0; }
}


// variable threshold test
int runConvert(int argc, char** argv) {
  error("conversion temporarily disabled!");
  return -1;
  // Parse the input arguments
  bool outvcf = false;
  bool outplink = false;
  bool outmatrix = true;
  double maxR2 = 1;
  int winR2 = 0;
  pVPHArgs arg;
  arg.field = "GT";
  std::string markerId;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("VCF Input Options")
      LONG_STRINGPARAMETER("vcf",&arg.vcf)
      LONG_STRINGPARAMETER("indf",&arg.indf)
      LONG_STRINGPARAMETER("bedf",&arg.bedf)
      LONG_STRINGPARAMETER("field",&arg.field)
      LONG_STRINGPARAMETER("region",&arg.region)
      LONG_STRINGPARAMETER("marker-id",&markerId)
      LONG_STRINGPARAMETER("rule",&arg.rule)
      LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
      LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
      LONG_INTPARAMETER("minMAC",&arg.minMAC)
      LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
      LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)
      LONG_PARAMETER("sepchr",&arg.sepchr)

      LONG_PARAMETER_GROUP("LD pruning Options")
      LONG_INTPARAMETER("win-r2",&winR2)
      LONG_DOUBLEPARAMETER("max-r2",&maxR2)

      LONG_PARAMETER_GROUP("Output format")
      EXCLUSIVE_PARAMETER("outvcf",&outvcf)
      EXCLUSIVE_PARAMETER("outplink",&outplink)
      EXCLUSIVE_PARAMETER("outmatrix",&outmatrix)

      LONG_PARAMETER_GROUP("Output Options")
      LONG_STRINGPARAMETER("out",&arg.outf)
      LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.vcf.empty() || arg.outf.empty()  ) {
    error("--vcf, --out are required parameters (--indf are also recommended)");
  }

  if ( !( outvcf || outplink || outmatrix ) ) {
    warning("Default option --outmatrix was turned off. Turning in again..");
    outmatrix = true;
  }

  bool gtFlag = (arg.field == "GT");
  fVcf tvcf;

  if ( ( !arg.region.empty() ) && ( !markerId.empty() ) ) {
    error("--region and --marker-id cannot be combined together");
  }

  std::string markerIdRef, markerIdAlt;
  if ( !markerId.empty() ) {
    std::vector<std::string> tokens;
    pFile::tokenizeLine(markerId.c_str(),":_/",tokens);
    arg.region = tokens[0] + ":" + tokens[1] + "-" + tokens[1];
    if ( tokens.size() > 2 ) {
      markerIdRef = tokens[2];
      markerIdAlt = tokens[3];
    }
  }

  if ( ( !arg.region.empty() ) && ( !arg.bedf.empty() ) ) {
    error("--region and --bedf cannot be combined together");
  }

  tvcf.load(arg.vcf.c_str(), arg.region.c_str(), arg.field.c_str(), arg.rule.c_str(), !arg.ignoreFilter, arg.indf.empty() ? NULL : arg.indf.c_str());

  int n = tvcf.nInds;
  int i,m,j,k;
  wFile *bimf = NULL;
  if ( outplink) {
    if ( !gtFlag ) error("--field GT must be set with --outplink option");
    //error("--outplink is not implemented yet");
    bimf = new wFile((arg.outf+".bim").c_str());
    wFile famf((arg.outf+".fam").c_str());
    for(i=0; i < n; ++i) {
      famf.printf("%d\t%s\t0\t0\t0\t-9\n",tvcf.inds[i].c_str(),tvcf.inds[i].c_str());
    }
    famf.close();
    arg.outf += ".bed";
  }
  wFile wf(arg.outf.c_str());
  if ( outplink ) {
    //error("--outplink is not implemented yet");
    // write headers
    char magicNumbers[3] = {0x6c,0x1b,0x01};
    wf.write(magicNumbers,3);
  }
  else if ( outmatrix ) {
    //wf.printf("#MARKER");
    for(j=0; j < n; ++j) {
      if ( j > 0 ) wf.printf("\t");
      wf.printf("%s",tvcf.inds[j].c_str());
    }
    wf.printf("\n");
  }
  else if ( outvcf ) {
    //error("--outvcf is not implemented yet");
    // write headers
    for(j=0; j < (int)tvcf.headers.size(); ++j) {
      wf.printf("%s\n",tvcf.headers[j].c_str());
    }
    wf.printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for(j=0; j < n; ++j) {
      wf.printf("\t%s",tvcf.inds[j].c_str());
    }
    wf.printf("\n");
    arg.unit = 1;
  }

  float geno;
  int ngeno;
  double af, maf;
  int nbytes = (n+3)/4;
  char* genos = outplink ? new char[nbytes]() : NULL;

  genomeLoci loci;
  if ( !arg.region.empty() ) {
    loci.add(arg.region.c_str());
  }
  else if ( !arg.bedf.empty() ) {
    const char* line = NULL;
    pFile bedf(arg.bedf.c_str());
    std::vector<std::string> tokens;
    while( (line = bedf.getLine()) != NULL ) {
      pFile::tokenizeLine(line, " \t\r\n", tokens);
      if ( tokens.size() < 3 )
        error("bed file must have three columns");

      if ( tokens[0].compare(0,3,"chr") == 0 ) {
        tokens[0] = tokens[0].substr(3);
      }

      int beg1 = atoi(tokens[1].c_str()) + 1;
      int end0 = atoi(tokens[2].c_str());

      loci.add(tokens[0].c_str(),beg1,end0);
    }
    loci.resolveOverlaps(); // make the bed files non-overlapping
  }

  if ( !loci.empty() ) {
    loci.rewind();
  }

  std::vector<double> r2G;
  std::vector<bool> r2P;
  if ( winR2 > 0 ) {
    r2G.resize(winR2 * n, 0);
    r2P.resize(winR2, false);
  }

  int M = 0, nout = 0;
  double mu, sigma, r2 = 0;
  do {
    if ( !loci.empty() ) {
      tvcf.updateRegion(loci.currentLocus().toString(), arg.sepchr);
      //notice("Extracing region %s..",loci.currentLocus().toString());
    }

    for(m = 0; tvcf.readMarkers(arg.unit); ) {
      if ( M / 10000 < (M + tvcf.nMarkers) / 10000 ) {
        fprintf(stderr,"Reading %d of %d markers and writing %d markers across %d individuals..\n", m, M + tvcf.nMarkers, nout, tvcf.nInds);
      }
      M += tvcf.nMarkers;
      //if ( ( arg.unit > 1 ) || ( M % 10000 == 0 ) )
      for(i=0; i < tvcf.nMarkers; ++i) { // for each marker
        af = tvcf.alleleFreq(i);
        maf = af > 0.5 ? 1-af : af;
        if ( maf < 0 ) maf = 0;
        mu = 2*af;
        sigma = tvcf.alleleSD(i,true);

        if ( !markerIdRef.empty() ) {
          if ( ( tvcf.refs[i] != markerIdRef ) || ( tvcf.alts[i] != markerIdAlt ) ) {
            continue;
          }
        }

        if ( ( maf >= arg.minMAF ) &&
             ( tvcf.callRate(i) >= arg.minCallRate ) &&
             ( tvcf.sumAlleles[i] >= arg.minAC ) &&
             ( tvcf.sumAlleles[i] <= arg.maxAC ) )
        { // if pass the criteria
          if ( winR2 > 0 ) {
            // fill the genotype matrix
            int offset = (m % winR2) * n;
            for(j=0; j < n; ++j) {
              geno = tvcf.genos[j + n*i];
              if ( std::isnan(geno) ) {
                r2G[offset+j] = 0;
              }
              else {
                r2G[offset+j] = (geno-mu)/sigma;
              }
            }
            // calculate sliding-window r2
            r2 = 0;
            for(j=0; j < winR2; ++j) {
              if ( ( r2P[j] ) && ( j != (m % winR2) ) ) {
                double r = 0;
                for(k=0; k < n; ++k) {
                  r += ((r2G[offset+k] * r2G[j*n + k])/n);
                }
                if ( r*r > r2 ) r2 = r*r;
              }
            }
            if ( r2 < maxR2 ) { r2P[m % winR2] = true; }
            else { r2P[m % winR2] = false; }
            //fprintf(stderr,"r2=%lf\n",r2);
          }

          if ( ( winR2 == 0 ) || ( r2P[m % winR2] ) ) {
            if ( outplink ) {
              //error("--outplink is not implemented yet");
              for(j=0; j < n; ++j) {
                geno = tvcf.genos[j + n*i];
                ngeno = std::isnan(geno) ? 0 : ((int)floor(geno)+1);
                switch(ngeno) {
                  case 0:
                    genos[j/4] |=  (0x1 << ((j%4)*2));
                    break;
                  case 1:
                    genos[j/4] |=  (0x0 << ((j%4)*2));
                    break;
                  case 2:
                    genos[j/4] |=  (0x2 << ((j%4)*2));
                    break;
                  case 3:
                    genos[j/4] |=  (0x3 << ((j%4)*2));
                    break;
                  default:
                    warning("Cannot correctly parse non-biallelic SNP, treating as HOMALT");
                    genos[j/4] |=  (0x3 << ((j%4)*2));
                }
              }
              wf.write(genos,nbytes);
              memset(genos,0,nbytes);
              bimf->printf("%d\t%s\t0\t%d\t%s\t%s\n",schr2nchr(tvcf.chroms[i].c_str()),tvcf.markers[i].c_str(),tvcf.pos1s[i],tvcf.refs[i].c_str(),tvcf.alts[i].c_str());
            }
            else if ( outmatrix ) {
              wf.printf("%s",tvcf.markers[i].c_str());
              for(j=0; j < n; ++j) {
                geno = tvcf.genos[j + n*i];
                if ( std::isnan(geno) ) {
                  wf.printf("\tNA");
                }
                else {
                  if ( gtFlag ) {
                    wf.printf("\t%d",(int)floor(geno));
                  }
                  else {
                    wf.printf("\t%.3lf",geno);
                  }
                }
              }
              wf.printf("\n");
            }
            else if ( outvcf ) {
              //error("--outvcf is not implemented yet");
              //fprintf(stderr,"foo %d\n",M);
              tvcf.writeSubsetMarker(wf, ""); //tvcf.tf.peekLine()); TODO: Reimplement
              //fprintf(stderr,"bar %d\n",M);
              //if ( arg.verbose)
              //fprintf(stderr,"fVcf::writeSubsetMarker() finished, M= %d",M);
            }
            ++nout;
          }
          ++m;
        }
      }
    }
  }
  while ( (!loci.empty()) && loci.next() );

  wf.close();
  if ( outplink ) {
    bimf->close();
    delete bimf;
    delete [] genos;
  }
  return 0;

}

// variable threshold test
int runPairLD(int argc, char** argv) {
  double minR2 = 1;
  int winR2 = 1000;
  pVPHArgs arg;
  bool ignorePhase = false;
  arg.field = "GT";
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("VCF Input Options")
      LONG_STRINGPARAMETER("vcf",&arg.vcf)
      LONG_STRINGPARAMETER("indf",&arg.indf)
      LONG_STRINGPARAMETER("field",&arg.field)
      LONG_STRINGPARAMETER("region",&arg.region)
      LONG_STRINGPARAMETER("rule",&arg.rule)
      LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
      LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
      LONG_INTPARAMETER("minMAC",&arg.minMAC)
      LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
      LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)
      LONG_PARAMETER("ignorePhase",&ignorePhase)

      LONG_PARAMETER_GROUP("LD pruning Options")
      LONG_INTPARAMETER("win-r2",&winR2)
      LONG_DOUBLEPARAMETER("min-r2",&minR2)

      LONG_PARAMETER_GROUP("Output Options")
      LONG_STRINGPARAMETER("out",&arg.outf)
      LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.vcf.empty() || arg.outf.empty()  ) {
    error("--vcf, --out are required parameters (--indf are also recommended)");
  }

  bool gtFlag = (arg.field == "GT");
  fVcf tvcf;

  tvcf.load(arg.vcf.c_str(), arg.region.c_str(), arg.field.c_str(), arg.rule.c_str(), !arg.ignoreFilter, arg.indf.empty() ? NULL : arg.indf.c_str());

  int n = tvcf.nInds;
  std::vector<double> r2G;
  std::vector<double> r2AF;
  std::vector<std::string> r2ID;
  if ( winR2 > 0 ) {
    r2G.resize(winR2 * n * 2, 0);
    r2ID.resize(winR2);
    r2AF.resize(winR2);
  }

  int i, j, k, m, M = 0, ngeno, phase, offset, offset2;
  double mu, af, maf, sigma, sigmaA, r, r2 = 0, h0, h1;
  bool phased = false;
  float geno;

  wFile wf(arg.outf.c_str());
  wf.printf("#MARKER1\tMARKER2\tAF1\tAF2\tR2\tR\n");

  for(m = 0; tvcf.readMarkers(arg.unit); ) {
    if ( M / 10000 < (M + tvcf.nMarkers) / 10000 ) {
      fprintf(stderr,"Reading %d / %d markers across %d individuals..\n", m, M + tvcf.nMarkers, tvcf.nInds);
    }
    M += tvcf.nMarkers;

    //if ( m == 0 ) { phased = (gtFlag && (tvcf.phases[0] > 0) && (!ignorePhase)); }
    //if ( ( arg.unit > 1 ) || ( M % 10000 == 0 ) )
    for(i=0; i < tvcf.nMarkers; ++i) { // for each marker
      af = tvcf.alleleFreq(i);
      maf = af > 0.5 ? 1-af : af;
      mu = 2*af;

      //notice("i=%d, phased = %d",i,phased);

      if ( ( maf >= arg.minMAF ) &&
           ( tvcf.callRate(i) >= arg.minCallRate ) &&
           ( tvcf.sumAlleles[i] >= arg.minAC ) &&
           ( tvcf.sumAlleles[i] <= arg.maxAC ) )
      {

        offset = (m % winR2) * (n*2);
        r2ID[m % winR2] = tvcf.markers[i];
        r2AF[m % winR2] = af;
        sigma = sqrt(af*(1.-af));
        sigmaA = tvcf.alleleSD(i,true);
        h0 = (0-af)/sigma;
        h1 = (1.-af)/sigma;
        for(j=0; j < n; ++j) {
          phased = (gtFlag && (tvcf.phases[j + n*i] > 0) && (!ignorePhase));
          if ( phased ) {
            geno = tvcf.genos[j + n*i];
            if ( std::isnan(geno) ) { error("Missing genotype was observed in phased VCF"); }
            ngeno = (int)geno;
            phase = tvcf.phases[j + n*i];
            switch(ngeno) {
              case 0:
                r2G[offset+j+j] = r2G[offset+j+j+1] = h0;
                break;
              case 1:
                if ( phase == 2 ) { r2G[offset+j+j] = h0; r2G[offset+j+j+1] = h1; }
                else if ( phase == 3 ) { r2G[offset+j+j] = h1; r2G[offset+j+j+1] = h0; }
                else { error("Incompatible phase information %d at heterozygous genotypes",phase); }
                break;
              case 2:
                r2G[offset+j+j] = r2G[offset+j+j+1] = h1;
                break;
              default:
                error("Incompatible genotype %d in phased VCF",ngeno);
            }
          }
          else {
            geno = tvcf.genos[j + n*i];
            if ( std::isnan(geno) ) {
              r2G[offset+j+j] = r2G[offset+j+j+1] = 0;
            }
            else {
              r2G[offset+j+j] = r2G[offset+j+j+1] = (geno-mu)/sigmaA;
            }
          }
        }

// 	  offset = (m % winR2) * (n*2);
// 	  r2ID[m % winR2] = tvcf.markers[i];
// 	  r2AF[m % winR2] = af;
// 	  if ( phased ) {
// 	    sigma = sqrt(af*(1.-af));
// 	    for(j=0; j < n; ++j) {
// 	      geno = tvcf.genos[j + n*i];
// 	      if ( std::isnan(geno) ) { error("Missing genotype was observed in phased VCF"); }
// 	      ngeno = (int)geno;
// 	      phase = tvcf.phases[j + n*i];
// 	      switch(ngeno) {
// 	      case 0:
// 		r2G[offset+j+j] = r2G[offset+j+j+1] = (0-af)/sigma;
// 		break;
// 	      case 1:
// 		if ( phase == 2 ) { r2G[offset+j+j] = (0-af)/sigma; r2G[offset+j+j+1] = (1.-af)/sigma; }
// 		else if ( phase == 3 ) { r2G[offset+j+j] = (1.-af)/sigma; r2G[offset+j+j+1] = (0-af)/sigma; }
// 		else { error("Incompatible phase information %d at heterozygous genotypes",phase); }
// 		break;
// 	      case 2:
// 		r2G[offset+j+j] = r2G[offset+j+j+1] = (1.-af)/sigma;
// 		break;
// 	      default:
// 		error("Incompatible genotype %d in phased VCF",ngeno);
// 	      }
// 	    }
// 	  }
// 	  else {
// 	    sigma = tvcf.alleleSD(i,true);
//       	    for(j=0; j < n; ++j) {
// 	      geno = tvcf.genos[j + n*i];
// 	      if ( std::isnan(geno) ) {
// 		r2G[offset+j] = 0;
// 	      }
// 	      else {
// 		r2G[offset+j] = (geno-mu)/sigma;
// 	      }
// 	    }
// 	  }

        // calculate sliding-window r2 matrix
        r2 = r = 0;
        for(j=1; j < winR2; ++j) {
          //notice("foo %d %d %d %d",j,offset2,k);
          offset2 = ((m + j) % winR2) * (n*2);
          if ( r2AF[(m+j) % winR2] > 0 ) {
            for(k=0; k < n + n; ++k) {
              //for(k=0; k < (phased ? n+n : n); ++k) {
              r += (r2G[offset + k] * r2G[offset2 + k]);
            }
            r /= (double)(n+n);

            // handle boundary condition
            if ( r > 1 ) r = 1.;
            else if ( r < -1 ) r = -1.;

            r2 = r*r;
            if ( r2 >= minR2 ) {
              wf.printf("%s\t%s\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n",r2ID[(m+j) % winR2].c_str(),r2ID[m % winR2].c_str(),r2AF[(m+j) % winR2],r2AF[m % winR2],r2,r);
            }
          }
        }
        ++m;
      }
    }
  }
  wf.close();
  return 0;
}

// create the list of r2 and 
int runIndexLD(int argc, char** argv) {
  // Parse the input arguments
  std::string index;
  double minR2 = 0;
  int win = 1000000;
  pVPHArgs arg;
  arg.field = "GT";
  bool ignorePhase = false;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("VCF Input Options")
      LONG_STRINGPARAMETER("vcf",&arg.vcf)
      LONG_STRINGPARAMETER("indf",&arg.indf)
      LONG_STRINGPARAMETER("field",&arg.field)
      LONG_STRINGPARAMETER("region",&arg.region)
      LONG_STRINGPARAMETER("rule",&arg.rule)
      LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
      LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
      LONG_INTPARAMETER("minMAC",&arg.minMAC)
      LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
      LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)
      LONG_PARAMETER("ignorePhase",&ignorePhase)
      LONG_PARAMETER("sepchr",&arg.sepchr)

      LONG_PARAMETER_GROUP("Index SNP position")
      LONG_STRINGPARAMETER("index",&index)
      LONG_INTPARAMETER("win",&win)
      LONG_DOUBLEPARAMETER("min-r2",&minR2)

      LONG_PARAMETER_GROUP("Output Options")
      LONG_STRINGPARAMETER("out",&arg.outf)
      LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.vcf.empty() || index.empty() || arg.outf.empty()  ) {
    error("--vcf, --index, --out are required parameters (--indf are also recommended)");
  }

  fVcf tvcf;

  // parse index SNP info and convert it into regions
  std::vector<std::string> tokens;
  pFile::tokenizeLine(index.c_str(), ":", tokens);
  if ( tokens.size() != 2 ) { error("Cannot parse --indx %s parameter",index.c_str()); }
  std::string indexChr = tokens[0];
  std::string indexSPos = tokens[1];
  int indexNPos = atoi(indexSPos.c_str());
  std::string indexRegion = index + "-" + indexSPos;

  // load index SNP
  tvcf.load(arg.vcf.c_str(), NULL, arg.field.c_str(), arg.rule.c_str(), !arg.ignoreFilter, arg.indf.empty() ? NULL : arg.indf.c_str());
  int n = tvcf.nInds;

  tvcf.updateRegion(indexRegion.c_str(), arg.sepchr);

  if ( tvcf.readMarkers(1) == 0 ) {
    error("Cannot find a SNP at marker position %s",index.c_str());
  }

  int j, k, m;
  double af = 0, maf, r, mu, sigma, sigmaA;
  int ngeno, phase;
  float geno;
  std::vector<double> iG(n*2,0);
  std::vector<double> tG(n*2,0);
  bool phased; // = (( tvcf.phases[0] > 0 ) && ( !ignorePhase ));
  double h0, h1;

  af = tvcf.alleleFreq(0);
  sigma = sqrt(af*(1.-af));
  mu = 2*af;
  sigmaA = tvcf.alleleSD(0,true);
  h0 = (0.0-af)/sigma;
  h1 = (1.0-af)/sigma;
  //notice("af=%lf sigma=%lf sigmaA=%lf",af,sigma,sigmaA);
  for(j=0; j < n; ++j) {
    geno = tvcf.genos[j];
    phased = (( tvcf.phases[j] > 0 ) && ( !ignorePhase ));
    if ( phased ) {
      if ( std::isnan(geno) ) { error("Missing genotype was observed in phased VCF"); }
      ngeno = (int)geno;
      phase = tvcf.phases[j];
      switch(ngeno) {
        case 0:
          iG[j+j] = iG[j+j+1] = h0;
          break;
        case 1:
          if ( phase == 2 ) { iG[j+j] = h0; iG[j+j+1] = h1; }
          else if ( phase == 3 ) { iG[j+j] = h1; iG[j+j+1] = h0; }
          else { error("Incompatible phase information %d at heterozygous genotypes",phase); }
          break;
        case 2:
          iG[j+j] = iG[j+j+1] = h1;
          break;
        default:
          error("Incompatible genotype %d in phased VCF",ngeno);
      }
      //af = (double)cd/(double)(ab+cd);
    }
    else {
      for(j=0; j < n; ++j) {
        geno = tvcf.genos[j];
        if ( std::isnan(geno) ) { iG[j+j] = iG[j+j+1] = 0; }
        else { iG[j+j] = iG[j+j+1] = (geno-mu)/sigmaA; }
      }
    }
  }

  notice("Found index SNP %s in the VCF at AF=%.5lf",index.c_str(),af);

  if ( arg.region.empty() ) {
    char buf[255];
    sprintf(buf,"%s:%d-%d",indexChr.c_str(),indexNPos-win,indexNPos+win);
    arg.region = buf;
  }

  notice("Retrieving the genomic region %s",arg.region.c_str());
  tvcf.updateRegion(arg.region.c_str(), arg.sepchr);

  wFile wf(arg.outf.c_str());
  wf.printf("##INDEX=%s\n",index.c_str());
  wf.printf("##AF=%.5lf\n",af);
  wf.printf("#CHROM\tPOS\tID\tREF\tALT\tAF\tR2\tR\n");

  int nout = 0;
  int M = 0;
  for(m = 0; tvcf.readMarkers(arg.unit); ) {
    if ( M / 10000 < (M + tvcf.nMarkers) / 10000 ) {
      fprintf(stderr,"Reading %d of %d markers and writing %d markers across %d individuals..\n", m, M + tvcf.nMarkers, nout, tvcf.nInds);
    }
    M += tvcf.nMarkers;

    for(int i=0; i < tvcf.nMarkers; ++i) { // for each marker
      af = tvcf.alleleFreq(i);
      maf = af > 0.5 ? 1-af : af;

      if ( ( maf >= arg.minMAF ) &&
           ( tvcf.callRate(i) >= arg.minCallRate ) &&
           ( tvcf.sumAlleles[i] >= arg.minAC ) &&
           ( tvcf.sumAlleles[i] <= arg.maxAC ) ) { // if pass the criteria
        mu = 2*af;
        sigma = sqrt(af*(1.-af));
        sigmaA = tvcf.alleleSD(0,true);
        h0 = (0.0-af)/sigma;
        h1 = (1.0-af)/sigma;
        for(j=0; j < n; ++j) {
          geno = tvcf.genos[j + n*i];
          phased = (( tvcf.phases[j + n*i] > 0 ) && ( !ignorePhase ));
          if ( phased ) { // assume phased, non-missing genotypes
            if ( std::isnan(geno) ) { error("Missing genotype was observed in phased VCF"); }
            ngeno = (int)geno;
            phase = tvcf.phases[j + n*i];
            switch(ngeno) {
              case 0:
                tG[2*j] = tG[2*j+1] = h0;
                break;
              case 1:
                if ( phase == 2 ) { tG[j+j] = h0; tG[j+j+1] = h1; }
                else if ( phase == 3 ) { tG[j+j] = h1; tG[j+j+1] = h0; }
                else { error("Incompatible phase information %d at heterozygous genotypes",phase); }
                break;
              case 2:
                tG[j+j] = tG[j+j+1] = h1;
                break;
              default:
                error("Incompatible genotype %d in phased VCF",ngeno);
            }
          }
          else {
            geno = tvcf.genos[j + n*i];
            if ( std::isnan(geno) ) {
              tG[j+j] = tG[j+j+1] = 0;
            }
            else {
              tG[j+j] = tG[j+j+1] = (geno-mu)/sigmaA;
            }
          }
        }

        // calculate r2 with index SNP
        r = 0;
        for(k=0; k < n+n; ++k) {
          r += (tG[k] * iG[k]);
          //notice("%d\t%lf\t%lf",k,iG[k],tG[k]);
        }
        r /= (double)(n+n);
        //error("r=%lf",r);
        if ( r*r >= minR2 ) { // print out R2
          wf.printf("%s\t%d\t.\t%s\t%s\t%.5lf\t%.5lf\t%.5lf\n",tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i].c_str(), tvcf.alts[i].c_str(), af, r*r, r);
          ++nout;
        }
      }
    }
  }

  notice("Successfully wrote %d variants at r2 threshold %f",nout,minR2);
  wf.close();
  return 0;
}

// create the list of r2 and 
int runPeakShift(int argc, char** argv) {
  // Parse the input arguments
  std::string score;  // score file contains [CHROM] [POS1] [SCOREF]
  std::string bed;    // BED file contains [CHROM] [BEG-0] [END-0] [SCORE or CATEGORY]
  std::string mpu;    // DEPTH file contains [CHROM] [POS] [REF] [DEPTH]
  int win = 1000000;  // Window to extend from the index bp
  std::string region; // Region to be explorered
  std::string chrpos; // chromposomal position of index SNPs
  //bool quantitiative = false; // use quantitative scores
  int maxShift = 100000; // maximum shift 
  int minDepth = 1;   // minimum depth
  std::string mask;   // genome mask
  std::string out;    // output files
  int seed = 0;       // random seed
  int nperm = 1000;   // number of permutations
  bool ucsc = false;
  char buf[65535];
  bool noAvgFlag = false;

  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Input Files")
      LONG_STRINGPARAMETER("score",&score)
      LONG_STRINGPARAMETER("bed",&bed)
      LONG_STRINGPARAMETER("mpu",&mpu)
      LONG_STRINGPARAMETER("mask",&mask)
      LONG_PARAMETER("ucsc",&ucsc)

      LONG_PARAMETER_GROUP("Input Options")
      LONG_PARAMETER("no-avg",&noAvgFlag)
      LONG_INTPARAMETER("seed",&seed)
      LONG_INTPARAMETER("nperm",&nperm)

      LONG_PARAMETER_GROUP("Regions to focus")
      LONG_STRINGPARAMETER("region",&region)
      LONG_STRINGPARAMETER("chrpos",&chrpos)
      LONG_INTPARAMETER("win",&win)
      LONG_INTPARAMETER("max-shift",&maxShift)
      LONG_INTPARAMETER("min-depth",&minDepth)

      LONG_PARAMETER_GROUP("Output Options")
      LONG_STRINGPARAMETER("out",&out)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( score.empty() || out.empty() ) {
    error("--score and --out are required parameters");
  }
  if ( mpu.empty() && bed.empty() ) {
    error("Either --mpu or --bed are required parameters");
  }
  if ( !mpu.empty() && !bed.empty() ) {
    error("Only one of --mpu or --bed are required parameters");
  }
  if ( !region.empty() && !chrpos.empty() ) {
    error("Only one of --region or --chrpos can be specified");
  }

  // open the score file first
  std::vector<std::string> tokens;
  if ( !chrpos.empty() ) {
    pFile::tokenizeLine(chrpos.c_str(),":",tokens);
    std::string chr = tokens[0];
    int pos = atoi(tokens[1].c_str());
    int beg = pos - win;
    int end = pos + win;
    sprintf(buf,"%s%s:%d-%d",ucsc ? "chr" : "",chr.c_str(),beg,end);
    region = buf;
  }

  pFile scoref(score.c_str());
  //if ( !region.empty() ) scoref.updateRegion(region.c_str());

  std::map<int,double> pos2score;
  std::string chr;
  std::string bedchr;
  int minpos = 1e9;
  int maxpos = 0;
  int firstpos, newpos;
  const char* line = NULL;
  while( ( line = scoref.getLine() ) != NULL ) {
    pFile::tokenizeLine(line, " \t\r\n", tokens);
    if ( tokens.size() < 3 ) error("bed file must have three columns");

    if ( chr.empty() ) { chr = tokens[0]; }
    else if ( chr != tokens[0] ) error("Multiple chromosomes - %s and %s - are not supported",chr.c_str(),tokens[0].c_str());
    int bp = atoi(tokens[1].c_str());
    if ( minpos > bp ) minpos = bp;
    if ( maxpos < bp ) maxpos = bp;
    pos2score[bp] = atof(tokens[2].c_str());
  }
  bedchr = ucsc ? "chr"+chr : chr;
  int length = (maxpos-minpos+1);
  sprintf(buf,"%s:%d-%d",bedchr.c_str(),minpos,maxpos);
  region = buf;

  //error("%s:%d-%d",bedchr.c_str(),minpos,maxpos);

  // read BED file 
  std::map<std::string, genomeLoci> mLoci;  // BED regions with predicted category
  if ( !bed.empty() ) {
    notice("Reading BED file...%s",bed.c_str());
    pFile bedf(bed.c_str(), region.c_str());
    int sz = 0;
    int nregions = 0;
    while( (line = bedf.getLine()) != NULL ) {
      pFile::tokenizeLine(line, " \t\r\n", tokens);
      std::string key("ALL");
      //double f = 1.;
      if ( tokens.size() < 3 ) {
        error("BED file must have at least three columns");
      }
      else if ( tokens.size() > 3 ) {
        key = tokens[3];
        //if ( tokens.size() > 4 ) {
        //f = atof(tokens[4].c_str());
        //}
      }
      int beg = atoi(tokens[1].c_str())+1;
      int end = atoi(tokens[2].c_str());
      if ( beg < minpos ) beg = minpos;
      if ( end > maxpos ) end = maxpos;
      sz += (end-beg+1);
      mLoci[key].add(tokens[0].c_str(), beg, end);
      ++nregions;
    }
    notice("Identified %d regions of %d base pairs across %d categories", nregions, sz, mLoci.size());
  }
  else { // PEAK file should be available
    notice("Reading peak files...\n");
    pFile mpuf(mpu.c_str(), region.c_str());
    while( ( line = mpuf.getLine() ) != NULL ) {
      if ( tokens.size() < 4 ) {
        error("MPU file must have at least four columns");
      }
      std::string key("ALL");
      int pos = atoi(tokens[1].c_str());
      if ( ( pos < minpos ) || ( pos > maxpos ) ) continue;
      int depth = atoi(tokens[3].c_str());
      if ( depth >= minDepth )
        mLoci[key].add(tokens[0].c_str(), pos, pos);
    }
    notice("Identified region of %d base pairs with depth %d or greater",mLoci.size(),minDepth);
  }
  notice("Resolving Overlaps...");
  for(std::map<std::string,genomeLoci>::iterator it = mLoci.begin(); it != mLoci.end(); ++it) {
    it->second.resolveOverlaps();
    //notice("%s = %d",it->first.c_str(), it->second.totalLength());
  }

  notice("Generating random offsets...\n");
  // Generate random offsets
  if ( seed == 0 ) srand(std::time(0));
  else srand(seed);
  std::vector<int> offsets;
  offsets.push_back(0);
  for(int i=1; i < nperm; ++i) {
    offsets.push_back( rand() % maxShift - maxShift / 2 );
  }

  // For each offset, evaluate the scores
  // what would be the best scores? we will simply take the summation here
  notice("Performing peak-shifting test...");
  std::map<std::string, std::vector<double> > sumScores;
  std::map<std::string, std::vector<int> > numScores;
  for(int i=0; i < nperm; ++i) {
    if ( i % 100 == 0 ) notice("Iteration %d",i+1);
    for(std::map<std::string,genomeLoci>::iterator it = mLoci.begin(); it != mLoci.end(); ++it) {
      double sum = 0;
      int num = 0;
      genomeLoci& loci = it->second;

      std::map<int,double>::iterator jt = pos2score.begin();
      firstpos = newpos = ( jt->first - minpos + offsets[i] ) % length + minpos;
      // move to the position we need
      loci.moveTo(bedchr.c_str(), newpos);
      //notice("++ %s %d %d %d %d %d %d",it->first.c_str(),loci.isend(),loci.it == loci.loci.begin(),loci.loci.size(),loci.it->beg1,firstpos,loci.it->end0);
      while ( jt != pos2score.end() && firstpos <= newpos ) {
        while( !loci.isend() && loci.it->end0 < newpos ) { loci.next(); }
        if ( ( loci.it->beg1 <= newpos ) && ( loci.it->end0 >= newpos ) ) {
          //if ( jt->second > 0.01 ) notice("** %d %d %d %d %.5lf",loci.it->beg1,newpos,loci.it->end0,i, jt->second);
          sum += jt->second;
          ++num;
        }
        ++jt;
        newpos = ( jt->first - minpos + offsets[i] ) % length + minpos;
      }
      loci.rewind();
      while( jt != pos2score.end() ) {
        while( !loci.isend() && loci.it->end0 < newpos ) { loci.next(); }
        if ( ( loci.it->beg1 <= newpos ) && ( loci.it->end0 >= newpos ) ) {
          //if ( jt->second > 0.01 ) notice("** %d %d %d %d %.5lf",loci.it->beg1,newpos,loci.it->end0,i, jt->second);
          sum += jt->second;
          ++num;
        }
        ++jt;
        newpos = ( jt->first - minpos + offsets[i] ) % length + minpos;
      }
      /*
      for(std::map<int,double>::iterator jt = pos2score.begin(); jt != pos2score.end(); ++jt) {
	newpos = ( jt->first - minpos + offsets[i] ) % length + minpos;
	if ( it->second.contains1(bedchr.c_str(), newpos) ) { // ???
	  sum += jt->second;
	  ++num;
	}
      }
      */
      //notice("sum = %lf, num = %d",sum,num);
      sumScores[it->first].push_back(noAvgFlag ? sum : sum/(num+1e-10));
      numScores[it->first].push_back(num);
    }
    //error("offset = %d",offsets[i]);
  }

  notice("Computing permutation-based p-values...");
  wFile wf(out.c_str());  // output file
  wf.printf("#KEY\tFRACBP\tFRACSNP\tSCORE\tSUM\tMIDP\tMAXP\tRANDP\n");
  for(std::map<std::string, std::vector<double> >::iterator it = sumScores.begin(); it != sumScores.end(); ++it) {
    std::vector<double>& s = it->second;
    int countGreater = 0;
    int countSame = 0;
    for(int i=0; i < (int)s.size(); ++i) {
      if ( s[i] == s[0]) ++countSame;
      else if ( s[i] > s[0] ) ++countGreater;
      //printf("%s\t%d\t%.3lg\t%d\t%d\t%d\n",it->first.c_str(), i, s[i], s[i] == s[0], s[i] > s[0], offsets[i]);
    }
    double mpmid = (double)(countGreater + countSame/2.)/(double)nperm;
    double mpmax = (double)(countGreater + countSame)/(double)nperm;
    double mprnd = (double)(countGreater + ((double)rand()/(double)(RAND_MAX+0.5))*countSame)/(double)nperm;
    wf.printf("%s\t%.5lf\t%.5lf\t%.2lg\t%.2lg\t%.5lf\t%.5lf\t%.5lf\n",it->first.c_str(), (double)mLoci[it->first].totalLength()/(double)length,(double)numScores[it->first][0]/(double)pos2score.size(), s[0], noAvgFlag ? s[0] : s[0] * numScores[it->first][0], mpmid, mpmax, mprnd);

    if ( it->first == "5" ) {
      for(int i=1; i < 100; ++i) {
        //wf.printf("--%s\t%.5lf\t%.5lf\t%.2lg\t%.2lg\t%d\t%d\n",it->first.c_str(), (double)mLoci[it->first].totalLength()/(double)length,(double)numScores[it->first][i]/(double)pos2score.size(), s[i], noAvgFlag ? s[i] : s[i] * numScores[it->first][i], s[i] > s[0], s[i] >= s[0]);
      }
    }
  }
  wf.close();
  return 0;
}

// Enrichment test criteria
// 1. How to summarize the distribution of p-values within a region
//   1) z-score, posterior, ..
//   2) accounting for LD
// 2. How to create the null distribution
//   1) contiguity
//   2) number of variants
//   3) size of regions
//   3) number of variants after LD correction
//   4) MAF distribution
//   5) distance to TSS
//   6) genomic features, such as mappability
//   7) background region
/*
int runBedEnrich2(int argc, char** argv) {
  bool ldIgnore = false, ldWeight = true, ldPrune = false;
  bool matchBaseCounts = true, matchVariantCounts = false;
  bool averageScore = true, sumScore = false;
  bool oddsRatio = false, stouffer = true, posterior = false;
  double maxP = 1, minP = 0;
  std::string score;  // EPACTS file containing PVALUE, MAF, LD-adjusted WEIGHTS
  std::string fgBed;  // foreground BED file
  std::string bgBed;  // background BED file
  std::string maskFa; // mask FASTA file
  int nperm = 1000;   // number of permutations
  bool ucsc = false;  
  double minMAF = 1e-10; // minimum MAF;
  double maxMAF = 1;     // maximum MAF
  int nperm;
}
*/

// run enrichment test
int runBedEnrich(int argc, char** argv) {
  // Parse the input arguments
  std::string score;  // EPACTS file containing PVALUE and MAF information
  std::string bed;    // BED file contains [CHROM] [BEG-0] [END-0] [SCORE or CATEGORY]
  int win = 1000000;  // Window to extend from the index bp
  std::string region; // Region to be focused
  std::string chrpos; // chromposomal position of index SNPs
  int maxShift = 1000000;   // maximum shift 
  int minDepth = 1;   // minimum depth
  std::string mask;   // genome mask
  std::string out;    // output files
  std::string outDetail; // detailed output files
  std::string outPost;   // output posterior file
  int rawz = 0;       // use raw z scores instead of converting to posterior probability
  int seed = 0;       // random seed
  int nperm = 1000;   // number of permutations
  bool ucsc = false;
  char buf[65535];
  bool noAvgFlag = false;
  double minMAF = 1e-9;      // minimum MAF
  double maxMAF = 1;         // maximum MAF
  double thresP = 1;         // p-value threshold to consider
  double maxP = 1e-5;        // maximum p-value when there is no significant SNPs around (e.g. 1e-5)
  double minP = 0;           // fix min P
  int binLD = 100000;        // bin for defining LD
  int winLD = 1000000;       // bin for defining maximum LD with causal SNPs

  std::string pvalcol("PVALUE");
  std::string weightcol("WEIGHT");
  std::string mafcol("MAF");
  std::string grp;

  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Input Files")
      LONG_STRINGPARAMETER("score",&score)
      LONG_STRINGPARAMETER("bed",&bed)
      LONG_STRINGPARAMETER("grp",&grp)
      LONG_STRINGPARAMETER("mask",&mask)
      LONG_PARAMETER("ucsc",&ucsc)
      LONG_STRINGPARAMETER("pval-col",&pvalcol)
      LONG_STRINGPARAMETER("maf-col",&mafcol)
      LONG_STRINGPARAMETER("weight-col",&weightcol)

      LONG_PARAMETER_GROUP("Input Options")
      LONG_PARAMETER("no-avg",&noAvgFlag)
      LONG_INTPARAMETER("seed",&seed)
      LONG_INTPARAMETER("nperm",&nperm)
      LONG_DOUBLEPARAMETER("min-maf",&minMAF)
      LONG_DOUBLEPARAMETER("max-maf",&maxMAF)
      LONG_PARAMETER("raw",&rawz)
      LONG_DOUBLEPARAMETER("min-pval",&minP)
      LONG_DOUBLEPARAMETER("thres-pval",&thresP)
      LONG_DOUBLEPARAMETER("max-peak-pval",&maxP)

      LONG_PARAMETER_GROUP("Regions to focus")
      LONG_STRINGPARAMETER("region",&region)
      LONG_STRINGPARAMETER("chrpos",&chrpos)
      LONG_INTPARAMETER("win",&win)
      LONG_INTPARAMETER("max-shift",&maxShift)
      LONG_INTPARAMETER("min-depth",&minDepth)

      LONG_PARAMETER_GROUP("Output Options")
      LONG_STRINGPARAMETER("out",&out)
      LONG_STRINGPARAMETER("out-detail",&outDetail)
      LONG_STRINGPARAMETER("out-post",&outPost)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( score.empty() || out.empty() ) {
    error("--score and --out are required parameters");
  }
  if ( bed.empty() ) {
    error("--bed is required parameters");
  }
  if ( !region.empty() && !chrpos.empty() ) {
    error("Only one of --region or --chrpos can be specified");
  }

  // open the score file first
  std::vector<std::string> tokens;
  if ( !chrpos.empty() ) {
    pFile::tokenizeLine(chrpos.c_str(),":",tokens);
    std::string chr = tokens[0];
    int pos = atoi(tokens[1].c_str());
    int beg = pos - win;
    int end = pos + win;
    sprintf(buf,"%s:%d-%d",chr.c_str(),beg,end);
    region = buf;
  }

  std::map<std::string, std::vector<std::string> > keyGroups;
  if ( ! grp.empty() ) {
    pFile grpf;
    grpf.load(grp.c_str(), NULL, true);
    const char* line = NULL;
    while( ( line = grpf.getLine() ) != NULL ) {
      pFile::tokenizeLine(line," \t\r\n",tokens);
      for(int i=1; i < (int)tokens.size(); ++i) {
        keyGroups[tokens[i]].push_back(tokens[0]);
      }
    }
  }


  // open EPACTS output 
  pFile scoref;
  scoref.load(score.c_str(), region.c_str(), true);

  // read the header line
  const char* line = scoref.getLine();
  if ( line[0] != '#' ) {
    error("Cannot find the header line in %s",score.c_str());
  }

  // determine the columns of PVALUE and MAF
  int ipval = -1;
  int imaf = -1;
  int iweight = -1;
  pFile::tokenizeLine(line," \t\r\n",tokens);

  for(int i=0; i < (int)tokens.size(); ++i) {
    if ( tokens[i].compare(pvalcol) == 0 ) {
      if ( ipval == -1 ) ipval = i;
      else error("Multiple %s column is defined",pvalcol.c_str());
    }
    if ( tokens[i].compare(mafcol) == 0 ) {
      if ( imaf == -1 ) imaf = i;
      else error("Multiple %s column is defined",mafcol.c_str());
    }
    if ( tokens[i].compare(weightcol) == 0 ) {
      if ( iweight == -1 ) iweight = i;
      else error("Multiple %s column is defined",weightcol.c_str());
    }
  }

  if ( ipval < 0 ) error("Cannot find PVALUE column from %s",score.c_str());
  if ( imaf < 0 ) error("Cannot find MAF column from %s",score.c_str());
  if ( iweight < 0 ) warning("WEIGHT column is not found.. Assume 1 for every variant");

  // read the p-value files and convert them into posterior probability scale
  std::map<uint64_t,double> upos2score;
  std::map<uint64_t,double> upos2pvalue;
  std::map<uint64_t,double> upos2weight;
  uint64_t minupos = 18446744073709551615ULL;
  uint64_t maxupos = 0;
  double maf, pval, zscore, weight;
  //double maxz = 0;
  std::map<uint64_t,double> maxzs;  // stores maximum z in LD window
  //if ( minP > 0 ) maxz = hDist::ltqnorm(minP/2.);
  double minPZ = hDist::pval2z(minP);
  double maxPZ = hDist::pval2z(maxP);
  while( ( line = scoref.getLine() ) != NULL ) { // Read EPACTS output file
    pFile::tokenizeLine(line, " \t\r\n", tokens);
    if ( tokens[ipval].compare("NA") != 0 ) {
      uint64_t upos = genomePos.fromPos(tokens[0],atoi(tokens[1].c_str()));
      if ( upos == 0 ) continue;
      if ( minupos > upos ) minupos = upos;
      if ( maxupos < upos ) maxupos = upos;
      maf = atof(tokens[imaf].c_str());
      if ( ( maf >= minMAF ) && ( maf <= maxMAF ) ) {
        pval = atof(tokens[ipval].c_str());
        //zscore = 0-hDist::ltqnorm(pval/2.);
        zscore = hDist::pval2z(pval);

        if ( pval > thresP ) zscore = 0;
        if ( zscore > minPZ ) zscore = minPZ;

        if ( !rawz ) {
          // set max-z for every possible bins
          for( int64_t win = 0-winLD; win <= winLD; win += binLD ) {
            uint64_t key = (upos + win)/binLD;
            if ( maxzs.find(key) == maxzs.end() ) {  // not found
              maxzs[key] = maxPZ; // default maxz is defined by maxP
            }
            if ( maxzs[key] < zscore )
              maxzs[key] = zscore;
          }
        }
        /*
        if ( pval > thresP ) zscore = 0;
        if ( maxz < zscore ) {
          if ( minP == 0 ) {
            maxz = zscore;  // max z is based on min P
          }
          else {
            zscore = maxz;  // max z is a fixed value
          }
        }
        */
        upos2score[upos] = zscore;
        upos2pvalue[upos] = pval;
        if ( iweight >= 0 ) {
          weight = atof(tokens[iweight].c_str());
          upos2weight[upos] = weight;
        }
      }
    }
  }

  wFile* postf = outPost.empty() ? NULL : new wFile(outPost.c_str());
  if ( postf != NULL ) {
    postf->printf("#CHR\tPOS\tSCORE\tWEIGHT\n");
  }

  double psum = 0;
  double maxz = 0;
  if ( !rawz ) {
    for(std::map<uint64_t,double>::iterator it = upos2score.begin();
        it != upos2score.end(); ++it) {
      //it->second = exp(0.5 * (it->second * it->second - maxz * maxz));
      maxz = maxzs[it->first/binLD];
      it->second = exp(0.5 * (it->second * it->second - maxz * maxz));
      psum += (it->second);
    }

    for(std::map<uint64_t,double>::iterator it = upos2score.begin();
        it != upos2score.end(); ++it) {
      it->second /= psum;
      if ( postf != NULL ) {
        std::pair<std::string,int> pos = genomePos.toPos(it->first);
        postf->printf("%s\t%d\t%.4lg\t%.4lf\n",pos.first.c_str(),pos.second,it->second,iweight >= 0 ? upos2weight[it->first] : 1.0);
      }
    }
  }

  // check whether genome-wide or not
  std::pair<std::string,int> minpos = genomePos.toPos(minupos);
  std::pair<std::string,int> maxpos = genomePos.toPos(maxupos);
  uint64_t length = 0;
  if ( minpos.first == maxpos.first ) { // same chromosome
    std::string bedchr = ucsc ? "chr"+minpos.first : minpos.first;
    length = (maxpos.second-minpos.second+1);
    sprintf(buf,"%s:%d-%d",bedchr.c_str(),minpos.second,maxpos.second);
    region = buf;
    notice("Detected single chromosomal region %s",buf);
  }
  else {
    region = "";
    length = genomePos.genomeSize();
    notice("Detected whole genome region spanning multiple chromosomes");
  }

  // read BED file 
  uint64_t firstupos, newupos;
  std::map<std::string, genomePosLoci> mLoci;  // BED regions with predicted category
  notice("Reading BED file...%s on %s",bed.c_str(),region.c_str());
  pFile bedf(bed.c_str(), region.c_str());
  int sz = 0;
  int nregions = 0;
  int ngroups = 0;
  while( (line = bedf.getLine()) != NULL ) {
    pFile::tokenizeLine(line, " \t\r\n", tokens);
    std::string key("ALL");
    if ( tokens.size() < 3 ) {
      error("BED file must have at least three columns");
    }
    else if ( tokens.size() > 3 ) {
      key = tokens[3];
    }
    std::string chr = ucsc ? tokens[0].substr(3) : tokens[0];
    uint64_t ubeg = genomePos.fromPos(chr,atoi(tokens[1].c_str())+1);
    uint64_t uend = genomePos.fromPos(chr,atoi(tokens[2].c_str()));
    if ( ( ubeg == 0 ) || ( uend == 0 ) ) continue;
    if ( ubeg < minupos ) ubeg = minupos;
    if ( uend > maxupos ) uend = maxupos;
    if ( ubeg > uend ) {
      continue;
      //error("Invalid interval %s",line);
    }
    sz += (uend-ubeg+1);

    if ( grp.empty() ) {
      mLoci[key].add(ubeg, uend);
      ++ngroups;
    }
    else {
      std::vector<std::string>& v = keyGroups[key];
      for(int i=0; i < (int)v.size(); ++i) {
        mLoci[v[i]].add(ubeg, uend);
        ++ngroups;
      }
    }
    ++nregions;
  }
  notice("Identified %d regions over %d groups of %llu base pairs across %d categories", nregions, ngroups, length,  mLoci.size());

  notice("Resolving Overlaps...");
  for(std::map<std::string,genomePosLoci>::iterator it = mLoci.begin(); it != mLoci.end(); ++it) {
    //notice("%s\t%llu",it->first.c_str(),it->second.totalLength());
    it->second.resolveOverlaps();
    //notice("%s\t%llu",it->first.c_str(),it->second.totalLength());
    //error("%s\t%llu",it->first.c_str(),it->second.totalLength());
  }

  if ( maxShift == 0 ) {
    maxShift = (maxupos-minupos+1);
  }

  notice("Generating random offsets... with max-shift=%d\n", maxShift);
  if ( seed == 0 ) srand(std::time(0));
  else srand(seed);
  std::vector<int> offsets;
  offsets.push_back(0);
  for(int i=1; i < nperm; ++i) {
    offsets.push_back( rand() % maxShift - maxShift / 2 );
  }

  // For each offset, evaluate the scores
  // what would be the best scores? we will simply take the summation here
  notice("Performing peak-shifting test...");
  std::map<std::string, std::vector<double> > sumScores;
  std::map<std::string, std::vector<int> > numScores;
  std::map<std::string, double> minPvalue;
  for(int i=0; i < nperm; ++i) {
    if ( i % 100 == 0 ) notice("Iteration %d",i+1);
    for(std::map<std::string,genomePosLoci>::iterator it = mLoci.begin(); it != mLoci.end(); ++it) {
      double minpval = 1;
      double sum = 0;
      double num = 0;
      genomePosLoci& loci = it->second;

      std::map<uint64_t,double>::iterator jt = upos2score.begin();
      firstupos = newupos = ( jt->first - minupos + offsets[i] ) % length + minupos;
      // move to the position we need
      loci.moveTo(newupos);
      while ( jt != upos2score.end() && firstupos <= newupos ) {
        while( !loci.isend() && loci.it->uend0 < newupos ) { loci.next(); }
        if ( ( loci.it->ubeg1 <= newupos ) && ( loci.it->uend0 >= newupos ) ) {
          sum += (iweight < 0) ? jt->second : (jt->second * upos2weight[jt->first]);
          num += (iweight < 0) ? 1.0 : (upos2weight[jt->first]);
          if ( i == 0 ) {
            double p = upos2pvalue[jt->first];
            if ( p < minpval ) minpval = p;
          }
        }
        ++jt;
        newupos = ( jt->first - minupos + offsets[i] ) % length + minupos;
      }
      loci.rewind();
      while( jt != upos2score.end() ) {
        while( !loci.isend() && loci.it->uend0 < newupos ) { loci.next(); }
        if ( ( loci.it->ubeg1 <= newupos ) && ( loci.it->uend0 >= newupos ) ) {
          sum += (iweight < 0) ? jt->second : (jt->second * upos2weight[jt->first]);
          num += (iweight < 0) ? 1.0 : (upos2weight[jt->first]);
          if ( i == 0 ) {
            double p = upos2pvalue[jt->first];
            if ( p < minpval ) minpval = p;
          }
        }
        ++jt;
        newupos = ( jt->first - minupos + offsets[i] ) % length + minupos;
      }
      sumScores[it->first].push_back(noAvgFlag ? sum : sum/(num+1e-10));
      numScores[it->first].push_back(num);
      if ( i == 0 ) {
        minPvalue[it->first] = minpval;
      }
    }
  }

  notice("Computing permutation-based p-values...");
  wFile wf(out.c_str());  // output file
  wf.printf("#KEY\tFRACBP\tFRACSNP\tMINPVAL\tSCORE\tSUM\tMIDP\tMAXP\tRANDP\n");
  wFile* pwdf = outDetail.empty() ? NULL : new wFile(outDetail.c_str());

  for(std::map<std::string, std::vector<double> >::iterator it = sumScores.begin(); it != sumScores.end(); ++it) {
    std::vector<double>& s = it->second;
    int countGreater = 0;
    int countSame = 0;

    if ( pwdf ) pwdf->printf("%s",it->first.c_str());

    for(int i=0; i < (int)s.size(); ++i) {
      if ( s[i] == s[0]) ++countSame;
      else if ( s[i] > s[0] ) ++countGreater;

      if ( pwdf ) pwdf->printf("\t%.6lg",s[i]);
    }

    if ( pwdf ) pwdf->printf("\n");

    double mpmid = (double)(countGreater + countSame/2.)/(double)nperm;
    double mpmax = (double)(countGreater + countSame)/(double)nperm;
    double mprnd = (double)(countGreater + ((double)rand()/(double)(RAND_MAX+0.5))*countSame)/(double)nperm;
    wf.printf("%s\t%.5lf\t%.5lf\t%.2lg\t%.2lg\t%.2lg\t%.5lf\t%.5lf\t%.5lf\n",it->first.c_str(), (double)mLoci[it->first].totalLength()/(double)length,(double)numScores[it->first][0]/(double)upos2score.size(), minPvalue[it->first], s[0], noAvgFlag ? s[0] : s[0] * numScores[it->first][0], mpmid, mpmax, mprnd);
  }
  wf.close();
  if ( pwdf ) {
    pwdf->close();
    delete pwdf;
  }
  return 0;
}

int runBGZIndex(int argc, char** argv) {
  ParameterList pl;
  std::string gzf;
  std::string outf("-");
  std::string sbyte("0");
  int nelems = 1;

  bool indexFlag = false;
  bool readFlag = true;
  bool forceFlag = false;

  bool uint8_t_Flag = true;
  bool uint16_t_Flag = false;
  bool uint32_t_Flag = false;
  bool uint64_t_Flag = false;
  bool int8_t_Flag = false;
  bool int16_t_Flag = false;
  bool int32_t_Flag = false;
  bool int64_t_Flag = false;
  bool char_t_Flag = false;
  bool float_t_Flag = false;
  bool double_t_Flag = false;

  bool verbose = false;

  BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Input Options")
      LONG_STRINGPARAMETER("in",&gzf)
      LONG_STRINGPARAMETER("offset",&sbyte)
      LONG_INTPARAMETER("nelems",&nelems)

      LONG_PARAMETER_GROUP("Read or Index")
      EXCLUSIVE_PARAMETER("index",&indexFlag)
      EXCLUSIVE_PARAMETER("read",&readFlag)

      LONG_PARAMETER_GROUP("Read Type")
      EXCLUSIVE_PARAMETER("uint8",&uint8_t_Flag)
      EXCLUSIVE_PARAMETER("uint16",&uint16_t_Flag)
      EXCLUSIVE_PARAMETER("uint32",&uint32_t_Flag)
      EXCLUSIVE_PARAMETER("uint64",&uint64_t_Flag)
      EXCLUSIVE_PARAMETER("int8",&int8_t_Flag)
      EXCLUSIVE_PARAMETER("int16",&int16_t_Flag)
      EXCLUSIVE_PARAMETER("int32",&int32_t_Flag)
      EXCLUSIVE_PARAMETER("int64",&int64_t_Flag)
      EXCLUSIVE_PARAMETER("char",&char_t_Flag)
      EXCLUSIVE_PARAMETER("float",&float_t_Flag)
      EXCLUSIVE_PARAMETER("double",&double_t_Flag)

      LONG_PARAMETER_GROUP("Output Options")
      LONG_PARAMETER("force",&forceFlag)
      LONG_PARAMETER("verbose",&verbose)
      LONG_STRINGPARAMETER("out",&outf)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc, argv);
  if ( verbose ) pl.Status();

  // sanity check of input arguments
  if ( gzf.empty() ) {
    error("--in are required parameters");
  }
  if ( ( forceFlag ) && ( readFlag ) ) {
    error("--force option is unavailable when --read is on");
  }

  if ( indexFlag ) {
    FILE* fp = fopen(gzf.c_str(), "r");
    if ( fp == NULL ) error("Error in reading BGZF file %s",gzf.c_str());

    std::string idxf = gzf + ".bgi";
    if ( !forceFlag ) {
      FILE* ftmp = fopen(idxf.c_str(), "r");
      if ( ftmp != NULL ) error("BGZF Index %s already exists. Use --force option to overwrite",idxf.c_str());
    }
    FILE* wp = fopen(idxf.c_str(), "w");
    if ( wp == NULL ) error("Cannot write to %s",idxf.c_str());

    uint8_t u8;
    uint16_t u16;
    uint32_t u32;

    if ( ( fread(&u8, sizeof(uint8_t), 1, fp) != 1 ) || ( u8 != 31 ) ) error("Error in reading BGZF file %s",gzf.c_str());
    if ( ( fread(&u8, sizeof(uint8_t), 1, fp) != 1 ) || ( u8 != 139 ) ) error("Error in reading BGZF file %s",gzf.c_str());
    if ( ( fread(&u8, sizeof(uint8_t), 1, fp) != 1 ) || ( u8 != 8 ) ) error("Error in reading BGZF file %s",gzf.c_str());
    if ( ( fread(&u8, sizeof(uint8_t), 1, fp) != 1 ) || ( u8 != 4 ) ) error("Error in reading BGZF file %s",gzf.c_str());
    if ( ( fread(&u32, sizeof(uint32_t), 1, fp) != 1 ) ) error("Error in reading BGZF file %s",gzf.c_str());
    if ( ( fread(&u8, sizeof(uint8_t), 1, fp) != 1 ) ) error("Error in reading BGZF file %s",gzf.c_str());
    if ( ( fread(&u8, sizeof(uint8_t), 1, fp) != 1 ) ) error("Error in reading BGZF file %s",gzf.c_str());
    if ( ( fread(&u16, sizeof(uint16_t), 1, fp) != 1 ) ) error("Error in reading BGZF file %s",gzf.c_str());
    //if ( u16 > 0 ) fseek(fp, u16, SEEK_CUR);

    uint64_t sz = 0;
    for(int i=0;; ++i) {
      if ( fwrite(&sz, sizeof(uint64_t), 1, wp) != 1 ) error("Error in writing BGZF Index %s",idxf.c_str());

      if ( fread(&u8, sizeof(uint8_t), 1, fp) != 1 ) break;
      else if ( u8 != 66 ) error("Error1 in reading BGZF file %s at block %d %u",gzf.c_str(),i,u8);
      if ( ( fread(&u8, sizeof(uint8_t), 1, fp) != 1 ) || ( u8 != 67 ) ) error("Error3 in reading BGZF file %s at block %d",gzf.c_str(),i);
      if ( ( fread(&u16, sizeof(uint16_t), 1, fp) != 1 ) || ( u16 != 2 ) ) error("Error4 in reading BGZF file %s at block %d",gzf.c_str(),i);
      if ( ( fread(&u16, sizeof(uint16_t), 1, fp) != 1 ) ) error("Error5 in reading BGZF file %s at block %d",gzf.c_str(),i);
      sz += (u16+1);
      //fprintf(stderr,"sz=%lu\n",sz);
      if ( fseek(fp, u16-5, SEEK_CUR) != 0 ) break;
    }
    fclose(wp);
    fclose(fp);
  }
  else {
    uint64_t byte = strtoul(sbyte.c_str(),NULL,0);
    uint64_t block = (byte >> 16);
    uint32_t offset = (byte & 0x0FFFF);

    uint64_t cblock = 0;
    if ( block > 0 ) {
      std::string idxf = gzf + ".bgi";
      FILE* fp = fopen(idxf.c_str(), "r");
      if ( fp == NULL ) error("Cannot open %s",idxf.c_str());
      fseek(fp, block * sizeof(uint64_t), SEEK_SET);
      if ( fread(&cblock, sizeof(uint64_t), 1, fp) != 1 )
        error("Cannot read from %s",idxf.c_str());
      fclose(fp);
    }

    int64_t start = (cblock << 16 | offset);
    int32_t szelem = ( uint8_t_Flag | int8_t_Flag | char_t_Flag ) ? sizeof(uint8_t) :
                     ( (uint16_t_Flag | int16_t_Flag ) ? sizeof(uint16_t) :
                       ( (uint32_t_Flag | int32_t_Flag | float_t_Flag) ? sizeof(uint32_t) :
                         ( ( uint64_t_Flag | int64_t_Flag | double_t_Flag ) ? sizeof(uint64_t) : -1 ) ) );
    if ( szelem == -1 ) {
      uint8_t_Flag = true;
      szelem = sizeof(uint8_t);
    }
    int64_t end = start + nelems * szelem;

    char* buffer = new char[WINDOW_SIZE];
    BGZF* bp = bgzf_open(gzf.c_str(),"r");

    wFile wf(outf.c_str());
    int c;

    if (bgzf_seek(bp, start, SEEK_SET) < 0) error("Cannot perform bgzf_seek in %s to %ld",gzf.c_str(),start);
    while (1) {
      if (end < 0) c = bgzf_read(bp, buffer, WINDOW_SIZE);
      else c = bgzf_read(bp, buffer, (end - start > WINDOW_SIZE)? WINDOW_SIZE:(end - start));
      if (c == 0) break;
      if (c < 0) error("Cannot perform bgzf_seek in %s",gzf.c_str());
      start += c;
      for(int i=0; i < c; i += szelem) {
        if ( uint8_t_Flag ) {
          wf.printf("%u\n",*(const uint8_t*)(buffer+i));
        }
        else if ( int8_t_Flag ) {
          wf.printf("%d\n",*(const int8_t*)(buffer+i));
        }
        else if ( char_t_Flag ) {
          wf.printf("%c\n",*(buffer+i));
        }
        else if ( uint16_t_Flag ) {
          wf.printf("%u\n",*(const uint16_t*)(buffer+i));
        }
        else if ( int16_t_Flag ) {
          wf.printf("%d\n",*(const int16_t*)(buffer+i));
        }
        else if ( uint32_t_Flag ) {
          wf.printf("%lu\n",*(const uint16_t*)(buffer+i));
        }
        else if ( int32_t_Flag ) {
          wf.printf("%ld\n",*(const int16_t*)(buffer+i));
        }
        else if ( float_t_Flag ) {
          wf.printf("%g\n",*(const float*)(buffer+i));
        }
        else if ( uint64_t_Flag ) {
          wf.printf("%llu\n",*(const uint64_t*)(buffer+i));
        }
        else if ( int32_t_Flag ) {
          wf.printf("%lld\n",*(const int64_t*)(buffer+i));
        }
        else if ( double_t_Flag ) {
          wf.printf("%g\n",*(const double*)(buffer+i));
        }
        else {
          error("Unrecognized output type");
        }
      }
      //write(f_dst, buffer, c);
      if (end >= 0 && start >= end) break;
    }
    delete[] buffer;
    wf.close();
  }
  return 0;
}


int main(int argc, char** argv) {
  if ( argc < 2 ) {
    printf("vcfast -- Fast analytic tools for analyzing and manipulating VCF\n");
    printf("Copyright (c) 2012 Hyun Min Kang\n");
    printf("Usage : %s [gen-hash] [options]\n",argv[0]);
    printf("\tType one of the following commands below to get detailed usage\n");
    printf("\t%s summary      [options] : summarize vcf\n",argv[0]);
    printf("\t%s convert      [options] : convert to different format\n",argv[0]);
    printf("\t%s index-LD     [options] : compute LD with index SNP\n",argv[0]);
    printf("\t%s pair-LD      [options] : compute pairwise LD\n",argv[0]);
  }
  else {
    std::string cmd(argv[1]);
    if ( cmd == "summary" ) {
      return runSummary(argc-1,argv+1);
    }
    else if ( cmd == "convert" ) {
      return runConvert(argc-1,argv+1);
    }
    else if ( cmd == "index-LD" ) {
      return runIndexLD(argc-1,argv+1);
    }
    else if ( cmd == "pair-LD" ) {
      return runPairLD(argc-1,argv+1);
    }
    else if ( cmd == "peak-shift" ) {
      return runPeakShift(argc-1,argv+1);
    }
    else if ( cmd == "bed-enrich" ) {
      return runBedEnrich(argc-1,argv+1);
    }
    else if ( cmd == "bgzf-index" ) {
      return runBGZIndex(argc-1,argv+1);
    }
    else {
      error("Unrecognized command %s\n",argv[0]);
    }
  }
  return 0;
}
