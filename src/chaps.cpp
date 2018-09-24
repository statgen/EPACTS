#include "Parameters.h"
#include "Error.h"
#include "fVcf.h"
#include "pcHMM.h"
#include "wFile.h"
#include "icGraph.h"
#include "sparseHMM.h"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <getopt.h>
#include <ctime>

class pCHAPSArgs {
public:
  // VCF-related string arguments
  std::string vcf;
  std::string inf;
  std::string region;
  std::string rule;
  std::string field;

  // Other input files
  std::string indf;
  std::string outf;

  int unit;
  bool verbose;
  bool ignoreFilter;
  bool ignoreMissing;

  int minAC;
  int minMAC;
  int maxAC;
  double minMAF;
  double maxMAF;
  double minCallRate;
  bool sepchr;

  static int const DEFAULT_NSNPS = 10000L;
  constexpr static double const DEFAULT_MIN_MAF = 1e-6;
  constexpr static double const DEFAULT_MAX_MAF = 1;

  pCHAPSArgs() :
    unit(DEFAULT_NSNPS), verbose(false), ignoreFilter(false), ignoreMissing(false), minAC(0), minMAC(0), maxAC(INT_MAX), minMAF(DEFAULT_MIN_MAF), maxMAF(DEFAULT_MAX_MAF), minCallRate(DEFAULT_MIN_MAF), sepchr(false)
    {}
};

static bool writeIDs(wFile& wf, std::vector<std::string>& ids) {
  char* p = (char*) malloc(65536); // maximum ID length is 64k
  
  // write individual IDs to file
  for(int i=0; i < (int)ids.size(); ++i) {
    *(int*)p = (int)ids[i].size();
    ids[i].copy(p+sizeof(int),(int)ids[i].size());
    wf.write(p,(int)ids[i].size() + sizeof(int));
  }
  free(p);
  return true;
}

void writeChapsHeader(wFile& wf, std::vector<std::string>& markers, std::vector<double>& AFs, std::vector<std::string>& inds, uint8_t* PLs, int maxM) {
    const char* MAGIC = "CHAPS_HMM";
    int N_MAGIC = strlen(MAGIC);
    int m = (int)markers.size();
    int n = (int)inds.size();
    int i;

    wf.write((void*)MAGIC, sizeof(char) * N_MAGIC);
    wf.write(&n, sizeof(int));
    wf.write(&m, sizeof(int));
    writeIDs(wf, inds);
    writeIDs(wf, markers);
    wf.write(AFs.data(), sizeof(double) * m);
    for(i=0; i < n; ++i) {  // write PLs
      wf.write(&PLs[maxM * 3 * i], m * 3);
    }
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

// run pairwise HMM to estimate the number of shared chromosomes
int runPairHMM(int argc, char** argv) {
  pCHAPSArgs arg;
  arg.field = "PL";
  double mu = 1e-8;
  double theta = 1e-8; // .01 per Mb
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
    LONG_INTPARAMETER("minMAC",&arg.minMAC)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)
    LONG_PARAMETER("sepchr",&arg.sepchr)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_DOUBLEPARAMETER("mu",&mu)
    LONG_DOUBLEPARAMETER("theta",&theta)

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

  fVcf tvcf;
  tvcf.load(arg.vcf.c_str(), arg.region.c_str(), arg.field.c_str(), arg.rule.c_str(), !arg.ignoreFilter, arg.indf.empty() ? NULL : arg.indf.c_str());
  int n = tvcf.nInds;
  int m = 0;
  int maxM = 1000;
  std::vector<double> AFs;
  std::vector<std::string> markers;
  std::vector<int> pos1s;
  uint8_t* PLs = new uint8_t[n * maxM * 3];
  double af, maf;

  for(int M = 0; tvcf.readMarkers(1); ++M) {
    if ( M % 10000 == 0 )
      notice("Processing %d/%d markers across %d individuals..\n", m, M, tvcf.nInds);
    af = tvcf.alleleFreq(0);
    maf = af > 0.5 ? 1-af : af;
    if ( maf < 0 ) maf = 0;
    if ( ( maf >= arg.minMAF ) && ( tvcf.callRate(0) >= arg.minCallRate ) && ( tvcf.sumAlleles[0] >= arg.minAC ) && ( tvcf.sumAlleles[0] <= arg.maxAC ) ) {
      // see whether we reached maxM
      if ( maxM == m ) {
        notice("Expanding maxM from %d to %d", maxM, maxM*2);
        uint8_t* cPLs = PLs;
        PLs = new uint8_t[n * maxM * 6];
        for(int i=0; i < n; ++i) {
          memcpy( &PLs[i*maxM*6], &cPLs[i*maxM*3], sizeof(uint8_t) * (maxM *3) );
        }
        maxM *= 2;
        delete [] cPLs;
      }
      pos1s.push_back(tvcf.pos1s[0]);
      AFs.push_back(af);
      markers.push_back(tvcf.markers[0]);
      for(int i=0; i < n; ++i) {
        PLs[i * maxM * 3 + m * 3 + 0] = tvcf.PLs[i*3];
        PLs[i * maxM * 3 + m * 3 + 1] = tvcf.PLs[i*3+1];
        PLs[i * maxM * 3 + m * 3 + 2] = tvcf.PLs[i*3+2];
      }
      ++m;
    }
  }
  notice("Finished loading the genotype likelihoods across %d markers",m);

  pcHMM pchmm(m, pos1s, AFs, mu, theta, .05, .001);
  wFile wf(arg.outf.c_str());

  writeChapsHeader(wf, markers, AFs, tvcf.inds, PLs, maxM);
  for(int i=0; i < n; ++i) {
    notice("Processing %s - index %d",tvcf.inds[i].c_str(),i+1);
    for(int j=0; j < i; ++j) {
      pchmm.viterbi( &PLs[j*maxM*3], &PLs[i*maxM*3] );
      pchmm.writeViterbi(wf);
      //wf.printf("%s\t%s",tvcf.inds[j].c_str(),tvcf.inds[i].c_str());
      //pchmm.writeViterbi(wf, markers, tvcf.inds);
    }
    //if ( i == 2 ) break;
  }
  delete [] PLs;
  wf.close();
  return 0;
}

int runBuildGraph(int argc, char** argv) {
  pCHAPSArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Options")
    LONG_STRINGPARAMETER("in",&arg.inf)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out",&arg.outf)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.inf.empty() || arg.outf.empty() ) {
    error("--in and --out required parameter");
  }

  icGraph graph;
  graph.loadPairHMM(arg.inf.c_str());

  wFile wf(arg.outf.c_str());
  graph.print(wf);

  return 0;
}

int runHaplotyper(int argc, char** argv) {
  pCHAPSArgs arg;
  int minState = 200;
  int maxState = 400;
  int rndState = 20;
  int seed = 0;
  int niter = 20;
  bool dosage = false;
  bool gprob = false;
  ParameterList pl;
  
  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Options")
    LONG_STRINGPARAMETER("in",&arg.inf)
    LONG_INTPARAMETER("min-state",&minState)
    LONG_INTPARAMETER("max-state",&maxState)
    LONG_INTPARAMETER("rand-state",&rndState)
    LONG_INTPARAMETER("seed",&seed)
    LONG_INTPARAMETER("niter",&niter)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out",&arg.outf)
    LONG_PARAMETER("dosage",&dosage)
    LONG_PARAMETER("gprob",&gprob)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.inf.empty() || arg.outf.empty() ) {
    error("--in and --out required parameter");
  }

  // load the graph structure into the HMM
  sparseHMM shmm(arg.inf.c_str(), minState, maxState, rndState, seed);

  // initialize Haplotypes
  shmm.initializeHaplotypes(true);

  notice("foo");

  for(int i=0; i < shmm.graph.n; ++i) {
    fprintf(stderr,"Running forward-backward algorithm for %d",i);
    shmm.betweenForwardBackward(i);
  }

  shmm.writeHaplotypes(arg.outf.c_str());

  return 0;
}

int main(int argc, char** argv) {
  if ( argc < 2 ) {
    printf("chaps -- Haplotype Reconstruction via Compatible Haplotype Graph VCF\n");
    printf("Copyright (c) 2012 Hyun Min Kang\n");
    printf("Usage : %s [gen-hash] [options]\n",argv[0]);
    printf("\tType one of the following commands below to get detailed usage\n");
    printf("\t%s pair-hmm     [options] : run pairwise hmm\n",argv[0]);
  }
  else {
    std::string cmd(argv[1]);
    if ( cmd == "pair-hmm" ) {
      return runPairHMM(argc-1,argv+1);
    }
    else if ( cmd == "build-graph") {
      return runBuildGraph(argc-1,argv+1);
    }
    else if ( cmd == "haplotype" ) {
      return runHaplotyper(argc-1,argv+1);
    }
    else {
      error("Unrecognized command %s\n",argv[0]);
    }
  }
  return 0;
}
