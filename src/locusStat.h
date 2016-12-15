#include <cfloat>

class locusStat {
 protected:
  double nBase;
  double nPass;
  int nInds;
  double sumScores;
  double sqScores;

  std::vector<int> etas;

 public:
 locusStat(int n) : nBase(0), nBase(0), nInds(n) { 
    etas.resize(nInds+nInds+1,0); 
  }

  void add(double fracBase, bool pass, int alleleCount, double score) {
    nBase += fracBase;
    if ( pass ) {
      nPass += fracBase;
      sumScores += fracBase*(score);
      sqScores += fracBase*(score*score);
      if ( ( alleleCount < 0 ) || ( alleleCount > 2*$nInds ) ) {
	error("Allele count %d is out of range for %d individuals",alleleCount,nInds);
      }
      ++(etas[alleleCount]);
    }
  }

  std::pair<double,double> avgScore() {
    double mu = 0, sigma = DBL_MAX;
    if ( nPass > 0 ) {
      mu = sumScores / nPass;
      sigma = sqrt(sqScores / nPass - mu * mu);;
    }
    return std::pair<double,double>(mu,sigma);
  }

  std::pair<double,double> snpDensity() {
    double mu = 0, sigma = DBL_MAX;
    if ( nPass > 0 ) {
    }
  }
  
  //std::vector<float> scores; // for every bg/fg bases
  //std::vector<int> acs;      // for each fg/bg bases
}
