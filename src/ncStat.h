#ifndef __NC_STAT_H__
#define __NC_STAT_H__

#include <cfloat>
#include <cmath>
#include <climits>

#include "codonHelper.h"
#include "Error.h"

class ncKey {
public:
  int scorebin;
  uint8_t cpg;
  uint8_t ref;
  uint16_t kmer;

  ncKey(int _scorebin, int _cpg, int _ref, int _kmer) {
    scorebin = _scorebin;
    cpg = (_cpg == -1) ? 255 : (uint8_t)_cpg;
    ref = (_ref == -1) ? 255 : (uint8_t)_ref;
    kmer = (_kmer == -1) ? 65535 : (uint16_t)_kmer;
  }

  bool operator< (const ncKey& x) const {
    if ( scorebin == x.scorebin ) {
      if ( cpg == x.cpg ) {
	if ( ref == x.ref ) {
	  return (kmer < x.kmer);
	}
	else return (ref < x.ref);
      }
      else return (cpg < x.cpg);
    }
    else return (scorebin < x.scorebin);
  }
};

class ncStat {
 public:
  static const double NAN_DBL;

  int nBase;
  int nPass;
  int nInds;
  double sumCons;
  double sqCons;
  std::vector<int> etas;
  std::vector<int> cnts;
  std::vector<int> coms;

  ncStat() : nBase(0), nPass(0), nInds(0), sumCons(0), sqCons(0) {}
  void setN(int n) { 
    nInds = n; 
    etas.clear();
    etas.resize(8*n,0); 
    cnts.resize(4,0); 
    coms.resize(4,0); 
  }

  void addBase(bool pass, float score) {
    ++nBase;
    if ( pass ) { 
      ++nPass;
      sumCons += score;
      sqCons += (score*score);
    }
  }

  void addVariant(uint8_t altAllele, int altCount) {
    if ( altCount > 0 ) {
      if ( ( altAllele > 0 ) && ( altAllele < 5 ) ) {
	++etas[(altAllele-1) * (2 * nInds) + (altCount-1)];
	if ( altCount < 2*nInds ) ++cnts[(altAllele-1)];
	if ( ( altCount >= 0.005*2*nInds ) && ( altCount <= 0.995*2*nInds ) ) ++coms[(altAllele-1)]; // AF>.5%
      }
      else {
	error("altAllele must be between 1 and 4");
      }
    }
  }

  // SE is estimated with one pseudo-count
  // mean = n/N
  // SD = sqrt((n+1)/(N+1)*((N-n)/(N+1))/(N+1))
  std::pair<double,double> allDensity(int altAllele = 0) {
    if ( nPass == 0 ) { return std::pair<double,double> (0,NAN_DBL); }
    else { 
      int n = (altAllele == 0) ? (cnts[0]+cnts[1]+cnts[2]+cnts[3]) : cnts[(altAllele-1)];
      double p = (double)n / (double)nPass;
      double sd = sqrt((n+1.) * (nPass-n) / (nPass+1.))/(nPass+1.);
      return std::pair<double,double> (p,sd);
    }
  }

  // SE is estimated with one pseudo-count
  // mean = n/N
  // SD = sqrt((n+1)/(N+1)*((N-n)/(N+1))/(N+1))
  std::pair<double,double> comDensity(int altAllele = 0) {
    if ( nPass == 0 ) { return std::pair<double,double> (0,NAN_DBL); }
    else { 
      int n = (altAllele == 0) ? (coms[0]+coms[1]+coms[2]+coms[3]) : coms[(altAllele-1)];
      double p = (double)n / (double)nPass;
      double sd = sqrt((n+1.) * (nPass-n) / (nPass+1.))/(nPass+1.);
      return std::pair<double,double> (p,sd);
    }
  }

  std::pair<double,double> avgCons() {
    if ( nPass == 0 ) { return std::pair<double,double> (0,NAN_DBL); }
    else { return std::pair<double,double> ( sumCons / nPass, sqrt( (sqCons/nPass - sumCons*sumCons/nPass/nPass) / nPass ) ); }
  }

  std::pair<double,double> theta(int altAllele = 0) {
    if ( nPass == 0 ) { return std::pair<double,double> (0,NAN_DBL); }
    else {
      double sum = 0;
      //double altMono = etas[0*2*nInds+2*nInds-1]+etas[1*2*nInds+2*nInds-1]+etas[2*2*nInds+2*nInds-1]+etas[3*2*nInds+2*nInds-1];
      for(int i=1; i < nInds*2; ++i) {
	sum += ((double)i*(nInds+nInds-i)*(double)(altAllele == 0 ? (etas[0*(2*nInds)+i-1]+etas[1*(2*nInds)+i-1]+etas[2*(2*nInds)+i-1]+etas[3*(2*nInds)+i-1]) : etas[(altAllele-1)*(2*nInds)+i-1]));
      }
      //notice("sum=%lf\n",sum);
      double theta = (sum == 0) ? ((nInds+nInds-1.)*(1./(2.0*nInds*(2.0*nInds-1.0)))/(double)nPass) : (sum * (2./(2.0*nInds*(2.0*nInds-1.0))/(double)nPass) );
      double seTheta = sqrt( (2.0*nInds+1.0)/(6.0*nInds-3.0)*theta + 2.0*(4.0*nInds*nInds+2.0*nInds+3.)/(18.0*nInds*(2.0*nInds-1.0))*theta*theta );
      return std::pair<double,double> (theta,seTheta);
    }
  }
};

const double ncStat::NAN_DBL = sqrt(-1.);  // assign double NAN value

#endif
