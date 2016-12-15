#ifndef __PC_HMM__H
#define __PC_HMM__H

#include <vector>
#include "Error.h"
#include "PhredHelper.h"
#include "wFile.h"

class pcHMM {
 public:
  int m;
  std::vector<double> trans;
  std::vector<double> afs;
  std::vector<double> deltas;
  std::vector<int> phis;
  std::vector<int> paths;
  std::vector<double> emis;
  double mu;
  double theta;
  double pi1, pi2;

 pcHMM(int _m, std::vector<int>& _pos1s, std::vector<double>& _afs, double _mu, double _theta, double _pi1, double _pi2) : m(_m), afs(_afs), mu(_mu), theta(_theta), pi1(_pi1), pi2(_pi2) {
    if ( (int)_pos1s.size() != m ) error("pos1s.size() = %d != %d",(int)_pos1s.size(),m);
    if ( (int)afs.size() != m ) error("afs.size() = %d != %d",(int)afs.size(),m);
    if ( pi1 + pi2 >= 1 ) error("pi1 + pi2 = %lf",pi1+pi2);
    for(int i=0; i < m-1; ++i) {
      trans.push_back(log((_pos1s[i+1]-_pos1s[i]+1e-2)*theta));
    }
    deltas.resize(3*m);
    phis.resize(3*m);
    paths.resize(m);
    emis.resize(3*m);
  }

  void viterbi(uint8_t* PL1, uint8_t* PL2) {
    // calculate emis
    int i, j, h0, h1, h2, h3;
    uint8_t *pl1, *pl2;
    double p20, p21, p22, p1, tmp, l20, l21, l22, AF;
    pl1 = PL1; pl2 = PL2;

    for(i=0; i < m; ++i) {
      AF = afs[i];
      l20 = l21 = l22 = 0;
      for(h0 = 0; h0 < 2; ++h0) {
	for(h1 = 0; h1 < 2; ++h1) {
	  p1 = ((h0 ? AF : 1-AF) * (h1 ? AF : 1-AF)) * phredConv.phred2Err[pl1[h0+h1]];
	  for(h2 = 0; h2 < 2; ++h2) {
	    for(h3 = 0; h3 < 2; ++h3) {
	      tmp = p1 * phredConv.phred2Err[pl2[h2+h3]];
	      p20 = (h2 ? AF : 1-AF) * (h3 ? AF : 1-AF);
	      p21 = (h2 ? AF : 1-AF) * (h3 == h1 ? 1-mu : mu);
	      p22 = (h2 == h0 ? 1-mu : mu) * (h3 == h1 ? 1-mu : mu);
	      l20 += (tmp * p20);
	      l21 += (tmp * p21);
	      l22 += (tmp * p22);
	    }
	  }
	}
      }
      emis[3*i+0] = log(l20);
      emis[3*i+1] = log(l21);
      emis[3*i+2] = log(l22);
      pl1 += 3;
      pl2 += 3;
      //printf("%d\t%lf\t%lf\t%lf\t%lf\t%d,%d,%d\t%d,%d,%d\n",i,l20,l21,l22,AF,pl1[0],pl1[1],pl1[2],pl2[0],pl2[1],pl2[2]);
      //if (i == 100) error("foo");
    }

    // calculate delta
    deltas[0] = log(1-pi1-pi2) + emis[0];
    deltas[1] = log(pi1) + emis[1];
    deltas[2] = log(pi2) + emis[2];

    for(i=1; i < m; ++i) {
      for(j=0; j < 3; ++j) {
	p20 = deltas[(i-1)*3+0] + (j == 0 ? 0 : ( j == 1 ? trans[i-1] : trans[i-1]*2 ) );
	p21 = deltas[(i-1)*3+1] + (j == 0 ? trans[i-1] : ( j == 1 ? 0 : trans[i-1] ) );
	p22 = deltas[(i-1)*3+2] + (j == 0 ? trans[i-1] * 2 : ( j == 1 ? trans[i-1] : 0 ) );
	if ( p20 > p21 ) {
	  if ( p20 > p22 ) {
	    deltas[i*3+j] = p20 + emis[3*i+j];
	    phis[i*3+j] = 0;
	  }
	  else {
	    deltas[i*3+j] = p22 + emis[3*i+j];
	    phis[i*3+j] = 2;
	  }
	}
	else if ( p21 > p22 ) {
	  deltas[i*3+j] = p21 + emis[3*i+j];
	  phis[i*3+j] = 1;
	}
	else {
	  deltas[i*3+j] = p22 + emis[3*i+j];
	  phis[i*3+j] = 2;
	}
      }
    }

    if ( deltas[(m-1)*3] > deltas[(m-1)*3+1] ) {
      if ( deltas[(m-1)*3] > deltas[(m-1)*3+2] ) paths[m-1] = 0;
      else paths[m-1] = 2;
    }
    else if ( deltas[(m-1)*3+1] > deltas[(m-1)*3+2] ) paths[m-1] = 1;
    else paths[m-1] = 2;
    for(i=m-1; i > 0; --i) {
      paths[i-1] = phis[i*3+paths[i]];
    }
  }

  // write pairwise HMM viterbi outputs
  int writeViterbi(wFile& wf) {
    // magic + n + m
    std::vector<int> ints;
    std::vector<uint8_t> nsh;
    for(int i=1; i < m; ++i) {
      if ( paths[i] != paths[i-1] ) {
	ints.push_back(i);
	nsh.push_back(paths[i-1]);
      }
    }
    if ( ints.empty() || ( ints.back() != m-1 ) ) {
      ints.push_back(m-1);
      nsh.push_back(paths[m-1]);
    }

    int nints = (int)ints.size();
    wf.write(&nints,sizeof(int));
    wf.write(ints.data(), sizeof(int)*nints);
    wf.write(nsh.data(), sizeof(uint8_t)*nints);
    return nints;
  }

  void printViterbi(wFile& wf) {
    int i, b = 0;
    for(i=1; i < m; ++i) {
      if ( paths[i] != paths[i-1] ) {
	wf.printf("\t%d,%d:%d",b,i-1,paths[i-1]);
	b = i;
      }
    }
    wf.printf("\t%d,%d:%d\n",b,m-1,paths[m-1]);
  }
};

#endif
