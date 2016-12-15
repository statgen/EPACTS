#ifndef __SPARSE_HMM__H
#define __SPARSE_HMM__H

#include <vector>
//#include <random>
#include <ctime>
#include <algorithm>
#include <utility>
#include "Error.h"
#include "PhredHelper.h"
#include "wFile.h"
#include "icGraph.h"

class sparseHMM {
 public:
  icGraph graph;
  int minStates;
  int maxStates;
  int rndStates;
  int maxBlocks;
  int states;
  uint32_t seed;
  //std::mt19937 rng;

  std::vector<double> trans;
  std::vector<double> alphas;
  std::vector<double> alphaSums;
  std::vector<double> betas;
  std::vector<double> betaSums;
  std::vector<double> gammas;
  std::vector<double> emis;
  std::vector<bool> bsels;
  std::vector<int> isels;
  std::vector<int> shuf;
  double mu;
  double theta;
  //double theta0; // within blocks : 1e-9
  //double theta1; // between blocks : 1e-7

 sparseHMM(const char* file, int _min, int _max, int _rnd, int _seed) : minStates(_min/2), maxStates(_max/2), rndStates(_rnd/2), seed((uint32_t)_seed) {
    // read input file and construct graph
    graph.loadPairHMM(file);
    //if ( seed == 0 ) seed = (uint32_t)std::time(0);
    //rng.seed(seed);
    if ( seed == 0 ) srand(std::time(0));
    else srand(seed);
    mu = 1e-3;
    theta = 1e-8;
    //theta0 = 1e-9;
    //theta1 = 1e-7;

    bsels.resize(graph.n);
    isels.resize(maxStates);
    for(int i=0; i < graph.n; ++i) shuf.push_back(i);
    std::random_shuffle(shuf.begin(), shuf.end());

    maxBlocks = 0; 
    for(int i=0; i < graph.n; ++i) {
      std::map<int,vecint>& b = graph.inds[i].blocks;
      if ( maxBlocks < (int)b.size() ) maxBlocks = (int)b.size();
      int prev = -1;
      for(std::map<int,vecint>::iterator it = b.begin(); it != b.end(); ++it) {
	if ( maxBlocks < it->first - prev + 1 ) {
	  maxBlocks = it->first - prev + 1;
	}
	prev = it->first;
      }
    }

    alphas.resize(maxStates * ( 2 * maxStates + 1 ) * maxBlocks );
    betas.resize(maxStates * ( 2 * maxStates + 1 ) * maxBlocks );
    emis.resize(maxStates * ( 2 * maxStates + 1 ) * maxBlocks );
    alphaSums.resize(maxStates * 2 * maxBlocks );
    betaSums.resize(maxStates * 2 * maxBlocks );
    trans.resize( maxBlocks );
  }

  std::pair<int,int> blockBound(int iind, int isnp) {
    std::map<int,vecint> &b = graph.inds[iind].blocks;
    std::map<int,vecint>::iterator lo = b.lower_bound(isnp);
    if ( lo == b.begin() ) {
      return std::pair<int,int>(-1, lo->first);
    }
    else {
      int r = lo->first;
      --lo;
      return std::pair<int,int>(lo->first, r);
    }
  }


  std::pair<int,int> vecSample(std::vector<double>& prods, int nstates, double sum) {
    int i,j,k;
    double u = (rand()+.5)/(RAND_MAX+1.) * sum;
    for(i=0, k=0; i < nstates + nstates; ++i) {
      for(j=0; j < i; ++j, ++k) {
	u -= (2*prods[k]);
	if ( u < 0 ) return std::pair<int,int>(i,j);
      }
      u -= prods[k];
      if ( u < 0 ) return std::pair<int,int>(i,j);
      ++k;
    }
    error("marginalSample() reached at the end");
    return std::pair<int,int>(nstates+nstates-1,nstates+nstates-2);
  }

  std::pair<int,int> marginalSample(std::vector<double>& alphas, std::vector<double>& betas, int nstates, int ib) {
    int i,j,k;
    double sum = 0;
    int maxPair = maxStates * ( 2 * maxStates + 1);
    std::vector<double> prods;
    for(i=0, k=0; i < nstates + nstates; ++i) {
      for(j=0; j <= i; ++j, ++k) {
	prods.push_back(alphas[k + ib*maxPair] * betas[k + ib*maxPair]);
	sum += (prods.back() * ( j == i ? 1 : 2 ));
      }
    }
    return vecSample(prods, nstates, sum);
  }

  // select individuals
  void selectEdges(int iind, int isnp) {
    // if the number of individuals is too small, do not worry
    std::fill(isels.begin(), isels.end(), -1);
    if ( minStates >= graph.n -1 ) {
      std::fill(bsels.begin(), bsels.end(), true);
      bsels[iind] = false;
      for(int i=0; i < graph.n; ++i) {
	if ( i != iind ) {
	  isels[i] = (i > iind) ? i-1 : i;
	}
      }
      states = graph.n - 1;
    }
    else {
      int i, nrnd, nrmv;

      std::fill(bsels.begin(), bsels.end(), false);

      // obtain block boundary
      //std::pair<int,int> bnds = blockBound(iind,isnp);
      //vecint& v = graph.inds[iid].blocks[bnds.second];
      vecint& v = graph.inds[iind].blocks[isnp];

      // mark all the connected individuals
      for(i=0; i < (int)v.size(); ++i) {
	if ( v[i] == iind ) error("v[i] = %d = iind",iind);
	bsels[i] = true;
      }
      bsels[iind] = true;

      //std::random_shuffle(shuf.begin(), shuf.end());

      if ( minStates > (int)v.size() + rndStates ) nrnd = minStates - v.size();
      else if ( maxStates > (int)v.size() + rndStates ) nrnd = rndStates;
      else { // exceed maxStates limit;
	nrnd = rndStates;
	nrmv = (int)v.size() + rndStates - maxStates;
      }

      if ( nrmv > 0 ) {
	for(i=0; i < graph.n; ++i) {
	  if ( bsels[v[i]] && ( nrmv > 0 ) ) {
	    bsels[v[i]] = false;
	    --nrmv;
	  }
	}
	//std::random_shuffle(shuf.begin(), shuf.end());
      }
      if ( nrnd > 0 ) {
	for(i=0; i < graph.n; ++i) {
	  if ( ( bsels[v[i]] == false ) && ( nrnd > 0 ) ) {
	    bsels[v[i]] = true;
	    --nrnd;
	  }
	}
      }

      bsels[iind] = false;
      states = 0;
      for(i=0; i < graph.n; ++i) {
	if ( bsels[i] ) isels[states++] = i;
      }
    }
  }

  void writeHaplotypes(const char* outf) {
    int i, j;
    //notice("writeHaplotype(%s) called",outf);
    wFile wf(outf);
    //notice("foo");
    wf.printf("##fileformat=VCFv4.0\n");
    wf.printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for(i=0; i < graph.n; ++i) wf.printf("\t%s",graph.inds[i].id.c_str());
    wf.printf("\n");
    for(j=0; j < graph.m; ++j) {
      wf.printf("%s\t%d\t.\t%s\t%s\t100\tPASS\tAF=%.6lf\tGT",graph.chrom.c_str(), graph.pos1s[j], graph.refs[j].c_str(), graph.alts[j].c_str(), graph.AFs[j]);
      for(i=0; i < graph.n; ++i) {
	switch(graph.inds[i].haps[j]) {
	case 0:
	  wf.printf("\t0|0");
	  break;
	case 1:
	  wf.printf("\t0|1");
	  break;
	case 16:
	  wf.printf("\t1|0");
	  break;
	case 17:
	  wf.printf("\t1|1");
	  break;
	default:
	  error("Cannot recognize haplotype %d",graph.inds[i].haps[j]);
	}
      }
      wf.printf("\n");
      notice("bar %d %d",graph.m,j);
    }
    wf.close();
  }

  void initializeHaplotypes(bool random = false) {
    int i;
    for(i=0; i < graph.n; ++i) {
      graph.inds[i].randomPhase();
    }

    //for(i=0; i < graph.n; ++i) graph.inds[i].printHaplotypes();
    //printf("\n");
    if ( !random ) {
      //improveHaplotypes();
      //for(i=0; i < graph.n; ++i) graph.inds[i].printHaplotypes();
      //printf("\n");
      //improveHaplotypes();
      //for(i=0; i < graph.n; ++i) graph.inds[i].printHaplotypes();
      //printf("\n");
    }
  }

  inline void countHapMismatch(uint8_t h1, uint8_t h2, int* cnts) {
    bool a = (h1 & 0x00f0) ? true : false;
    bool b = (h1 & 0x000f) ? true : false;
    bool c = (h2 & 0x00f0) ? true : false;
    bool d = (h2 & 0x000f) ? true : false;
    if ( a != c ) ++cnts[0];
    if ( a != d ) ++cnts[1];
    if ( b != c ) ++cnts[2];
    if ( b != d ) ++cnts[3];
  }

  void improveHaplotypes() { // needs to change to E-M phasing
    int i, k, l, r, lo, hi;
    int c[4];
    
    // re-assign HETs by examining compatible haplotypes
    for(i=0; i < graph.n; ++i) {
      std::map<int,vecint>& b = graph.inds[i].blocks;
      lo = -1;
      for(std::map<int,vecint>::iterator it = b.begin(); it != b.end(); ++it) {
	hi = it->first;
	vecint& v = it->second;

	// spot the het positions
	std::vector<int> ihets;
	std::vector<uint8_t>& h1 = graph.inds[i].haps;
	for(k=lo+1; k <= hi; ++k) {
	  if ( ( h1[k] == 1 ) || ( h1[k] == 16 ) )
	    ihets.push_back(k);
	}
	// between the two possible haplotypes, find a closer pair
	if ( !ihets.empty() ) {
	  std::vector<int> votes(4*ihets.size(),0);
	  for(k=0; k < (int)v.size(); ++k) {
	    c[0] = c[1] = c[2] = c[3] = 0;
	    for(l = lo+1; l <= hi; ++l) {
	      countHapMismatch(graph.inds[i].haps[l],graph.inds[v[k]].haps[l],c);
	    }
	    if ( c[0] < c[1] ) {
	      if ( c[2] < c[3] ) {
		if ( c[0] < c[2] ) { // a-c
		  for(r=0; r < (int)ihets.size(); ++r) {
		    if ( graph.inds[v[k]].haps[ihets[r]] & 0xf0 ) ++votes[1+4*r];
		    else ++votes[0+4*r];
		  }
		}
		else { // b-c
		  for(r=0; r < (int)ihets.size(); ++r) {
		    if ( graph.inds[v[k]].haps[ihets[r]] & 0xf0 ) ++votes[3+4*r];
		    else ++votes[2+4*r];
		  }
		}
	      }
	      else {
		if ( c[0] < c[3] ) { // a-c
		  for(r=0; r < (int)ihets.size(); ++r) {
		    if ( graph.inds[v[k]].haps[ihets[r]] & 0xf0 ) ++votes[1+4*r];
		    else ++votes[0+4*r];
		  }
		}
		else { // b-d
		  for(r=0; r < (int)ihets.size(); ++r) {
		    if ( graph.inds[v[k]].haps[ihets[r]] & 0x0f ) ++votes[3+4*r];
		    else ++votes[2+4*r];
		  }
		}
	      }
	    }
	    else {
	      if ( c[2] < c[3] ) {
		if ( c[1] < c[2] ) { // a-d
		  for(r=0; r < (int)ihets.size(); ++r) {
		    if ( graph.inds[v[k]].haps[ihets[r]] & 0x0f ) ++votes[1+4*r];
		    else ++votes[0+4*r];
		  }
		}
		else { // b-c
		  for(r=0; r < (int)ihets.size(); ++r) {
		    if ( graph.inds[v[k]].haps[ihets[r]] & 0xf0 ) ++votes[3+4*r];
		    else ++votes[2+4*r];
		  }
		}
	      }
	      else {
		if ( c[1] < c[3] ) { // a-d
		  for(r=0; r < (int)ihets.size(); ++r) {
		    if ( graph.inds[v[k]].haps[ihets[r]] & 0x0f ) ++votes[1+4*r];
		    else ++votes[0+4*r];
		  }
		}
		else { // b-d
		  for(r=0; r < (int)ihets.size(); ++r) {
		    if ( graph.inds[v[k]].haps[ihets[r]] & 0x0f ) ++votes[3+4*r];
		    else ++votes[2+4*r];
		  }
		}
	      }
	    }
	  }
	  // re-assign hets
	  for(r=0; r < (int)ihets.size(); ++r) {
	    if ( votes[0+4*r] < votes[1+4*r] ) {
	      if ( votes[2+4*r] > votes[3+4*r] ) { // 0-1
		h1[ihets[r]] = (uint8_t)1;
	      }
	    }
	    else if ( votes[2+4*r] < votes[3+4*r] ) {
		h1[ihets[r]] = (uint8_t)16;
	    }
	  }
	}
      }
    }
  }

  inline int getGenotype(int s1, int s2, int isnp) {
    return ( ((graph.inds[s1/2].haps[isnp] >> ((s1 % 2) * 4)) & 0x0f) + 
	     ((graph.inds[s2/2].haps[isnp] >> ((s2 % 2) * 4)) & 0x0f) );
  }

  inline double getLogEmis(uint8_t* PLs, int h1, int h2, double mu) {
    return log(phredConv.phred2Err[PLs[h1+h2]] + mu * phredConv.phred2Err[PLs[1-h1+h2]] + mu * phredConv.phred2Err[PLs[1-h2+h1]] + mu * mu * phredConv.phred2Err[PLs[2-h1-h2]]);
  }

  inline double getEmis(uint8_t* PLs, int h1, int h2, double mu) {
    return phredConv.phred2Err[PLs[h1+h2]] + mu * phredConv.phred2Err[PLs[1-h1+h2]] + mu * phredConv.phred2Err[PLs[1-h2+h1]] + mu * mu * phredConv.phred2Err[PLs[2-h1-h2]];
  }

  void betweenForwardBackward(int iind) {
    int i,j,k,l;

    std::map<int,vecint>& b = graph.inds[iind].blocks;
    std::vector<uint8_t>& PLs = graph.inds[iind].PLs;
    std::vector<int> nstates;
    vecvecint edges; // list of edges selected for each block
    vecint isnps;
    int nb = b.size();
    //vecint nedges;
    int maxPair = maxStates * ( 2 * maxStates + 1);
    int maxHap = maxStates * 2;
    std::vector<int> indices(nb * graph.n, -1);
    
    isnps.push_back(-1);
    int ib = 0;
    int prevsz = 0, sz = 0;
    for(std::map<int,vecint>::iterator it = b.begin(); it != b.end(); ++it, ++ib) {
      notice("Processing block %d",ib);
      selectEdges(iind, it->first);
      edges.push_back(vecint(states));
      for(i=0; i < (int)isels.size(); ++i) {
	indices[ isels[i] + ib * nb ] = i;
      }
      std::copy( isels.begin(), isels.begin() + states, edges.back().begin());
      // calculating emission probability within each block
      notice("Calculating emission probability for block %d",ib);
      double maxLogE = -1e99;
      nstates.push_back(states);
      for(i=0, l=0; i < states * 2; ++i) {
	for(j=0; j <= i; ++j, ++l) {
	  double& e = emis[l + ib * maxPair];
	  e = 0;
	  for(k=isnps.back()+1; k <= it->first; ++k) {
	    e += getLogEmis(PLs.data()+3*k, (graph.inds[i/2].haps[k] >> ((i % 2) * 4)) & 0x0f, (graph.inds[j/2].haps[k] >> ((j % 2) * 4)) & 0x0f, mu);
	  }
	  if ( maxLogE < e ) maxLogE = e;
	}
      }

      for(i=0, k=0; i < states * 2; ++i) {
	for(j=0; j <= i; ++j, ++k) {
	  double& e = emis[k + ib * maxPair];
	  e = exp(e-maxLogE);  // make normalized emission probability
	}
      }

      sz = (ib == 0) ? (graph.pos1s[it->first]-graph.pos1s.front()+1) : ( graph.pos1s[it->first] - graph.pos1s[isnps.back()] );
      isnps.push_back(it->first);
      //nedges.push_back(states);
      //npairs.push_back(states*(2*states+1));
      if ( ib > 0 ) trans[ib-1] = theta * (sz + prevsz) / 2.; // linear approx.
    }

    // calculate forward probability
    // start with uniform priors
    double sum = 0;
    for(i=0, k=0; i < nstates[0]+nstates[0]; ++i) {
      for(j=0; j < i; ++j, ++k) {
	sum += (2*(alphas[k] = emis[k]));
      }
      sum += (alphas[k] = emis[k]);
      ++k;
    }
    // normalize the probabilities
    std::fill(alphaSums.begin(), alphaSums.end(), 0);
    for(i=0, k=0; i < nstates[0]+nstates[0]; ++i) {
      for(j=0; j < i; ++j, ++k) {
	alphaSums[i] += (alphas[k] /= sum);
	alphaSums[j] += alphas[k];
      }
      alphaSums[i] += (alphas[k] /= sum);
      ++k;
    }
    // make transition along the blocks
    // transition is made from different number of states
    int iip, ijp, kp;
    for(ib = 1; ib < nb; ++ib) {
      sum = 0;
      for(i=0, k=0; i < nstates[ib]+nstates[ib]; ++i) {
	iip = indices[ (ib-1) * graph.n + edges[ib][i/2] ] * 2 + ( i % 2 );
	//iin = indices[ ib * graph.n + edges[ib][i/2] ];
	for(j=0; j <= i; ++j, ++k) {
	  ijp = indices[ (ib-1) * graph.n + edges[ib][j/2] ] * 2 + ( j % 2 );
	  //ijn = indices[ ib * graph.n + edges[ib][j/2] ];
	  if ( iip > 0 ) {
	    if ( ijp > 0 ) {
	      kp = iip * ( iip + 1 ) / 2 + ijp;
	      alphas[k + ib * maxPair] = ( alphas[ kp + (ib-1) * maxPair] * (1-trans[ib-1]) * (1-trans[ib-1]) + ( alphaSums[ iip + (ib-1) * maxHap] + alphaSums[ ijp + (ib-1) * maxHap ] ) * trans[ib-1] * (1 - trans[ib-1]) / ( nstates[ib-1] * 2 ) + trans[ib-1] * trans[ib-1] / ( nstates[ib-1] * nstates[ib-1] * 4 ) ) * emis[ k + ib * maxPair ];
	    }
	    else {
	      alphas[k + ib * maxPair] = ( alphaSums[ iip + (ib-1) * maxHap] * trans[ib-1] * (1 - trans[ib-1]) / ( nstates[ib-1] * 2 ) + trans[ib-1] * trans[ib-1] / ( nstates[ib-1] * nstates[ib-1] * 4 ) ) * emis[ k + ib * maxPair ];
	    }
	  }
	  else if ( ijp > 0 ) {
	    alphas[k + ib * maxPair] = ( alphaSums[ ijp + (ib-1) * maxHap] * trans[ib-1] * (1 - trans[ib-1]) / ( nstates[ib-1] * 2 ) + trans[ib-1] * trans[ib-1] / ( nstates[ib-1] * nstates[ib-1] * 4 ) ) * emis[ k + ib * maxPair ];
	  }
	  else {
	    alphas[k + ib * maxPair] = ( trans[ib-1] * trans[ib-1] / ( nstates[ib-1] * nstates[ib-1] * 4 ) ) * emis[ k + ib * maxPair ];
	  }
	  sum += (alphas[k + ib * maxPair] * (j == i ? 1 : 2));
	  alphaSums[i + ib*maxHap] += alphas[k + ib*maxPair];
	  if ( i != j ) alphaSums[j + ib*maxHap] += alphas[k + ib*maxPair];
	}
      }
      for(i=0, k=0; i < nstates[ib]+nstates[ib]; ++i) {
	for(j=0; j <= i; ++j, ++k) {
	  alphas[k + ib*maxPair] /= sum;
	}
	alphaSums[i + ib * maxHap] /= sum;
      }
    }

    // calculate the initial backward probability
    sum = 0;
    std::fill(betaSums.begin(), betaSums.end(), 0);
    for(i=0, k=0; i < nstates[nb-1]+nstates[nb-1]; ++i) {
      for(j=0; j < i; ++j, ++k) {
	double& e = emis[ k + (nb-1) * maxPair ];
	sum += (2 * e);
	betaSums[i + (nb-1) * maxHap ] += e;
	betaSums[j + (nb-1) * maxHap ] += e;
      }
      double& e = emis[ k + (nb-1) * maxPair ];
      sum += e;
      betaSums[i] += e;
      ++k;
    }
    for(i=0, k=0; i < nstates[nb-1]+nstates[nb-1]; ++i) {
      for(j=0; j < i; ++j, ++k) {
	betas[k + (nb-1) * maxPair] = 1/sum;
      }
      betaSums[i + (nb-1) * maxHap] /= sum;
      ++k;
    }

    // iterate backward algorithm
    int iin, ijn, kn;
    double betaE;
    for(ib = nb-1; ib > 0; ++ib) {
      sum = 0;
      for(i=0, k=0; i < nstates[ib-1]+nstates[ib-1]; ++i) {
	iin = indices[ ib * graph.n + edges[ib-1][i/2] ] * 2 + ( i % 2 );
	for(j=0; j <= i; ++j, ++k) {
	  ijn = indices[ ib * graph.n + edges[ib-1][j/2] ] * 2 + ( j % 2 );
	  if ( iin > 0 ) {
	    if ( ijn > 0 ) {
	      kn = iin * ( iin + 1 ) / 2 + ijn;
	      betas[k + (ib-1) * maxPair] = betas[ kn + ib * maxPair] * (1-trans[ib-1]) * (1-trans[ib-1]) * emis[ kn + ib * maxPair] + ( betaSums[ iin + ib * maxHap] + betaSums[ ijn + ib * maxHap ] ) * trans[ib-1] * (1 - trans[ib-1]) / ( nstates[ib] * 2 ) + trans[ib-1] * trans[ib-1] / ( nstates[ib] * nstates[ib] * 4 );
	    }
	    else {
	      betas[k + (ib-1) * maxPair] = ( betaSums[ iin + ib * maxHap] ) * trans[ib-1] * (1 - trans[ib-1]) / ( nstates[ib] * 2 ) + trans[ib-1] * trans[ib-1] / ( nstates[ib] * nstates[ib] * 4 );
	    }
	  }
	  else if ( ijn > 0 ) {
	    betas[k + (ib-1) * maxPair] = ( betaSums[ ijn + ib * maxHap ] ) * trans[ib-1] * (1 - trans[ib-1]) / ( nstates[ib] * 2 ) + trans[ib-1] * trans[ib-1] / ( nstates[ib] * nstates[ib] * 4 );
	  }
	  else {
	    betas[k + (ib-1) * maxPair] = trans[ib-1] * trans[ib-1] / ( nstates[ib] * nstates[ib] * 4 );
	  }

	  betaE = betas[k + (ib-1)*maxPair] * emis[ k + (ib-1) * maxPair];
	  sum += (betaE * (j==i ? 1 : 2));
	  betaSums[i + (ib-1)*maxHap] += betaE;
	  if ( j != i ) betaSums[j + (ib-1)*maxHap] += betaE;
	}
      }
      for(i=0, k=0; i < nstates[ib-1]+nstates[ib-1]; ++i) {
	for(j=0; j <= i; ++j, ++k) {
	  betas[k + (ib-1)*maxPair] /= sum;
	}
	betaSums[i + (ib-1)*maxHap] /= sum;
      }
    }

    // sample haplotypes
    std::vector<int> hap1s;
    std::vector<int> hap2s;
    std::pair<int,int> rp = marginalSample(alphas, betas, nstates[0], 0);
    iip = indices[ edges[0][rp.first/2] ] * 2 + rp.first % 2;
    ijp = indices[ edges[0][rp.second/2] ] * 2 + rp.second % 2;
    double t;
    std::vector<double> prods;
    hap1s.push_back(iip); hap2s.push_back(ijp);
    for(ib=1; ib < nb; ++ib) {
      sum = 0;
      prods.clear();
      for(i=0, k=0; i < nstates[ib]+nstates[ib]; ++i) {
	iin = indices[ ib * graph.n + edges[ib][rp.first/2] ] * 2 + rp.first % 2;
	for(j=0; j <= i; ++j, ++k) {
	  ijn = indices[ ib * graph.n + edges[ib][rp.second/2] ] * 2 + rp.second  % 2;
	  if ( ( iip == iin ) && ( ijp == ijn ) ) {
	    t = (1-trans[ib-1]*(nstates[ib]*2-1)/(nstates[ib]*2));
	    t = t * t;
	  }
	  else if ( ( iip == iin ) || ( ijp == ijn ) || ( iip == ijn ) || ( ijp == iin ) ) {
	    t = (1-trans[ib-1]*(nstates[ib]*2-1)/(nstates[ib]*2)) * trans[ib-1] / (nstates[ib]*2);
	  }
	  else {
	    t = trans[ib-1] * trans[ib-1] / (nstates[ib]*nstates[ib]*4);
	  }
	  prods.push_back(t * emis[k + ib*maxPair] * betas[k + ib*maxPair]);
	  sum += (prods.back() * ( i == j ? 1 : 2 ));
	}
      }
      rp = vecSample(prods, nstates[ib],sum);
      iip = indices[ edges[0][rp.first/2] ] * 2 + rp.first % 2;
      ijp = indices[ edges[0][rp.second/2] ] * 2 + rp.second % 2;
      hap1s.push_back(iip); hap2s.push_back(ijp);
    }

    for(ib=0; i < nb; ++ib) {
      for(i=isnps[ib]+1; i <= isnps[ib+1]; ++i) {
	graph.inds[iind].haps[i] = ((graph.inds[hap1s[ib]/2].haps[i] >> ((hap1s[ib] % 2) * 4) & 0x0f) << 4) | (graph.inds[hap2s[ib]/2].haps[i] >> ((hap2s[ib] % 2)*4) & 0x0f);
      }
    }
  }

  // run forward-backward within a block
  void withinForwardBackward(int iind, int isnp) {
    int i,j,k,ib;
    int beg, end;

    std::map<int,vecint>& b = graph.inds[iind].blocks;
    std::map<int,vecint>::iterator lo = b.lower_bound(isnp);
    end = isnp;
    if ( lo == b.begin() ) beg = 0;
    else { --lo; beg = lo->first + 1; }
    selectEdges(iind, isnp);
    std::vector<uint8_t>& PLs = graph.inds[iind].PLs;
    int maxPair = maxStates * ( 2 * maxStates + 1);
    int maxHap = maxStates * 2;
    
    // calculate emission probability for each SNP
    for(ib = 0; ib <= end-beg; ++ib) {
      for(i=0, k=0; i < states * 2; ++i) {
	for(j=0; j <= i; ++j, ++k) {
	  emis[k + ib * maxPair] = getEmis(PLs.data()+3*k, (graph.inds[i/2].haps[isnp] >> ((i % 2) * 4)) & 0x0f, (graph.inds[j/2].haps[isnp] >> ((j % 2) * 4)) & 0x0f, mu);
	}
      }
      
    }
    int nb = end-beg+1;
    
    // calculate forward probability
    // start with uniform priors
    double sum = 0;
    for(i=0, k=0; i < states+states; ++i) {
      for(j=0; j < i; ++j, ++k) {
	sum += (2*(alphas[k] = emis[k]));
      }
      sum += (alphas[k] = emis[k]);
      ++k;
    }
    // normalize the probabilities
    std::fill(alphaSums.begin(), alphaSums.end(), 0);
    for(i=0, k=0; i < states+states; ++i) {
      for(j=0; j < i; ++j, ++k) {
	alphaSums[i] += (alphas[k] /= sum);
	alphaSums[j] += alphas[k];
      }
      alphaSums[i] += (alphas[k] /= sum);
      ++k;
    }

    // make transition along the blocks
    for(ib = 1; ib < nb; ++ib) {
      sum = 0;
      for(i=0, k=0; i < states + states; ++i) {
	for(j=0; j <= i; ++j, ++k) {
	  alphas[k + ib * maxPair] = ( alphas[ k + (ib-1) * maxPair] * (1-trans[ib-1]) * (1-trans[ib-1]) + ( alphaSums[ i + (ib-1) * maxHap] + alphaSums[ j + (ib-1) * maxHap ] ) * trans[ib-1] * (1 - trans[ib-1]) / ( states * 2 ) + trans[ib-1] * trans[ib-1] / ( states * states * 4 ) ) * emis[ k + ib * maxPair ];
	}
      }
      for(i=0, k=0; i < states + states; ++i) {
	for(j=0; j <= i; ++j, ++k) {
	  alphas[k + ib*maxPair] /= sum;
	}
	alphaSums[i + ib * maxHap] /= sum;
      }
    }

    // calculate the initial backward probability
    sum = 0;
    std::fill(betaSums.begin(), betaSums.end(), 0);
    for(i=0, k=0; i < states + states; ++i) {
      for(j=0; j < i; ++j, ++k) {
	double& e = emis[ k + (nb-1) * maxPair ];
	sum += (2 * e);
	betaSums[i + (nb-1) * maxHap ] += e;
	betaSums[j + (nb-1) * maxHap ] += e;
      }
      double& e = emis[ k + (nb-1) * maxPair ];
      sum += e;
      betaSums[i] += e;
      ++k;
    }
    for(i=0, k=0; i < states + states; ++i) {
      for(j=0; j < i; ++j, ++k) {
	betas[k + (nb-1) * maxPair] = 1/sum;
      }
      betaSums[i + (nb-1) * maxHap] /= sum;
      ++k;
    }

    // iterate backward algorithm
    double betaE;
    for(ib = nb-1; ib > 0; ++ib) {
      sum = 0;
      for(i=0, k=0; i < states + states; ++i) {
	for(j=0; j <= i; ++j, ++k) {
	  betas[k + (ib-1) * maxPair] = betas[ k + ib * maxPair] * (1-trans[ib-1]) * (1-trans[ib-1]) * emis[ k + ib * maxPair] + ( betaSums[ i + ib * maxHap] + betaSums[ j + ib * maxHap ] ) * trans[ib-1] * (1 - trans[ib-1]) / ( states * 2 ) + trans[ib-1] * trans[ib-1] / ( states * states * 4 );

	  betaE = betas[k + (ib-1)*maxPair] * emis[ k + (ib-1) * maxPair];
	  sum += (betaE * (j==i ? 1 : 2));
	  betaSums[i + (ib-1)*maxHap] += betaE;
	  if ( j != i ) betaSums[j + (ib-1)*maxHap] += betaE;
	}
      }
      for(i=0, k=0; i < states + states; ++i) {
	for(j=0; j <= i; ++j, ++k) {
	  betas[k + (ib-1)*maxPair] /= sum;
	}
	betaSums[i + (ib-1)*maxHap] /= sum;
      }
    }
  }
  /*
  void forward(int i);
  void backward(int i);
  void forwardBackward(int i);
  void sampleHaplotype(int i);
  */
};

#endif
