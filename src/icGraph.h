#ifndef __IC_GRAPH__H
#define __IC_GRAPH__H

#include <vector>
//#include <random>
#include "Error.h"
#include "PhredHelper.h"
#include "pFile.h"
#include "wFile.h"

typedef std::vector<int> vecint;
typedef std::vector< std::vector<int> > vecvecint;

class icInd;

class icGraph {
 public:
  int m;
  int n;
  std::vector<icInd> inds;
  std::vector<double> AFs;
  std::vector<std::string> markers;
  std::vector<int> pos1s;
  std::vector<std::string> refs;
  std::vector<std::string> alts;
  std::string chrom;

  static bool readIDs(pFile& tf, int n, std::vector<std::string>& ids);
  static bool readID(pFile& tf, std::string& id);
  void loadPairHMM(const char* filename);
  void print(wFile& wf);
};

class icInd {
 protected:
  icInd() {}
  icGraph* pGraph;

 public:
  std::string id;
  int nid;
  std::map<int,vecint> blocks;
  std::vector<uint8_t> PLs;  // size m * 3
  std::vector<uint8_t> haps; // size m

 icInd(int m, icGraph* graph) : pGraph(graph), PLs(3*m), haps(m) {
    blocks[m-1].clear(); // create one element
  }

  inline uint8_t makeHap(uint8_t h1, uint8_t h2) {
    return (uint8_t)(((h1 << 4) & 0x00f0) | (h2 & 0x000f));
  }

  void printHaplotypes() {
    printf("%s\t%d\t",id.c_str(), nid);
    for(int i=0; i < 50; ++i) {
      printf("%d",haps[i] > 15 ? haps[i]-14 : haps[i]);
    }
    printf("\n");
  }

  //void randomPhase(std::vector<double>& AFs, std::mt19937& rng) {
  void randomPhase() {
    int m = (int)haps.size();
    double p0, p1, p2, f, u;
    //std::uniform_real_distribution<double> uni(0,1);
    for(int i=0; i < m; ++i) {
      f = pGraph->AFs[i];
      p0 = (1-f)*(1-f)*phredConv.phred2Err[PLs[i*3+0]];
      p1 = 2*f*(1-f)*phredConv.phred2Err[PLs[i*3+1]];
      p2 = f*f*phredConv.phred2Err[PLs[i*3+2]];
      u = (rand()+0.5)/(RAND_MAX+1.)*(p0+p1+p2);
      if ( u < p1 ) {
	if ( u < p1/2 ) haps[i] = makeHap(0,1);
	else haps[i] = makeHap(1,0);
      }
      else if ( u < p0+p1 ) haps[i] = makeHap(0,0);
      else haps[i] = makeHap(1,1);
    }

    //for(int i=0; i < 10; ++i) fprintf(stderr,"%d",haps[i]>15 ? haps[i]-14 : haps[i]); error("foo");
  }

  void addBlock(int beg, int end, int iind) {
    //notice("Starting to call addBlock(%d, %d, %d) for %d",beg,end,iind,nid);
    //printBlock();
    // lo <= beg <= end < up is guaranteed
    std::map<int,vecint>::iterator lo = blocks.upper_bound(beg); 
    if ( lo != blocks.end() ) --lo;
    std::map<int,vecint>::iterator up = blocks.upper_bound(end);
    std::map<int,vecint>::iterator prev = blocks.end();
    //error("%d %d",lo == blocks.end() ? -1 : lo->first, up == blocks.end() ? -1 : up->first);
    if ( lo == blocks.end() || lo->first > beg ) { // if lo is not found, have to add [0,beg] at first
      lo = blocks.begin();
      if ( beg >= 0 )
	blocks[beg] = lo->second; // add [0,beg] without iid
      //notice("foo");
    }
    else if ( lo->first <= beg ) {  
      prev = lo;
      ++lo;
      blocks[beg] = lo->second; // add (lo,beg] without iid
      //notice("bar");
    }
    else {
      //notice("joe %d %d",lo->first,beg);
    }
    while( lo != up ) {
      if ( lo->first <= beg ) error("lo->first = %d, beg = %d, end = %d",lo->first,beg,end);
      if ( ( !lo->second.empty() ) && ( lo->second.back() == iind ) ) {
	//printBlock();
	error("inconsistent state %d, trying to add %d - %d for %d to %d - %d",nid, beg,end,iind,prev == blocks.end() ? -1 : prev->first, lo->first);
      }
      lo->second.push_back(iind);
      prev = lo;
      ++lo;
    }
    if ( ( prev->first < end ) && ( lo != blocks.end() ) && ( lo == up ) ) {
      blocks[end] = up->second;
      blocks[end].push_back(iind);
    }
    //notice("Successfully finished addBlock(%d, %d, %d) for %d",beg,end,iind,nid);
    //printBlock();
    //fprintf(stderr,"\n");
  }

  void printBlock(wFile& wf) {
    wf.printf("%s\t%d",id.c_str(),nid);
    int i, prev = 0;
    std::vector<bool> bmap(pGraph->n,false); 
    for(std::map<int,vecint>::iterator it = blocks.begin(); it != blocks.end(); ++it) {
      wf.printf("\t%d-%d:%d:",prev,it->first,(int)it->second.size());
      std::vector<bool> bmap2(pGraph->n,false);
      for(i=0; i < (int)it->second.size(); ++i) {
	if ( i > 0 ) wf.printf(",");
	wf.printf("%d",it->second[i]);
	bmap2[it->second[i]] = true;
      }
      if ( it->second.size() == 0 ) wf.printf(".");

      wf.printf(":");
      for(i=0; i < pGraph->n; ++i) {
	if ( ( bmap[i] == false ) && ( bmap2[i] == true ) ) {
	  wf.printf("+%d",i);
	}
	else if ( ( bmap2[i] == false ) && ( bmap[i] == true ) ) {	
	  wf.printf("-%d",i);
	}
      }
      bmap = bmap2;
      prev = it->first+1;
    }
    wf.printf("\n");
  }
};

bool icGraph::readIDs(pFile& tf, int n, std::vector<std::string>& ids) {
  ids.clear();
  char s[65536];
  int sz;
  for(int i=0; i < n; ++i) {
    tf.read(&sz,sizeof(int));
    tf.read(s,sz);
    ids.push_back(std::string(s,sz));
  }
  return true;
}

bool icGraph::readID(pFile& tf, std::string& id) {
  char s[65536];
  int sz;
  tf.read(&sz,sizeof(int));
  tf.read(s,sz);
  id.assign(s,sz);
  return true;
}

void icGraph::loadPairHMM(const char* filename) {
  int i, j, k, l, nints;
  pFile tf(filename);
  const char* MAGIC = "CHAPS_HMM";
  int N_MAGIC = strlen(MAGIC);
  char buf[N_MAGIC];
  
  tf.read(buf,N_MAGIC);
  if ( strncmp(buf,MAGIC,N_MAGIC) != 0 ) error("MAGIC number for PairHMM does not match");
  tf.read(&n, sizeof(int));
  tf.read(&m, sizeof(int));
  
  int* ints = new int[m];
  uint8_t* nshs = new uint8_t[m];
  
  inds.resize(n, icInd(m, this));
  
  for(i=0; i < n; ++i) {
    readID(tf, inds[i].id);
    inds[i].nid = i;
  }
  readIDs(tf, m, markers);
  for(i=0; i < m; ++i) { // marker is in the form of [chr]:[pos]_[ref]/[alt]
    j = markers[i].find(':');
    k = markers[i].find('_');
    if ( i == 0 ) {
      chrom = markers[i].substr(0,j);
    }
    pos1s.push_back(atoi(markers[i].substr(j+1,k-j-1).c_str()));
    j = markers[i].find('/');
    refs.push_back(markers[i].substr(k+1,j-k-1));
    k = markers[i].find_first_of('_',j);
    if ( k == (int)std::string::npos ) alts.push_back(markers[i].substr(j+1));
    else alts.push_back(markers[i].substr(j+1,k-j-1));
  }
  
  AFs.resize(m);
  tf.read(AFs.data(), m * sizeof(double));
  for(i=0; i < n; ++i) {
    tf.read(inds[i].PLs.data(), m * 3);
  }

  std::string outf = filename;
  outf += ".hkin";
  wFile wf(outf.c_str());
  for(i=0; i < n; ++i) {
    if ( i > 0 ) wf.printf("\t");
    wf.printf("%s",inds[i].id.c_str());
  }
  // construct blocks
  std::vector<double> K;
  for(i=0; i < n; ++i) {
    for(j=0; j < i; ++j) {
      tf.read(&nints,sizeof(int));
      tf.read(ints, sizeof(int)*nints);
      tf.read(nshs, nints);
      
      // compact the intervals to merge all positive values
      int n0 = 0;
      int n1 = 0;

      for(k=0; k < nints; ++k) {
	if ( nshs[k] > 0 ) {
	  n1 += (k > 0 ? ints[k]-ints[k-1] : ints[k]);
	}
	else {
	  n0 += (k > 0 ? ints[k]-ints[k-1] : ints[k]);
	}
      }
      K.push_back((double)n1/(double)(n0+n1));

      for(k=0, l=0; k < nints;) {
	if ( nshs[k] > 0 ) {
	  while( ( nshs[k] > 0 ) && ( k < nints ) ) ++k;
	  nshs[l] = 1;
	  ints[l] = ints[k-1];
	}
	else {
	  ++k;
	  nshs[l] = 0;
	  ints[l] = ints[k-1];
	}
	++l;
      }
      nints = l;
      
      if ( ints[0] == 0 ) printf("%d %d %d %d %d %d\n",i,j,ints[0],ints[1],nshs[0],nshs[1]);
      
      for(k=0; k < nints; ++k) {
	if ( nshs[k] > 0 ) {
	  inds[i].addBlock(k == 0 ? -1 : ints[k-1], ints[k], j);
	  inds[j].addBlock(k == 0 ? -1 : ints[k-1], ints[k], i);
	}
      }
    }
  }
  
  for(i=0; i < n; ++i) {
    wf.printf("%s",inds[i].id.c_str());
    for(j=0; j < n; ++j) {
      if ( i == j ) {
	wf.printf("\t%.5lf",1);
      }
      else if ( i > j ) {
	wf.printf("\t%.5lf",K[i*(i+1)/2+j]);
      }
      else {
	wf.printf("\t%.5lf",K[j*(j+1)/2+i]);
      }
    }
    wf.printf("\n");
  }
  wf.close();
  //error("foo");
  delete [] ints;
  delete [] nshs;
}

void icGraph::print(wFile& wf) {
  //fprintf(stderr,"m = %d, n = %d\n",m,n);
  for(int i=0; i < n; ++i) {
    inds[i].printBlock(wf);
  }
}

#endif
