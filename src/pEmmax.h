#ifndef __EPACTS_EMMAX_H
#define __EPACTS_EMMAX_H

#include <iostream>
#include <cstring>
#include <cmath>
#include <complex>
#include <cstdio>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "pFile.h"
#include "wFile.h"
#include "Error.h"
#include "genomeScore.bak.h"

#define MAGIC_KIN "EMMA_KIN"
#define MAGIC_EIG "EMMA_EIG"

using namespace Eigen;

class pEmmaxArgs {
public:
  // VCF-related string arguments
  std::string vcf;
  std::string region;
  std::string weight;
  std::string rule;
  std::string groupf;
  std::string markerset;

  // Other input files
  std::string pedf;
  std::string phenof;
  std::string scoref;
  std::string covf;
  std::string indf;
  std::string betaf;

  // Other output or intermediate files
  std::vector<std::string> kinfs;
  std::string ineigf;
  std::string outeigf;
  std::string kinf;
  std::string outf;
  std::string remlf;
  std::string tryf;
  std::string assocf;
  std::string scoredir;

  int unit;
  int bins;
  int nsimul;
  std::string field;
  bool IBS, BN, YV, adjDiag;
  bool verbose;
  bool ignoreFilter;
  bool sepchr;
  bool raw;
  bool cov;
  bool normalize;
  bool noIntercept;
  bool binary;
  bool recessive;
  bool summarize;
  bool compact;

  double minMAF;
  double maxMAF;
  int minMAC;
  int maxMAC;
  double minCallRate;
  double minRSQ;
  double h2;
  unsigned int seed;
  int maxperm;
  int minperm;
  int pca;
  double relpair;
  bool dump;
  bool maxkin;
  int digits;
  double maxP;

  static int const DEFAULT_UNIT = 10000L;
  constexpr static double const DEFAULT_MIN = 1e-6;
  constexpr static double const DEFAULT_MAX = 1;
  static int const DEFAULT_MIN_MAC = 1;
  static int const DEFAULT_MAX_MAC = 1000000000;
  static int const DEFAULT_MAX_PERM = 10000000;
  static int const DEFAULT_MIN_PERM = 100;
  static int const DEFAULT_DIGITS = 5;

 pEmmaxArgs() :
  unit(DEFAULT_UNIT), bins(4), nsimul(1), field("GT"), IBS(false), BN(true), YV(false), adjDiag(false), verbose(false), ignoreFilter(false), sepchr(false), raw(false), cov(false), normalize(false), noIntercept(false), binary(false), recessive(false), summarize(false), compact(false), minMAF(DEFAULT_MIN), maxMAF(DEFAULT_MAX), minMAC(DEFAULT_MIN_MAC), maxMAC(DEFAULT_MAX_MAC), minCallRate(DEFAULT_MIN), minRSQ(0), h2(0), seed(0), maxperm(DEFAULT_MAX_PERM), minperm(DEFAULT_MIN_PERM), pca(0), relpair(1.), dump(false), maxkin(false), digits(DEFAULT_DIGITS), maxP(1.)
    {}
};

class REML {
 public:
  double delta;
  double h2;
  double vg;
  double ve;
  double LLK1;
  double LLK0;
  double dLLK0;
  double ddLLK0;
 REML(double d, double h, double g, double e, double l1, double l0, double dl0, double ddl0) : delta(d), h2(h), vg(g), ve(e), LLK1(l1), LLK0(l0), dLLK0(dl0), ddLLK0(ddl0) {} 
};

class pEmmaxHelper {
public:
  static const int N_MAGIC = 8;
  static const char *magicKin;
  static const char *magicEig;
  constexpr static const double ZEPS = 1e-10;
  constexpr static const double TOL = 1e-6;

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

  static bool readIDs(pFile& tf, int n, std::vector<std::string>& ids) {
    ids.clear();
    char* s = new char[65536];
    int sz;
    for(int i=0; i < n; ++i) {
      tf.read(&sz,sizeof(int));
      tf.read(s,sz);
      ids.push_back(std::string(s,sz));
    }
    delete [] s;    
    return true;
  }

  static bool writeKinWithIDs(const char* outname, MatrixXd& K, double wsum, std::vector<std::string>& ids) {
    wFile wf(outname);

    int i,j;
    //size_t k;
    int n = K.rows();

    //size_t size = N_MAGIC + sizeof(int) + sizeof(int) + sizeof(double) + sizeof(double) * (size_t) n * (size_t) ( n + 1 ) / 2;

    //notice("Allocating a memory of size %llu", size);

    // write the header
    size_t size = N_MAGIC + sizeof(int) + sizeof(int) + sizeof(double);
    char* p = (char*) malloc(size);
    notice("writeKinWithIDs() - Allocating a memory of size %llu", size);    
    for(int i=0; i < N_MAGIC; ++i) { p[i] = magicKin[i]; }
    *(int*)(p + N_MAGIC) = n;
    *(int*)(p + N_MAGIC + sizeof(int)) = 0;
    *(double*)(p + N_MAGIC + sizeof(int) + sizeof(int)) = wsum;
    wf.write(p,size);    

    /*char* p = (char*) malloc(size);
    for(int i=0; i < N_MAGIC; ++i) { p[i] = magicKin[i]; }
    *(int*)(p + N_MAGIC) = n;
    *(int*)(p + N_MAGIC + sizeof(int)) = 0;
    *(double*)(p + N_MAGIC + sizeof(int) + sizeof(int)) = wsum;
    double* a = (double*)(p + N_MAGIC + sizeof(int) + sizeof(int) + sizeof(double));
    */

    free(p);
    size = sizeof(double) * (size_t)(n+1);
    notice("writeKinWithIDs() - Allocating a memory of size %llu", size);
    p = (char*) malloc(size);        
    double* a = (double*)p;
    //for(i=0, k=0; i < n; ++i) {
    for(i=0; i < n; ++i) {      
      for(j=0; j <= i; ++j) { //, ++k) {
	//a[k] = K(j,i);
	a[j] = K(j, i);
      }
      wf.write(p,sizeof(double) * (i+1)); // write kinship coefficient for each row
    }
    //wf.write(p,size);
    free(p);

    writeIDs(wf, ids);

    return true;
  }

  static bool writeTxtKinWithIDs(const char* outname, MatrixXd& K, double wsum, std::vector<std::string>& ids) {
    wFile wf(outname);
    int i,j;
    int n = K.rows();
    for(i=0; i < n; ++i) {
      if ( i > 0 ) wf.printf("\t");
      wf.printf("%s",ids[i].c_str());
    }
    wf.printf("\n");

    if ( fabs(wsum-1) > 1e-3 ) K *= wsum;

    for(i=0; i < n; ++i) {
      wf.printf("%s",ids[i].c_str());
      for(j=0; j < n; ++j) {
	wf.printf("\t%.10lf",K(i,j));
      }
      wf.printf("\n");
    }
    wf.close();
    return true;
  }
  
  static bool readKinWithIDs(const char* inname, MatrixXd& K, int& flag, double& wsum, std::vector<std::string>& ids) {
    pFile tf(inname);
    char buf[N_MAGIC];
    tf.read(buf,N_MAGIC);
    int n, fl;
    double ws;
    if ( strncmp(buf,MAGIC_KIN,N_MAGIC) != 0 ) error("MAGIC number for kinship matrix does not match");
    tf.read(&n,sizeof(int));
    tf.read(&fl,sizeof(int)); flag = fl;
    tf.read(&ws,sizeof(double)); wsum = ws;
    //size_t size = (size_t)n * ( (size_t)n + 1 ) / 2;
    //double* p = new double[size];
    double* p = new double[(size_t)n];

    notice("readKinWithIDs(): Allocating size %zu bytes.. \n", n * sizeof(double));    

    //notice("Allocating size %zu bytes.. \n", size * sizeof(double));
	   
    //size_t nr = tf.read(p, sizeof(double)*size);
    int i, j;
    size_t k;
    
    //if ( nr != sizeof(double)*size ) error("Kinship matrix size is truncated");
    K.resize(n,n);
    for(i=0, k=0; i < n; ++i) {
      // read (i+1) byte at a fime
      size_t nr = tf.read(p, sizeof(double)*(i+1));
      if ( nr != sizeof(double)*(i+1) ) error("Kinship matrix size is truncated at i=%d", i);
      for(j=0; j <= i; ++j, ++k) {
	if ( i == j ) {
	  //K(i,j) = p[k];
	  K(i,j) = p[j];	  
	}
	else {
	  //K(i,j) = K(j,i) = p[k];
	  K(i,j) = K(j,i) = p[j];  
	}
      }
    }
    delete[] p;

    readIDs(tf, n, ids);
    return true;    
  }

  static bool readTxtKinWithIDs(const char* inname, MatrixXd& K, int& flag, double& wsum, std::vector<std::string>& ids) {
    pFile tf(inname);  // read kinship matrix as a plain or gzipped file
    const char* line = NULL;
    std::vector<std::string> tokens;
    // the kinship file format is expected to have a header line (R compatible format)
    for( int i=0; ( line = tf.getLine() ) != NULL; ++i ) {
      if ( i == 0 ) {
	if ( line[0] == '#' ) pFile::tokenizeLine(line+1," \t\r\n",ids);
	else pFile::tokenizeLine(line," \t\r\n",ids);
	K.resize((int)ids.size(),(int)ids.size());
      }
      else {
	pFile::tokenizeLine(line," \t\r\n",tokens);
	// sanity checking
	if ( tokens[0] != ids[i-1] ) {
	  error("The ID %s in line %d in text-format kinship %s does not match to the ID %s in %d-th column in the first line",tokens[0].c_str(), i+1, inname, ids[i-1].c_str(), i);
	}
	if ( tokens.size() != ids.size() + 1 ) {
	  error("The number of columns in text-format kinship %s is inconsistent. Observed %d columns in the first line, and %d columns in line %d", inname, (int)ids.size(), (int) tokens.size(), i+1);
	}
	for( int j=0; j < (int)ids.size(); ++j) {
	  K(i-1,j) = atof(tokens[j+1].c_str());
	  if ( i-1 > j ) {
	    if ( fabs(K(i-1,j)-K(j,i-1)) > 1e-10 ) {
	      error("The kinship matrix is not symmetric. Observed %.10lf in (%d,%d) and %.10lf in (%d,%d), differing by %lg",
		    K(i-1,j), i-1, j, K(j,i-1), j, i-1, fabs(K(i-1,j)-K(j,i-1)));
	    }
	    K(i-1,j) = K(j,i-1);
	  }
	}
      }
    }
    flag = 0;
    wsum = 1.;
    return true;
  }

  static void cov2cor(MatrixXd& m) {
    int nr = m.rows();
    int nc = m.cols();
    if ( nc != nr ) { error("pEmmax::covcor() : input matrix must be square"); }
    
    VectorXd v(nr);
    for(int i=0; i < nr; ++i) {
      v(i) = sqrt(1./m(i,i));
    }

    for(int i=0; i < nr; ++i) {
      for(int j=0; j < nr; ++j) {
	m(i,j) = m(i,j) * v(i) * v(j);
      }
    }
  }

    /*
  static void vectorWiseProd(MatrixXd& m, VectorXd& v, bool rowwise = true) {
    int nv = v.size();
    int nr = m.rows();
    int nc = m.cols();
    if ( rowwise == false) {  // m[i,] *= v (component-wise)
      if ( nc != nv ) { throw pexception; }
      for(int i=0; i < nr; ++i) {
	for(int j=0; j < nc; ++j) {
	  m(i,j) *= v(j);
	}
      }
    }
    else {
      if ( nr != nv ) { throw pexception; }
      for(int i=0; i < nr; ++i) {
	for(int j=0; j < nc; ++j) {
	  m(i,j) *= v(i);
	}
      }
    }
  }

  static void standardize(MatrixXd& m, bool rowwise = true) {
    if ( rowwise ) {
      VectorXd cmean = m.rowwise().mean();
      VectorXd sqsum = m.rowwise().squaredNorm();
      VectorXd csd = (sqsum.array()/m.cols() - cmean.array().square()).sqrt().inverse().matrix();
      m.colwise() -= cmean;
      vectorWiseProd(m,csd,rowwise);
    }
    else {
      VectorXd cmean = m.colwise().mean();
      VectorXd sqsum = m.colwise().squaredNorm();
      VectorXd csd = (sqsum.array()/m.rows() - cmean.array().square()).sqrt().inverse().matrix();
      m.rowwise() -= cmean;
      vectorWiseProd(m,csd,rowwise);
    }
  }
    */

  static bool writeEigenWithIDs(const char* outname, MatrixXd& evecs, VectorXd& evals, double trK, std::vector<std::string>& ids) {
    wFile wf(outname);

    int i,j;
    //size_t k;
    int r = evecs.rows();
    int c = evecs.cols();
    
    size_t size = N_MAGIC + sizeof(int) + sizeof(int) + sizeof(double); // + sizeof(double) * ( (size_t)r * (size_t)c + (size_t)c );
    char* p = (char*) malloc(size);

    // write headers first
    notice("writeEigenWithIDs(): Allocating a memory of size %llu", size);    
    for(int i=0; i < N_MAGIC; ++i) { p[i] = magicEig[i]; }
    *(int*)(p + N_MAGIC) = r;
    *(int*)(p + N_MAGIC + sizeof(int)) = c;
    *(double*)(p + N_MAGIC + sizeof(int) + sizeof(int)) = trK;
    wf.write(p,size);        
    //double* a = (double*)(p + N_MAGIC + sizeof(int) + sizeof(int) + sizeof(double));
    free(p);

    //for(i=0, k=0; i < c; ++i) { // ++k) {

    // write eigenvvalues next
    size = (sizeof(double) * (size_t)c);
    p = (char*) malloc(size);
    double* a = (double*)p;
    notice("writeEigenWithIDs(): Allocating a memory of size %llu", size);        
    for(i=0; i < c; ++i) { 
      //a[k] = evals(k);
      a[i] = evals(i);      
    }
    wf.write(p,size);            
    free(p);    

    // write eigenvvectors next    
    size = (sizeof(double) * (size_t)r);
    p = (char*) malloc(size);
    a = (double*)p;
    notice("writeEigenWithIDs(): Allocating a memory of size %llu", size);          
    for(i=0; i < c; ++i) {
      for(j=0; j < r; ++j) { // ++k) {
      //a[k] = evecs(j,i);
	a[j] = evecs(j,i);	
      }
      wf.write(p,size);                
    }
    //wf.write(p,size);
    free(p);

    writeIDs(wf, ids);
    return true;
  }

  static bool readEigenWithIDs(const char* inname, MatrixXd& evecs, VectorXd& evals, double& trK, std::vector<std::string>& ids) {
    pFile tf(inname);
    char buf[N_MAGIC];
    tf.read(buf,N_MAGIC);
    int i, j, r, c;
    size_t k;
    if ( strncmp(buf,MAGIC_EIG,N_MAGIC) != 0 ) { error("pEmmax::readEigen() - magic number mismatch"); }
    tf.read(&r,sizeof(int));
    tf.read(&c,sizeof(int));
    tf.read(&trK,sizeof(double));
    size_t size = r * c + c;

    notice("Allocating a size %zu bytes\n", size);
    
    double* p = new double[size];
    size_t nr = tf.read(p, sizeof(double)*size);
    if ( nr != (size_t)sizeof(double)*size ) { error("pEmmax::readEigen() - file is truncated. Expected %zu bytes, but read %zu bytes", sizeof(double)*size, nr ); }
    evecs.resize(r,c);
    evals.resize(c);
    for(k=0; k < (size_t)c; ++k) {
      evals(k) = p[k];
    }
    for(i=0; i < c; ++i) {
      for(j=0; j < r; ++j, ++k) {
	evecs(j,i) = p[k];
      }
    }
    delete[] p;

    readIDs(tf, r, ids);

    return true;
  }

  static int computeEigenRestricted(MatrixXd& K, MatrixXd& X, MatrixXd& evecs, VectorXd& evals, double& trK) {
    int n = K.rows();
    int p = X.cols();
    MatrixXd I = MatrixXd::Identity(n,n);

    //error("%d %d %d",n,p,X.rows());
    MatrixXd S;
    if ( p == 0 ) {
      S = I;
    }
    else {
      JacobiSVD<MatrixXd> svd(X, ComputeThinU); // X = UDV'
      // S = I - X %*% solve(t(X)*X) * t(X) = I - UU'
      notice("Computing SVD of covariates (%d x %d) , n = %d",(int)svd.matrixU().rows(), (int)svd.matrixU().cols(),n);
      S = I - svd.matrixU() * svd.matrixU().transpose();
    }
    MatrixXd rK = S * (K + I) * S;
    trK = rK.trace();

    SelfAdjointEigenSolver<MatrixXd> eig(rK);

    VectorXd v = eig.eigenvalues();
    for(int i=0; i < p; ++i) {
      if ( fabs(v[i]) > ZEPS ) {
	error("eigenvalue at index %d is %lf",i,v[i]);
      }
    }

    evals.resize(n-p);
    for(int i=p; i < n; ++i) {
      evals(i-p) = v(i) - 1.;
    }

    evecs = eig.eigenvectors().block(0,p,n,n-p);

    return (n-p);
  }

  static bool exists(const char* filename) {
    FILE* fp = fopen(filename,"rb");
    if ( fp == NULL ) {
      return false;
    }
    else {
      fclose(fp);
      return true;
    }
  }


  static double RLL(double delta, const VectorXd& U, const VectorXd& evals) {
    double lXi = 0, R = 0;
    int c = U.rows();
    for(int i=0; i < c; ++i) {
      lXi += log(evals(i)+delta);
      R += ( U(i)*U(i)/(evals(i) + delta) );
    }
    return ( -0.5*c*log(2*M_PI)-0.5*lXi-.5*c*log(R/c)-c/2. );
  }

  static double dRLL(double delta, const VectorXd& U, const VectorXd& evals) {
    double B = 0, R = 0, E = 0;
    int c = U.rows();
    for(int i=0; i < c; ++i) {
      R += ( U(i)*U(i)/(evals(i) + delta) );
      B += ( U(i)*U(i)/(evals(i) + delta)/(evals(i) + delta) );
      E += ( 1/(evals(i) + delta) );
    }
    return (.5*delta*(c*B/R-E));
  }

  static double dRL(double delta, const VectorXd& U, const VectorXd& evals) {
    double B = 0, R = 0, E = 0;
    int c = U.rows();
    for(int i=0; i < c; ++i) {
      R += ( U(i)*U(i)/(evals(i) + delta) );
      B += ( U(i)*U(i)/(evals(i) + delta)/(evals(i) + delta) );
      E += ( 1/(evals(i) + delta) );
    }
    return (.5*(c*B/R-E));
  }

  static double ddRLL(double delta, const VectorXd& U, const VectorXd& evals) {
    double B = 0, R = 0, T = 0, E = 0, E2 = 0, num = 0, ed = 0;
    int c = U.rows();
    for(int i=0; i < c; ++i) {
      num = U(i)*U(i);
      ed = evals(i)+delta;
      R += ( num/ed );
      B += ( num/(ed*ed) );
      T += ( num/(ed*ed*ed) );
      E += ( 1/ed );
      E2 += ( 1/(ed*ed) );
    }
    return (.5*(c*B/R-E + delta*(c*(2*T/R - B*B/(R*R)) + E2) ));
  }

  static double ddRL(double delta, const VectorXd& U, const VectorXd& evals) {
    double B = 0, R = 0, T = 0, E = 0, E2 = 0, num = 0, ed = 0;
    int c = U.rows();
    for(int i=0; i < c; ++i) {
      num = U(i)*U(i);
      ed = evals(i)+delta;
      R += ( num/ed );
      B += ( num/(ed*ed) );
      T += ( num/(ed*ed*ed) );
      E += ( 1/ed );
      E2 += ( 1/(ed*ed) );
    }
    return (.5*(c*(2*T/R - B*B/(R*R)) + E2) );
  }

  static double bisect(double lo, double hi, double e, const VectorXd& U, const VectorXd& evals) {
    double d = hi - lo;
    double point = lo + d * 0.5;
    double fpoint = dRLL(pow(10.,point),U,evals);
    if ( fabs(d) < e || fpoint == 0.0 ) {
      return pow(10.,point);
    }
    else if ( fpoint < 0 ) {
      return bisect(point, hi, e, U, evals);
    }
    else {
      return bisect(lo, point, e, U, evals);
    }
  }
  
  static REML computeREML(const VectorXd& y, const MatrixXd& evecs, const VectorXd& evals, double trK) {
    VectorXd U = evecs.transpose() * y; // (n-p) * 1 matrix
    return computeREMLU(U, evecs, evals, trK);
  }

  static REML computeREMLU(const VectorXd& U, const MatrixXd& evecs, const VectorXd& evals, double trK) {
    int c = U.rows();

    //std::cout << U << std::endl;
    
    // set up grid search
    double minLogDelta = -5;
    double maxLogDelta = 5;
    int nUnits = 100;
    double binSize = (maxLogDelta - minLogDelta) / nUnits;
    
    std::vector<double> dLLs;
    for(int i=0; i <= nUnits; ++i) {
      double delta = pow(10.,minLogDelta + binSize * i);
      dLLs.push_back(dRLL(delta,U,evals));
      //std::cout << delta << "\t" << dLLs.back() << "\t" << RLL(delta,U,evals) << std::endl;
    }

    // find possible maxima
    double optDelta = pow(10.,maxLogDelta);
    double optLL = RLL(optDelta,U,evals);
    double LL0 = optLL;
    double dLL0 = dRL(optDelta,U,evals);
    double ddLL0 = ddRL(optDelta,U,evals);

    if ( dLLs.front() < 0 ) {
      double d = pow(10.,minLogDelta);
      double ll = RLL(d,U,evals);
      if ( ll > optLL ) {
	optLL = ll;
	optDelta = d;
      }
    }

    for(int i=0; i < nUnits; ++i) {
      if ( ( dLLs[i] > 0-TOL ) && ( dLLs[i+1] < TOL ) && ( dLLs[i] * dLLs[i+1] < 0-ZEPS ) ) {
	double d = bisect( minLogDelta + binSize * i, minLogDelta + binSize * (i+1), TOL, U, evals );
	double ll = RLL(d,U,evals);
	if ( ll > optLL ) {
	  optLL = ll;
	  optDelta = d;
	}
      }
    }

    double vg = 0;
    for(int i=0; i < c; ++i) {
      vg += ( U(i)*U(i)/(evals(i) + optDelta)/c );
    }
    double ve = vg * optDelta;

    double h2 = trK/(c * optDelta+trK);

    return REML(optDelta, h2, vg, ve, optLL, LL0, dLL0, ddLL0);
  }

  static double betacf(double a, double b, double x) {
    int m,m2; 
    double aa,c,d,del,h,qab,qam,qap;
    double FPMIN = 1e-30;
    int MAXIT = 1000;
    double EPS = 1e-10;
    
    qab=a+b; 
    qap=a+1.0; 
    qam=a-1.0; 
    c=1.0; 
    d=1.0-qab*x/qap; if (fabs(d) < FPMIN) d=FPMIN; d=1.0/d;
    h=d; 
    for (m=1;m<=MAXIT;m++) {
      m2=2*m; aa=m*(b-m)*x/((qam+m2)*(a+m2)); 
      d=1.0+aa*d; 
      if (fabs(d) < FPMIN) d=FPMIN; 
      c=1.0+aa/c; 
      if (fabs(c) < FPMIN) c=FPMIN; 
      d=1.0/d; 
      h *= d*c; 
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
      d=1.0+aa*d; 
      if (fabs(d) < FPMIN) d=FPMIN; 
      c=1.0+aa/c; 
      if (fabs(c) < FPMIN) c=FPMIN; 
      d=1.0/d; 
      del=d*c; 
      h *= del; 
      if (fabs(del-1.0) < EPS) break;
    } 
    if (m > MAXIT) {
      error("a or b too big, or MAXIT too small in betacf %lf %lf %lf",a,b,x);
    }
    return h;
  }
  
  static double gammln(double xx) {
    double x,y,tmp,ser; 
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5}; 
    int j;
    y=x=xx; 
    tmp=x+5.5; 
    tmp -= (x+0.5)*log(tmp); 
    ser=1.000000000190015; 
    for (j=0;j<=5;j++) 
      ser += cof[j]/++y; 
    return -tmp+log(2.5066282746310005*ser/x);
  }
  
  static double betai(double a, double b, double x) {
    double bt;
    if (x < 0.0 || x > 1.0) {
      error("Bad x in routine betai"); 
    }
    if (x == 0.0 || x == 1.0) bt=0.0; 
    else bt=exp((gammln(a+b)-gammln(a)-gammln(b))+(a*log(x))+(b*log(1.0-x)));
    if (x < (a+1.0)/(a+b+2.0))
      return bt*betacf(a,b,x)/a; 
    else
      return 1.0-bt*betacf(b,a,1.0-x)/b;
  }
  
  static double tcdf(double t, double nu) {
    if ( std::isnan(t) ) return 1.;
    else return betai(nu/2.,0.5,nu/(nu+t*t));
  }
};

class pEmmax {
 public:
  // below is for member variables and functions
  int flagK, ncovs;
  double wsumK, trK;
  bool normalize, noIntercept;
  VectorXd y, evalR; // phenotypes
  MatrixXd X, evecR, T; // covariates
  MatrixXd K, Ks, P; // kinship matrices
  fVcf tvcf;  // VCF file
  std::vector<std::string> inds;
  std::vector<std::string> pedcols;

 pEmmax() : flagK(0), ncovs(0), wsumK(0), trK(0), normalize(false), noIntercept(false) {}

 // load files from 
 void loadFiles(const char* _phef, const char* _covf, const char* _indf, const char* _kinf, const char* _eigf, const char* _vcf, const char* _region, const char* _rule, const char* _key, bool _pass) {
    //if ( _phef == NULL ) { error("pEmmax::loadFiles() - phef cannot be NULL");
   if ( ( _phef != NULL ) && ( _phef[0] == '\0' ) ) _phef = NULL;
   if ( ( _covf != NULL ) && ( _covf[0] == '\0' ) ) _covf = NULL;
   if ( ( _indf != NULL ) && ( _indf[0] == '\0' ) ) _indf = NULL;
   if ( ( _kinf != NULL ) && ( _kinf[0] == '\0' ) ) _kinf = NULL;
   if ( ( _eigf != NULL ) && ( _eigf[0] == '\0' ) ) _eigf = NULL;
   if ( ( _vcf != NULL ) && ( _vcf[0] == '\0' ) ) _vcf = NULL;
   if ( ( _region != NULL ) && ( _region[0] == '\0' ) ) _region = NULL;
   if ( ( _rule != NULL ) && ( _rule[0] == '\0' ) ) _rule = NULL;
   if ( ( _region != NULL ) && ( _region[0] == '\0' ) ) {
     error("pEmmax::loadFiles() : _key cannot be empty()");
   }
   
   if ( ( _kinf == NULL ) && ( _eigf == NULL ) && ( _vcf == NULL ) ) {
     error("pEmmax::loadFiles() - one of kinf, eigf, vcf cannot be NULL");
   }
   
   std::map<std::string,int> pheorder, covorder, indorder, kinorder, eigorder, idcnts;
   std::vector<std::string> kinids, eigids, tokens;
   std::vector<double> vphe, vcov;
   pFile tphef, tindf, tcovf, teigf;
   int i, j, nfiles = 0;
   int nMissing = 0;
   const char* line = NULL;
   
   if ( _phef != NULL ) {  // if phenotype exists, stored them into pheorder map
     tphef.load(_phef);
     while ( (line = tphef.getLine()) != NULL ) {
       if ( line[0] == '#' )  continue;
       pFile::tokenizeLine(line," \t\r\n", tokens);
       if ( tokens.size() != 2 ) 
	 error("Phenotype file must have two columns");
       pheorder[tokens[0]] = (int)vphe.size();
       ++idcnts[tokens[0]];
       vphe.push_back(atof(tokens[1].c_str()));
     }
     ++nfiles;
   }
   
   if ( _covf == NULL ) { // if covariate exists store them into covorder map
     covorder = pheorder;
     vcov.resize(vphe.size());
     for(i=0; i < (int)vcov.size(); ++i) { vcov[i] = 1.; }
     ncovs = 1;
   }
   else {
     tcovf.load(_covf);
     int r = 0, c = 0;
     while( (line = tcovf.getLine()) != NULL ) {
       if ( line[0] == '#' ) continue;
       pFile::tokenizeLine(line," \t\r\n", tokens);
       covorder[tokens[0]] = r;
       ++idcnts[tokens[0]];
       if ( c == 0 ) {
	 c = (int)tokens.size();
       }
       for(i=1; i < c; ++i) {
	 if ( tokens[i] == "NA" ) {
	   vcov.push_back(std::numeric_limits<double>::quiet_NaN());
	   ++nMissing;
	 }
	 else {
	   vcov.push_back(atof(tokens[i].c_str()));
	 }
       }
       if ( !noIntercept ) vcov.push_back(1.);
       ++r;
     }
     ncovs = noIntercept ? c-1 : c;
     ++nfiles;
   }
   
   if ( _indf != NULL ) {  // if indf exists, indorder is updated
     tindf.load(_indf);
     i = 0;
     while( (line = tindf.getLine()) != NULL ) {
       if ( line[0] == '#' ) continue;
       pFile::tokenizeLine(line," \t\r\n", tokens);
       if ( tokens.size() != 2 ) 
	 error("Individual index file must have two columns");
       indorder[tokens[0]] = i;
       ++idcnts[tokens[0]];
       ++i;
     }
     ++nfiles;
   }

   //notice("goo");
   
   if ( _kinf != NULL ) {  // if kinf exists kinorder is updated
     pEmmaxHelper::readKinWithIDs(_kinf, K, flagK, wsumK, kinids);
     
     notice("Reading kinship matrix %s with %d IDs",_kinf,(int)kinids.size());
     
     if ( wsumK > 1.001 ) {
       if ( normalize ) {
	 K *= (1./wsumK);
	 wsumK = 1.;
       }
       else {
	 warning("Kinship file %s is not normalized",_kinf);
       }
     }
     
     for(i=0; i < (int)kinids.size(); ++i) {
       kinorder[kinids[i]] = i;
       ++idcnts[kinids[i]];
     }
     ++nfiles;
   }

   //notice("foo");

   if ( _eigf != NULL ) {  // if eigf exists eigorder is updated
     notice("Reading eigenvectors");
     pEmmaxHelper::readEigenWithIDs(_eigf, evecR, evalR, trK, eigids);
     for(i=0; i < (int)eigids.size(); ++i) {
       eigorder[eigids[i]] = i;
       ++idcnts[eigids[i]];
     }
     ++nfiles;
   }

   //notice("bar");
   
   // identify matching IDs
   // eigids should be exact overlaps
   inds.clear();
   std::set<std::string> idset;
   for(std::map<std::string,int>::iterator it = idcnts.begin();
       it != idcnts.end(); ++it) {
     if ( it->second == nfiles ) {  // individuals must overlap to all 
       idset.insert(it->first);
       if ( _vcf == NULL ) inds.push_back(it->first); // if vcf does not exist, use alphabetical order
     }
     else if ( it->second > nfiles ) {
       error("Observed unexpected ID overlap counts for %s. Was individual ID duplicated?",it->first.c_str());
     }
   }

   //notice("foo2");
   
   if ( _vcf != NULL ) {  // if VCF exists
     //notice("_rule = %s",_rule);
     tvcf.load(_vcf, _region, _key, _rule, _pass, idset);
     inds = tvcf.inds;  // Respect individual orders in VCF
     // make sure eigf is in the same order with VCF
     if ( _eigf != NULL ) {  // check order of eigen
       if ( ( tvcf.nInds == (int)eigids.size() ) && ( std::equal( tvcf.inds.begin(), tvcf.inds.end(), eigids.begin() ) ) ) {
       }
       else if ( tvcf.nInds == (int)eigids.size() ) { // if size matches
	 int r = evecR.rows();
	 int c = evecR.cols();
	 MatrixXd evecsS(r, c);
	 std::vector<int> ioids;
	 std::map<std::string,int>::iterator it;
	 for(i=0; i < tvcf.nInds; ++i) {
	   it = eigorder.find(tvcf.inds[i]);
	   if ( it == eigorder.end() ) {
	     error("Cannot match consensus individual ID %s with eigf file",tvcf.inds[i].c_str());
	   }
	   ioids.push_back(it->second);
	 }
	 for(i=0; i < r; ++i) { 
	   for(j=0; j < c; ++j) {
	     evecsS(i,j) = evecR(ioids[i],j);
	   }
	 }
	 evecR = evecsS;
       }
       else {
	 error("The number of individuals between --in-eigf %s and --vcf %s are not identical",_eigf,_vcf);
       }
       
       //if ( !std::equal( tvcf.inds.begin(), tvcf.inds.end(), eigids.begin() ) ) {
       //notice("%d %d",(int)tvcf.inds.size(),(int)eigids.size());
       //for(i=0; i < (int)eigids.size(); ++i) {
       //  notice("%d\t%s\t%s",i,tvcf.inds[i].c_str(),eigids[i].c_str());
       //}
       //error("The individual IDs between --in-eigf %s and --vcf %s are not identical",_eigf,_vcf);
       //}
     }

     //notice("bar2");
   
     if ( _kinf != NULL ) {
       if ( ( tvcf.nInds == (int)kinids.size() ) && ( std::equal( tvcf.inds.begin(), tvcf.inds.end(), kinids.begin() ) ) ) {
	 //error("The individual IDs between --in-eigf %s and --kinf %s are not identical",_eigf,_kinf);
	 Ks = K;
       }
       else {
	 notice("Reordering kinship matrix..\n");
	 std::vector<int> ioids;
	 std::map<std::string,int>::iterator it;
	 for(i=0; i < tvcf.nInds; ++i) {
	   it = kinorder.find(tvcf.inds[i]);
	   if ( it == kinorder.end() ) {
	     error("Cannot match consensus individual ID %s with kinship file",tvcf.inds[i].c_str());
	   }
	   ioids.push_back(it->second);
	 }
	 Ks.resize(tvcf.nInds, tvcf.nInds);
	 for(i=0; i < tvcf.nInds; ++i) {
	   for(j=0; j < tvcf.nInds; ++j) {
	     Ks(j,i) = K(ioids[j],ioids[i]);
	   }
	 }
       }
     }
   }
   else if ( _eigf != NULL ) {  // if --in-eigf was given
     if ( idset.size() != eigids.size() )  {
       error("Number of consensus individual does not match between --in-eigf and phenotype files");
     }
     
     inds = eigids;
     
     if ( _kinf != NULL ) {
       if ( std::equal( kinids.begin(), kinids.end(), eigids.begin() ) ) {
	 Ks = K;
       }
       else {
	 std::vector<int> ioids;
	  std::map<std::string,int>::iterator it;
	  for(i=0; i < (int)eigids.size(); ++i) {
	    it = kinorder.find(eigids[i]);
	    if ( it == kinorder.end() ) {
	      error("Cannot match consensus individual ID %s with kinship file",eigids[i].c_str());
	    }
	    ioids.push_back(it->second);
	  }
	  
	  int n = (int)inds.size();
	  
	  Ks.resize(n, n);
	  for(i=0; i < n; ++i) {
	    for(j=0; j < n; ++j) {
	      Ks(j,i) = K(ioids[j],ioids[i]);
	    }
	  }
	}
      }
    }
    else if ( _kinf != NULL ) {  // if --in-eigf was given
      //notice("idset.size() = %d, kinids.size() = %d",(int)idset.size(),(int)kinids.size());
      
      if ( idset.size() == kinids.size() )  {
	Ks = K;
	inds = kinids;
      }
      else {
	// when subsetting the kinship matrix, the order must be preserved

	std::vector<int> ioids;
	inds.clear();
	for(i=0; i < (int)kinids.size(); ++i) {
	  if ( idset.find(kinids[i]) != idset.end() ) {
	    ioids.push_back(i);
	    inds.push_back(kinids[i]);
	  }
	  //notice("%d\t%s",i,kinids[i].c_str());
	}
	/*
	std::set<std::string>::iterator it;
	std::map<std::string,int>::iterator mit;
	for(it = idset.begin(); it != idset.end(); ++it) {
	  mit = kinorder.find(*it);
	  if ( mit == kinorder.end() ) {
	    error("Cannot match consensus individual ID %s with kinship file",it->c_str());
	  }
	  ioids.push_back(mit->second);
	}
	*/

	int n = idset.size();
	if ( n != (int)ioids.size() ) {
	  error("Number of individuals are inconsistent between files (%d vs %d). Is any individual ID duplicated?",n,(int)ioids.size());
	}
	
	Ks.resize(n, n);
	for(i=0; i < n; ++i) {
	  for(j=0; j < n; ++j) {
	    Ks(j,i) = K(ioids[j],ioids[i]);
	  }
	}
      }
    }

    if ( _phef != NULL ) {
      std::map<std::string,int>::iterator it;

      int n = (_vcf != NULL) ? tvcf.nInds : (int)inds.size();
      
      y.resize(n);
      X.resize(n,ncovs);
      
      for(i=0; i < n; ++i) {
	it = pheorder.find(inds[i]);
	if ( it == pheorder.end() ) {
	  error("Cannot match consensus individual ID %s with kinship file",inds[i].c_str());
	}
	y(i) = vphe[it->second];
	for(j=0; j < ncovs; ++j) {
	  X(i,j) = vcov[it->second * ncovs + j];
	}
      }
      
      // resolve missing
      if ( nMissing > 0 ) {
	notice("Detected missing covariates... Using sample average to resolve %d missing covariates",nMissing);
	int r = X.rows();
	int c = X.cols();
	
	std::vector<int> colCnts(c,0);
	std::vector<double> colSums(c,0.);
	
	for(int i=0; i < r; ++i) {
	  for(int j=0; j < c; ++j) {
	    if ( !std::isnan(X(i,j)) ) {
	      ++colCnts[j];
	      colSums[j] += X(i,j);
	    }
	  }
	}
	
	for(int i=0; i < r; ++i) {
	  for(int j=0; j < c; ++j) {
	    if ( std::isnan(X(i,j)) ) {
	      X(i,j) = (colCnts[j] > 0) ? (colSums[j]/colCnts[j]) : 0;
	    }
	  }
	}
      }
      notice("Finished loading phenotypes across %d individuals with %d covariates (including interception)", n, ncovs);
    }
    else if ( _covf != NULL ) {
      error("--covf option is given without --phenof");
    }
 }
};

class pEmmaxMulti {
 public:
  // below is for member variables and functions
  int flagK, ncovs, nphes;
  double wsumK, trK;
  bool normalize, noIntercept;
  VectorXd evalR;    // eigenvalues
  MatrixXd Y, X, evecR, T; // covariates
  MatrixXd K, P; // kinship matrices
  fVcf tvcf;  // VCF file
  std::vector<std::string> inds;
  std::vector<std::string> pedcols;

 pEmmaxMulti() : flagK(0), ncovs(0), nphes(0), wsumK(0), trK(0), normalize(false), noIntercept(false) {}

 void loadFiles(const char* _phef, const char* _covf, const char* _indf, const char* _kinf, const char* _eigf, const char* _vcf, const char* _region, const char* _rule, const char* _key, bool _pass) {
    //if ( _phef == NULL ) { error("pEmmax::loadFiles() - phef cannot be NULL");
   if ( ( _phef != NULL ) && ( _phef[0] == '\0' ) ) _phef = NULL;
   if ( ( _covf != NULL ) && ( _covf[0] == '\0' ) ) _covf = NULL;
   if ( ( _indf != NULL ) && ( _indf[0] == '\0' ) ) _indf = NULL;
   if ( ( _kinf != NULL ) && ( _kinf[0] == '\0' ) ) _kinf = NULL;
   if ( ( _eigf != NULL ) && ( _eigf[0] == '\0' ) ) _eigf = NULL;
   if ( ( _vcf != NULL ) && ( _vcf[0] == '\0' ) ) _vcf = NULL;
   if ( ( _region != NULL ) && ( _region[0] == '\0' ) ) _region = NULL;
   if ( ( _rule != NULL ) && ( _rule[0] == '\0' ) ) _rule = NULL;
   if ( ( _region != NULL ) && ( _region[0] == '\0' ) ) {
     error("pEmmaxMulti::loadFiles() : _key cannot be empty()");
   }

    if ( ( _kinf == NULL ) && ( _eigf == NULL ) && ( _vcf == NULL ) ) {
      error("pEmmaxMulti::loadFiles() - one of kinf, eigf, vcf cannot be NULL");
    }
    
    std::map<std::string,int> pheorder, covorder, indorder, kinorder, eigorder, idcnts;
    std::vector<std::string> kinids, eigids, tokens;
    std::vector<double> vphe, vcov;
    pFile tphef, tindf, tcovf, teigf;
    int i, j, nfiles = 0;
    int nMissing = 0;
    const char* line = NULL;
  
    if ( _phef != NULL ) {  // fill pheorder if phe exists
      tphef.load(_phef);
      int r = 0, c = 0;
      while( (line = tphef.getLine()) != NULL ) {
	if ( line[0] == '#' ) {
	  pFile::tokenizeLine(line+1," \t\r\n", tokens);
	  for(i=1; i < (int)tokens.size(); ++i) {
	    pedcols.push_back(tokens[i]);
	  }
	  continue;
	}
	pFile::tokenizeLine(line," \t\r\n", tokens);
	pheorder[tokens[0]] = r;
	++idcnts[tokens[0]];
	if ( c == 0 ) {
	  c = (int)tokens.size();
	}
	for(i=1; i < c; ++i) {
	  if ( tokens[i] == "NA" ) {
	    vphe.push_back(std::numeric_limits<double>::quiet_NaN());
	    ++nMissing;
	  }
	  else {
	    vphe.push_back(atof(tokens[i].c_str()));
	  }
	}
	++r;
      }
      nphes = c-1;
      ++nfiles;
    }

    if ( _covf == NULL ) { // fill covorder if exist, otherwise identical to pheorder
      covorder = pheorder;
      vcov.resize(vphe.size());
      for(i=0; i < (int)vcov.size(); ++i) { vcov[i] = 1.; }
      ncovs = 1;
    }
    else {
      tcovf.load(_covf);
      int r = 0, c = 0;
      while( (line = tcovf.getLine()) != NULL ) {
	if ( line[0] == '#' ) continue;
	pFile::tokenizeLine(line," \t\r\n", tokens);
	covorder[tokens[0]] = r;
	++idcnts[tokens[0]];
	if ( c == 0 ) {
	  c = (int)tokens.size();
	}
	for(i=1; i < c; ++i) {
	  if ( tokens[i] == "NA" ) {
	    vcov.push_back(std::numeric_limits<double>::quiet_NaN());
	    ++nMissing;
	  }
	  else {
	    vcov.push_back(atof(tokens[i].c_str()));
	  }
	}
	if ( !noIntercept ) vcov.push_back(1.);
	++r;
      }
      ncovs = noIntercept ? c-1 : c;
      ++nfiles;
    }

    if ( _indf != NULL ) {  // fill indorder if exists
      tindf.load(_indf);
      i = 0;
      while( (line = tindf.getLine()) != NULL ) {
	if ( line[0] == '#' ) continue;
	pFile::tokenizeLine(line," \t\r\n", tokens);
	if ( tokens.size() != 2 ) 
	  error("Individual index file must have two columns");
	indorder[tokens[0]] = i;
	++idcnts[tokens[0]];
	++i;
      }
      ++nfiles;
    }

    //notice("bar");

    if ( _kinf != NULL ) { // kill inorder if exists
      pEmmaxHelper::readKinWithIDs(_kinf, K, flagK, wsumK, kinids);

      notice("Reading kinship matrix %s with %d IDs",_kinf,(int)kinids.size());

      if ( wsumK > 1.001 ) {
	if ( normalize ) {
	  K *= (1./wsumK);
	  wsumK = 1.;
	}
	else {
	  warning("Kinship file %s is not normalized",_kinf);
	}
      }

      for(i=0; i < (int)kinids.size(); ++i) {
	kinorder[kinids[i]] = i;
	++idcnts[kinids[i]];
      }
      ++nfiles;
    }

    //notice("foo");

    if ( _eigf != NULL ) { // read eigorder if exists
      pEmmaxHelper::readEigenWithIDs(_eigf, evecR, evalR, trK, eigids);
      for(i=0; i < (int)eigids.size(); ++i) {
	eigorder[eigids[i]] = i;
	++idcnts[eigids[i]];
      }
      ++nfiles;
    }

    //notice("moo");

    // identify matching IDs
    // eigids should be exact overlaps
    inds.clear();
    std::set<std::string> idset;
    for(std::map<std::string,int>::iterator it = idcnts.begin();
	it != idcnts.end(); ++it) {
      if ( it->second == nfiles ) {
	idset.insert(it->first);
	if ( _vcf == NULL ) inds.push_back(it->first);
      }
      else if ( it->second > nfiles ) {
	error("Observed unexpected ID overlap counts for %s. Was individual ID duplicated?",it->first.c_str());
      }
    }

    //notice("noo");

    if ( _vcf != NULL ) {
      //notice("_rule = %s",_rule);
      tvcf.load(_vcf, _region, _key, _rule, _pass, idset);
      inds = tvcf.inds;  // assume that order of individuals are preserved
      // make sure eigf is in the same order with VCF
      if ( _eigf != NULL ) {
	if ( ( tvcf.nInds == (int)eigids.size() ) && ( std::equal( tvcf.inds.begin(), tvcf.inds.end(), eigids.begin() ) ) ) {
	}
	/*
	else if ( tvcf.nInds == (int)eigids.size() ) {
	  int r = evecR.rows();
	  int c = evecR.cols();
	  MatrixXd evecsS(r, c);
	  std::vector<int> ioids;
	  std::map<std::string,int>::iterator it;
	  for(i=0; i < tvcf.nInds; ++i) {
	    it = eigorder.find(tvcf.inds[i]);
	    if ( it == eigorder.end() ) {
	      error("Cannot match consensus individual ID %s with eigf file",tvcf.inds[i].c_str());
	    }
	    ioids.push_back(it->second);
	  }
	  for(i=0; i < r; ++i) {
	    for(j=0; j < c; ++j) {
	      evecsS(i,j) = evecR(ioids[i],j);
	    }
	  }
	  evecR = evecsS;
	}
	*/
	else {
	  error("The number of individuals between --in-eigf %s and --vcf %s are not identical",_eigf,_vcf);
	}
      }
      if ( _kinf != NULL ) {
	//notice("goo");
	if ( ( tvcf.nInds == (int)kinids.size() ) && ( std::equal( tvcf.inds.begin(), tvcf.inds.end(), kinids.begin() ) ) ) {
	  //Ks = K;
	}
	else {
	  notice("Reordering kinship matrix..\n");
	  MatrixXd Ks(tvcf.nInds, tvcf.nInds);

	  std::vector<int> ioids;
	  std::map<std::string,int>::iterator it;
	  for(i=0; i < tvcf.nInds; ++i) {
	    it = kinorder.find(tvcf.inds[i]);
	    if ( it == kinorder.end() ) {
	      error("Cannot match consensus individual ID %s with kinship file",tvcf.inds[i].c_str());
	    }
	    ioids.push_back(it->second);
	  }
	  for(i=0; i < tvcf.nInds; ++i) {
	    for(j=0; j < tvcf.nInds; ++j) {
	      Ks(j,i) = K(ioids[j],ioids[i]);
	    }
	  }

	  K = Ks;
	  //error("The individual IDs between VCF %s and --kinf %s are not identical",_vcf,_kinf);
	}
      }
    }
    else if ( _eigf != NULL ) {  // if --in-eigf was given
      if ( idset.size() != eigids.size() )  {
	error("Number of consensus individual does not match between --in-eigf and phenotype files");
      }

      inds = eigids;
      
      if ( _kinf != NULL ) {
	if ( std::equal( kinids.begin(), kinids.end(), eigids.begin() ) ) {
	  //Ks = K;
	}
	else {
	  std::vector<int> ioids;
	  std::map<std::string,int>::iterator it;
	  for(i=0; i < (int)eigids.size(); ++i) {
	    it = kinorder.find(eigids[i]);
	    if ( it == kinorder.end() ) {
	      error("Cannot match consensus individual ID %s with kinship file",eigids[i].c_str());
	    }
	    ioids.push_back(it->second);
	  }
	  
	  int n = (int)inds.size();
	  
	  MatrixXd Ks(n, n);
	  for(i=0; i < n; ++i) {
	    for(j=0; j < n; ++j) {
	      Ks(j,i) = K(ioids[j],ioids[i]);
	    }
	  }
	  K = Ks;
	  //error("The individual IDs between --eigf %s and --kinf %s are not identical",_eigf,_kinf);
	}
      }
    }
    else if ( _kinf != NULL ) {  // 
      if ( idset.size() == kinids.size() )  {
	//Ks = K;
	inds = kinids;
      }
      else {
	std::vector<int> ioids;
	inds.clear();
	for(i=0; i < (int)kinids.size(); ++i) {
	  if ( idset.find(kinids[i]) != idset.end() ) {
	    ioids.push_back(i);
	    inds.push_back(kinids[i]);
	  }
	}

	int n = idset.size();
	if ( n != (int)ioids.size() ) {
	  error("Number of individuals are inconsistent between files (%d vs %d). Is any individual ID duplicated?",n,(int)ioids.size());
	}
	
	MatrixXd Ks(n, n);
	for(i=0; i < n; ++i) {
	  for(j=0; j < n; ++j) {
	    Ks(j,i) = K(ioids[j],ioids[i]);
	  }
	}
	K = Ks;
	//error("The individual IDs between --VCF %s and --kinf %s are not identical",_vcf,_kinf);
      }
    }

    if ( _phef != NULL ) {
      std::map<std::string,int>::iterator it;

      int n = (_vcf != NULL) ? tvcf.nInds : (int)inds.size();
      
      Y.resize(n,nphes);
      X.resize(n,ncovs);
      
      for(i=0; i < n; ++i) {
	it = pheorder.find(inds[i]);
	if ( it == pheorder.end() ) {
	  error("Cannot match consensus individual ID %s with kinship file",inds[i].c_str());
	}
	for(j=0; j < nphes; ++j) {
	  Y(i,j) = vphe[it->second * nphes + j];
	}
	for(j=0; j < ncovs; ++j) {
	  X(i,j) = vcov[it->second * ncovs + j];
	}
      }
      
      // resolve missing
      if ( nMissing > 0 ) {
	notice("Detected missing covariates... Using sample average to resolve %d missing covariates",nMissing);
	int r = X.rows();
	int c = X.cols();
	
	std::vector<int> colCnts(c,0);
	std::vector<double> colSums(c,0.);
	
	for(int i=0; i < r; ++i) {
	  for(int j=0; j < c; ++j) {
	    if ( !std::isnan(X(i,j)) ) {
	      ++colCnts[j];
	      colSums[j] += X(i,j);
	    }
	  }
	}
	
	for(int i=0; i < r; ++i) {
	  for(int j=0; j < c; ++j) {
	    if ( std::isnan(X(i,j)) ) {
	      X(i,j) = (colCnts[j] > 0) ? (colSums[j]/colCnts[j]) : 0;
	    }
	  }
	}
      }
      notice("Finished loading phenotypes across %d individuals with %d covariates (including interception)", n, ncovs);
    }
    else if ( _covf != NULL ) {
      error("--covf option is given without --phenof");
    }
 }
};


const char* pEmmaxHelper::magicKin = "EMMA_KIN";
const char* pEmmaxHelper::magicEig = "EMMA_EIG";

#endif // __EIGEN_HELPER_H
