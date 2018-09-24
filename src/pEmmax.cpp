#include "Parameters.h"
#include "Error.h"
#include "fVcf.h"
#include "pEmmax.h"
#include "PhredHelper.h"

#include <iostream>
#include <string>
#include <getopt.h>
#include <ctime>
#include <algorithm>
#include <set>
#include <cmath>

#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Cholesky"

#define DS_THRES 0.5

using namespace Eigen;

// Options for generating kinship matrix
int runGenKin(int argc, char** argv) {
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_INTPARAMETER("unit",&arg.unit)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)

    LONG_PARAMETER_GROUP("Filtering Options")
    LONG_STRINGPARAMETER("weight",&arg.weight)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)

    LONG_PARAMETER_GROUP("Empircial Kinship Method")
    EXCLUSIVE_PARAMETER("BN",&arg.BN)
    EXCLUSIVE_PARAMETER("IBS",&arg.IBS)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out-kinf",&arg.kinf)
    LONG_PARAMETER("raw",&arg.raw)
    LONG_PARAMETER("adj-diag",&arg.adjDiag)
    LONG_PARAMETER("cov",&arg.cov)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  //int argstart = pl.ReadWithTrailer(argc, argv) + 1;
  pl.Read(argc, argv);
  pl.Status();

  if ( arg.vcf.empty() || arg.kinf.empty() ) {
    error("ERROR: --vcf and --out-kinf parameters are required");
  }

  pEmmax emx;
  emx.loadFiles(NULL, NULL, arg.indf.c_str(), NULL, NULL, arg.vcf.c_str(), arg.region.c_str(), arg.rule.c_str(), arg.field.c_str(), !arg.ignoreFilter);

  notice("nInds = %d : %s ...",emx.tvcf.nInds,emx.tvcf.inds[0].c_str());

  // read (#UNIT) markers
  int offset, m, i, j;
  double mu, sigma;
  MatrixXd genos;
  VectorXd adjs;
  for(int M=0; emx.tvcf.readMarkers(arg.unit); ) {
    if ( M == 0 ) {
      emx.K = MatrixXd::Zero( emx.tvcf.nInds, emx.tvcf.nInds );
      genos.resize( emx.tvcf.nInds, arg.unit );
      adjs = VectorXd::Zero(arg.unit);
    }
    M += emx.tvcf.nMarkers;

    notice("Processing %d markers across %d individuals...", M, emx.tvcf.nInds);

    for(i=0, m=0; i < emx.tvcf.nMarkers; ++i) {
      if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) &&
           ( emx.tvcf.callRate(i) >= arg.minCallRate ) &&
           ( ( arg.minRSQ == 0 ) || ( emx.tvcf.RSQ(i) >= arg.minRSQ ) ) ) {
        offset = i * emx.tvcf.nInds;
        mu = emx.tvcf.alleleFreq(i)*2.;
        sigma = sqrt(mu * (2. - mu) / 2.);
        if ( arg.BN ) {
          for(j=0; j < emx.tvcf.nInds; ++j) {
            if ( std::isnan(emx.tvcf.genos[j + offset]) ) {
              genos(j,m) = 0;
            }
            else {
              genos(j,m) = (emx.tvcf.genos[j + offset] - mu)/sigma;
            }
          }

          if ( arg.adjDiag ) {
            adjs(m) = (mu-1.)/sigma;
          }
        }
        else if ( arg.IBS ) {
          for(j=0; j < emx.tvcf.nInds; ++j) {
            if ( std::isnan(emx.tvcf.genos[j + offset]) ) {
              genos(j,m) = (mu - 1);
            }
            else {
              genos(j,m) = (emx.tvcf.genos[j + offset] - 1.);
            }
          }
        }
        else {
          // not implemented yet
          error("Kinship option is only available for --BN or --IBS\n");
        }
        ++m;
      }
    }

    fprintf(stderr,"%d markers included\n",m);

    if ( arg.BN ) {
      emx.wsumK += m;
      emx.K += ( genos.block(0,0,emx.tvcf.nInds,m) * genos.block(0,0,emx.tvcf.nInds,m).transpose() );
      if ( arg.adjDiag ) {
        emx.K += MatrixXd(genos.block(0,0,emx.tvcf.nInds,m) * adjs.block(0,0,m,1)).asDiagonal();
      }
    }
    else if ( arg.IBS ) {
      emx.wsumK += (2*m);
      emx.K += ( genos.block(0,0,emx.tvcf.nInds,m) * genos.block(0,0,emx.tvcf.nInds,m).transpose() );
    }
  }

  if ( arg.IBS ) {
    emx.K.array() += (emx.wsumK/2.);
  }

  if ( arg.raw == false ) {
    emx.K *= (1./emx.wsumK);
    emx.wsumK = 1;
    if ( arg.cov == false ) pEmmaxHelper::cov2cor(emx.K);
  }

  if ( emx.wsumK > 0 ) {
    notice("Started writing kinship matrix into %s",arg.kinf.c_str());
    pEmmaxHelper::writeKinWithIDs(arg.kinf.c_str(), emx.K, emx.wsumK, emx.tvcf.inds);
    notice("Finished writing kinship matrix into %s", arg.kinf.c_str());
  }
  else {
    notice("No SNPs were included in kinship calculation. Skipping writing kinship matrix...");
  }

  return 0;
}

// ln(A+B+C) = ln(exp(lnA) + exp(lnB) + exp(lnC)) = lnA + ln(1 + exp(lnB-lnA) + exp(lnC-lnA));
double logAdd3(double a, double b, double c) {
  double M = (a > b) ? (a > c ? a : c) : (b > c ? b : c);
  return ( M + log(exp(a-M) + exp(b-M) + exp(c-M)) );
}

// ln(A+B+C) = ln(exp(lnA) + exp(lnB) + exp(lnC)) = lnA + ln(1 + exp(lnB-lnA) + exp(lnC-lnA));
double logAdd8(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7) {
  double M = a0;
  if ( M < a1 ) M = a1;
  if ( M < a2 ) M = a2;
  if ( M < a3 ) M = a3;
  if ( M < a4 ) M = a4;
  if ( M < a5 ) M = a5;
  if ( M < a6 ) M = a6;
  if ( M < a7 ) M = a7;
  return ( M + log(exp(a0-M) + exp(a1-M) + exp(a2-M) + exp(a3-M) + exp(a4-M) + exp(a5-M) + exp(a6-M) + exp(a7-M) ) );
}

// ln(A+B+C) = ln(exp(lnA) + exp(lnB) + exp(lnC)) = lnA + ln(1 + exp(lnB-lnA) + exp(lnC-lnA));
double logAdd9(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8) {
  double M = a0;
  if ( M < a1 ) M = a1;
  if ( M < a2 ) M = a2;
  if ( M < a3 ) M = a3;
  if ( M < a4 ) M = a4;
  if ( M < a5 ) M = a5;
  if ( M < a6 ) M = a6;
  if ( M < a7 ) M = a7;
  if ( M < a8 ) M = a8;
  return ( M + log(exp(a0-M) + exp(a1-M) + exp(a2-M) + exp(a3-M) + exp(a4-M) + exp(a5-M) + exp(a6-M) + exp(a7-M) + exp(a8-M)) );
}

int runMergeKin(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Options")
    LONG_PARAMETER("raw",&arg.raw)
    LONG_PARAMETER("cov",&arg.cov)

    LONG_PARAMETER_GROUP("Empircial Kinship Method")
    EXCLUSIVE_PARAMETER("BN",&arg.BN)
    EXCLUSIVE_PARAMETER("IBS",&arg.IBS)
    EXCLUSIVE_PARAMETER("YV",&arg.YV)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out-kinf",&arg.kinf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  int argstart = pl.ReadWithTrailer(argc, argv) + 1;
  pl.Status();

  int ntrails = argc - argstart;
  argv += argstart;

  if ( ntrails == 0 )     
    error("ERROR: merge command must accompany multiple kinship files at the end");
  else if ( ntrails == 1 ) 
    warning("WARNING: only 1 kinship file is being merged");

  notice("Merging the following kinship matrices");
  for(int i=0; i < ntrails; ++i) {
    arg.kinfs.push_back(argv[i]);
    printf("%s\n",argv[i]);
  }

  if ( arg.kinf.empty() ) {
    error("ERROR: --out-kinf parameters is required");
  }

  MatrixXd Ks, Ki;
  double wss = 0, wsi = 0;
  int flag;

  std::vector<std::string> ids, tmp_ids;
  bool firstKin = true;

  // read multiple kinship files
  for(int i=0; i < (int)arg.kinfs.size(); ++i) {
    notice("Merging kinship file %s..", arg.kinfs[i].c_str());
    if ( pFile::fileType(arg.kinfs[i].c_str()) < 0 ) {
      warning("Kinship file %s does NOT exist. Skipping...",arg.kinfs[i].c_str());
      continue;
    }

    pEmmaxHelper::readKinWithIDs(arg.kinfs[i].c_str(), Ki, flag, wsi, tmp_ids);

    // check whether ID matches or not
    if ( ids.empty() ) { ids = tmp_ids; }
    else if ( ! std::equal( ids.begin(), ids.end(), tmp_ids.begin() ) ) {
      error("The individual IDs between the kinship matrices are not identical to the previous files in kinship file %s", arg.kinfs[i].c_str());
    }

    if ( wsi == 0 ) {
      warning("Kinship file %s has only wsum = 0, meaning it is empty. It won't be used in the analysis",arg.kinfs[i].c_str());
      continue;
    }
    else if ( wsi < 1.001 ) {
      warning("Kinship file %s has only wsum = %lf. It may already have been normalized. Use --raw option for merging multiple kinship matrices",arg.kinfs[i].c_str());
    }
    if ( firstKin ) {
      //notice("foo");
      wss = wsi;
      Ks = Ki;
      firstKin = false;
      //notice("bar");
    }
    else {
      wss += wsi;
      Ks += Ki;
    }
  }

  if ( !arg.raw ) {
    Ks *= (1./wss);
    wss = 1.;
    if ( !arg.cov ) {
      pEmmaxHelper::cov2cor(Ks);
    }
  }

  pEmmaxHelper::writeKinWithIDs(arg.kinf.c_str(), Ks, wss, ids);
  // std::cout << K << std::endl; // for debugging
  return 0;
}

// Perform REML fitting
int runReml(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("covf",&arg.covf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_PARAMETER("no-intercept",&arg.noIntercept)
    LONG_STRINGPARAMETER("kinf",&arg.kinf)
    LONG_PARAMETER("normalize",&arg.normalize)
    LONG_STRINGPARAMETER("in-eigf",&arg.ineigf)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out-eigf",&arg.outeigf)
    LONG_STRINGPARAMETER("out-remlf",&arg.remlf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.phenof.empty() || arg.kinf.empty() || arg.remlf.empty() ) {
    error("--phenof, --kinf, and --out-remlf are required parameters");
  }

  pEmmax emx;
  emx.normalize = arg.normalize;
  emx.loadFiles(arg.phenof.c_str(), arg.covf.c_str(), arg.indf.c_str(), arg.kinf.c_str(), arg.ineigf.c_str(), NULL, NULL, NULL, NULL, true);

  // read/compute eigR and save if needed
  // if eigen file exists, read them from file
  if ( !arg.ineigf.empty() && pEmmaxHelper::exists(arg.ineigf.c_str()) ) { 
    // eigendecomposition should is assumed to be loaded
    notice("Loading eigendecomposition matrix %s",arg.ineigf.c_str());
  }
  else {
    notice("Computing eigendecomposition of residual kinship matrix with size %d x %d, predictor matrix %d x %d",emx.Ks.rows(),emx.Ks.cols(),emx.X.rows(),emx.X.cols());
    //std::cout << X << std::endl;
    pEmmaxHelper::computeEigenRestricted(emx.Ks, emx.X, emx.evecR, emx.evalR, emx.trK); // compute eigL
    if ( ! arg.outeigf.empty() ) {
      notice("Writing eigendecomposition of kinship matrix",arg.outeigf.c_str());
      pEmmaxHelper::writeEigenWithIDs(arg.outeigf.c_str(), emx.evecR, emx.evalR, emx.trK, emx.inds); // save eigL
    }
    //abort();
  }

  // compute REML variance parameters for one variance component
  notice("Computing REML");
  REML reml = pEmmaxHelper::computeREML(emx.y, emx.evecR, emx.evalR, emx.trK);

  wFile wf(arg.remlf.c_str());
  wf.printf("delta\t%lf\n",reml.delta);
  wf.printf("vg\t%lf\n",reml.vg);
  wf.printf("ve\t%lf\n",reml.ve);
  wf.printf("LLK0\t%lf\n",reml.LLK0);
  wf.printf("LLK1\t%lf\n",reml.LLK1);
  wf.printf("h2\t%lf\n",reml.h2);
  return 0;
}

// y ~ X + Z + e
// y(t) ~ Xt + Zt + et
int runAssoc(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_INTPARAMETER("unit",&arg.unit)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
    LONG_INTPARAMETER("minMAC",&arg.minMAC)
    LONG_INTPARAMETER("maxMAC",&arg.maxMAC)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("eigf",&arg.ineigf)
    LONG_STRINGPARAMETER("remlf",&arg.remlf)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out-assocf",&arg.assocf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.phenof.empty() || arg.ineigf.empty() || arg.remlf.empty() ) {
    error("--phenof, --out-remlf and --kinf are required parameters");
  }

  pEmmax emx;
  emx.normalize = arg.normalize;
  emx.loadFiles(arg.phenof.c_str(), arg.covf.c_str(), arg.indf.c_str(), NULL, arg.ineigf.c_str(), arg.vcf.c_str(), arg.region.c_str(), arg.rule.c_str(), arg.field.c_str(), !arg.ignoreFilter);

  pFile treml(arg.remlf.c_str());
  double delta = 0;
  std::vector<std::string> tokens;
  const char* line = NULL;
  while( (line = treml.getLine()) != NULL ) {
    pFile::tokenizeLine(line," \t\r\n", tokens);
    if ( tokens.size() != 2 )
      error("Phenotype file must have two columns");
    delta = atof(tokens[1].c_str());
    notice("Reading delta = %lf",delta);
    break;
  }

  int c = emx.evalR.rows();
  int n = emx.inds.size();
  emx.T.resize(c,n);
  for(int i=0; i < c; ++i) {
    for(int j=0; j < n; ++j) {
      emx.T(i,j) = emx.evecR(j,i)/sqrt(emx.evalR(i)+delta);
    }
  }

  bool binaryFlag = true;
  for(int i=0; i < n; ++i) {
    if ( ( fabs(emx.y(i) - 1) > pEmmaxHelper::ZEPS ) && ( fabs(emx.y(i) - 2) > pEmmaxHelper::ZEPS ) ) {
      binaryFlag = false;
    }
  }
  std::vector<bool> isCases;
  if ( binaryFlag ) {
    for(int i=0; i < n; ++i) {
      isCases.push_back(emx.y(i)>1.5);
    }
  }

  wFile wf(arg.assocf.c_str());
  double mu = 0;
  int i,m,j;

  //wf.printf("#MARKER_ID\tNS\tAC\tCALLRATE\tGENOCNT\tMAF\tSTAT\tPVALUE\tBETA\tSEBETA\tR2\n");
  VectorXd yt = emx.T * emx.y; // (n-p) * 1 matrix
  double sy = 0, syy = 0, sx = 0, sxx = 0, sxy = 0, r = 0, t = 0, beta = 0, sebeta = 0, pval = 0, varE = 0;
  //double ysum = 0, ysq = 0, xsum = 0, xsq = 0, xy = 0, r = 0, t = 0, beta = 0, sebeta = 0, sdx = 0, sdy = 0, pval = 0;
  for(int i=0; i < c; ++i) {
    sy += yt(i);
    syy += (yt(i)*yt(i));
    //ysum += yt(i);
    //ysq += yt(i)*yt(i);
  }

  int genocnts[3] = {0,0,0};
  int casectrlcnts[6] = {0,0,0,0,0,0};
  VectorXd x(n), xt(n);
  if ( binaryFlag ) {
    wf.printf("#CHROM\tBEG\tEND\tMARKER_ID\tNS\tAC\tCALLRATE\tGENOCNT\tMAF\tSTAT\tPVALUE\tBETA\tSEBETA\tR2\tCTRLCNT\tCASECNT\n");
  }
  else {
    wf.printf("#CHROM\tBEG\tEND\tMARKER_ID\tNS\tAC\tCALLRATE\tGENOCNT\tMAF\tSTAT\tPVALUE\tBETA\tSEBETA\tR2\n");
  }

  for(int M=0; emx.tvcf.readMarkers(arg.unit); ) {
    M += emx.tvcf.nMarkers;

    notice("Processing %d markers across %d individuals...", M, emx.tvcf.nInds);
    for(i=0, m=0; i < emx.tvcf.nMarkers; ++i) {
      if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) &&
           ( emx.tvcf.MAF(i) <= arg.maxMAF ) &&
           ( emx.tvcf.MAC(i) >= arg.minMAC ) &&
           ( emx.tvcf.MAC(i) <= arg.maxMAC ) &&
           ( emx.tvcf.callRate(i) >= arg.minCallRate ) &&
           ( ( arg.minRSQ == 0 ) || ( emx.tvcf.RSQ(i) >= arg.minRSQ ) ) ) {
        //fprintf(stderr,"%d\n",i);
        int offset = i * emx.tvcf.nInds;
        mu = emx.tvcf.alleleFreq(i)*2.;
        for(j=0; j < emx.tvcf.nInds; ++j) {
          if ( std::isnan(emx.tvcf.genos[j + offset]) ) {
            x(j) = mu;
          }
          else {
            x(j) = emx.tvcf.genos[j + offset];
          }
        }

        xt = emx.T * x;
        sx = sxx = sxy = 0;
        //xy = xsum = xsq = 0;
        for(j=0; j < c; ++j) {
          //xy += xt(j)*yt(j);
          //xsum += xt(j);
          //xsq += xt(j)*xt(j);
          sxy += (xt(j)*yt(j));
          sx += xt(j);
          sxx += (xt(j)*xt(j));
        }

        // this is the estimates including the intercept always

        beta = ((c+1)*sxy-sx*sy)/((c+1)*sxx-sx*sx);
        varE = 1/(c+1.)/(c-1.)*((c+1)*syy-sy*sy-beta*beta*((c+1)*sxx-sx*sx));
        sebeta = sqrt((c+1)*varE/((c+1)*sxx-sx*sx));
        r = ((c+1)*sxy-sx*sy)/sqrt(((c+1)*sxx-sx*sx)*((c+1)*syy-sy*sy));
        t = r * sqrt((c-1)/(1-r*r+pEmmaxHelper::ZEPS));
        pval = pEmmaxHelper::tcdf(t, c-1);

        /*
        beta = sxy/sxx;
        sebeta = (syy*sxx-sxy*sxy)/(c*sxx*sxx);
        r = sxy/sqrt(sxx*sxy);
        t = r * sqrt((c-1)/(1-r*r+pEmmaxHelper::ZEPS));
        pval = pEmmaxHelper::tcdf(t, c-1);
        */

        //if ( ( fabs(sy) > 1e-6 ) || ( fabs(sx) > 1e-6 ) ) {
        //error("sy = %lg, sx = %lg\n",sy,sx);
        //}

        //sdx = sqrt(xsq/c-xsum*xsum/c/(c+1.));
        //sdy = sqrt(ysq/c-ysum*ysum/c/(c+1.));

        //r = (xy/c - xsum*ysum/c/(c+1))/(sdx*sdy);
        //beta = r * sdy / sdx;

        //sebeta = sdy * sqrt( (1-r*r+pEmmaxHelper::ZEPS)*c/(c-1) );

        // get the count of possible genotypes
        emx.tvcf.GENOCNT(i,genocnts);
        if ( binaryFlag ) { emx.tvcf.CASECTRLCNT(i,casectrlcnts,isCases); }

        // CHROM BEG END MARKER_ID NS AC CALLRATE GENOCNT MAF STAT PVALUE BETA SEBETA R2
        //fprintf(stderr,"bar %s:%d %s %s\n",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.markers[i].c_str(),emx.tvcf.refs[i].c_str());
        wf.printf("%s\t%d\t%d\t%s\t%d\t%.2lf\t%.5lf\t%d/%d/%d\t%.5lf\t%.4lf\t%.4lg\t%.4lg\t%.4lg\t%.4lg",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(), emx.tvcf.numAlleles[i]/2, emx.tvcf.sumAlleles[i], emx.tvcf.callRate(i), genocnts[0], genocnts[1], genocnts[2], emx.tvcf.MAF(i), t, pval, beta, sebeta, r*r);
        if ( binaryFlag ) {
          wf.printf("\t%d/%d/%d\t%d/%d/%d",casectrlcnts[0],casectrlcnts[1],casectrlcnts[2],casectrlcnts[3],casectrlcnts[4],casectrlcnts[5]);
        }
        wf.printf("\n");
        ++m;
      }
      else {
        mu = emx.tvcf.alleleFreq(i)*2.;
        emx.tvcf.GENOCNT(i,genocnts);
        if ( binaryFlag ) { emx.tvcf.CASECTRLCNT(i,casectrlcnts,isCases); }
        wf.printf("%s\t%d\t%d\t%s\t%d\t%.2lf\t%.5lf\t%d/%d/%d\t%.5lf\tNA\tNA\tNA\tNA\tNA",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1, emx.tvcf.markers[i].c_str(), emx.tvcf.numAlleles[i]/2, emx.tvcf.sumAlleles[i], emx.tvcf.callRate(i), genocnts[0], genocnts[1], genocnts[2], emx.tvcf.MAF(i));
        if ( binaryFlag ) {
          wf.printf("\t%d/%d/%d\t%d/%d/%d",casectrlcnts[0],casectrlcnts[1],casectrlcnts[2],casectrlcnts[3],casectrlcnts[4],casectrlcnts[5]);
        }
        wf.printf("\n");
      }
    }
    fprintf(stderr,"%d markers included\n",m);

  }
  return 0;
}


// y(t) ~ Xt + Zt + et
int runBurdenAssoc(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("groupf",&arg.groupf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
    LONG_INTPARAMETER("minMAC",&arg.minMAC)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)
    LONG_PARAMETER("sepchr",&arg.sepchr)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("eigf",&arg.ineigf)
    LONG_STRINGPARAMETER("remlf",&arg.remlf)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out-assocf",&arg.assocf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.vcf.empty() || arg.phenof.empty() || arg.groupf.empty() || arg.ineigf.empty() || arg.remlf.empty() ) {
    error("--vcf, --phenof, --groupf, --out-remlf and --kinf are required parameters");
  }

  pEmmax emx;
  emx.loadFiles( arg.phenof.c_str(), arg.covf.c_str(), arg.indf.c_str(), NULL, arg.ineigf.c_str(), arg.vcf.c_str(), NULL, NULL, arg.field.c_str(), !arg.ignoreFilter );

  pFile treml(arg.remlf.c_str());
  double delta = 0;
  const char* line = NULL;
  std::vector<std::string> tokens;
  while( (line = treml.getLine()) != NULL ) {
    pFile::tokenizeLine(line," \t\r\n", tokens);
    if ( tokens.size() != 2 )
      error("Phenotype file must have two columns");
    delta = atof(tokens[1].c_str());
    notice("Reading delta = %lf",delta);
    break;
  }

  int c = emx.evalR.rows();
  int n = (int)emx.inds.size();
  emx.T.resize(c,n);
  for(int i=0; i < c; ++i) {
    for(int j=0; j < n; ++j) {
      emx.T(i,j) = emx.evecR(j,i)/sqrt(emx.evalR(i)+delta);
    }
  }

  VectorXd yt = emx.T * emx.y; // (n-p) * 1 matrix
  VectorXd x, xt;
  x.resize(n);

  double sy = 0, syy = 0;
  for(int i=0; i < c; ++i) {
    sy += yt(i);
    syy += yt(i)*yt(i);
  }

  wFile wf(arg.assocf.c_str());
  double mu = 0;
  int i,m,j,nPass;

  wf.printf("#CHROM\tBEG\tEND\tMARKER_ID\tTOT_MARKERS\tPASS_MARKERS\tBURDEN_CNT\tFRAC_WITH_RARE\tSTAT\tPVALUE\tBETA\tSEBETA\tR2\n");

  // here comes genotype reading part
  // start reading groupf, line by, line
  pFile tgrp(arg.groupf.c_str());
  while( (line = tgrp.getLine()) != NULL ) {
    std::vector<std::string> tokens;
    //double xsum = 0, nx = 0, xsq = 0, xy = 0, r = 0, t = 0, beta = 0, sebeta = 0, sdx = 0, sdy = 0, pval = 0;
    double sx = 0, sxx = 0, sxy = 0, r = 0, t = 0, beta = 0, sebeta = 0, pval = 0, varE = 0, nx = 0;

    if ( line[0] == '#' ) continue;
    //notice("foo");
    //notice("Parsing line %s",line);
    pFile::tokenizeLine(line," \t\r\n", tokens);
    //notice("goo");

    //notice("tokens.size() = %d",(int)tokens.size());

    // tokens[0] is group name and others are markers
    emx.tvcf.clear();
    m = emx.tvcf.readMarkerGroup(tokens,1,arg.sepchr);
    //notice("Processing %d markers in group %s",m,tokens[0].c_str());

    // create burden test collapsing variables
    // initial version ignores missing genotypes
    nPass = 0;
    for(i=0; i < n; ++i) { x(i) = 0; }

    for(i=0; i < m; ++i) {
      if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) &&
           ( emx.tvcf.MAF(i) <= arg.maxMAF ) &&
           ( emx.tvcf.MAC(i) >= arg.minMAC ) &&
           ( emx.tvcf.callRate(i) >= arg.minCallRate ) &&
           ( ( arg.minRSQ == 0 ) || ( emx.tvcf.RSQ(i) >= arg.minRSQ ) ) ) {

        int offset = i * emx.tvcf.nInds;
        ++nPass;
        mu = emx.tvcf.alleleFreq(i)*2.;
        for(j=0; j < emx.tvcf.nInds; ++j) {
          if ( !std::isnan(emx.tvcf.genos[j + offset]) ) {
            if ( mu < 1 ) { // AF < .5
              if ( emx.tvcf.genos[j + offset] > DS_THRES ) x(j) = 1.;
            }
            else {
              if ( emx.tvcf.genos[j + offset] < 2-DS_THRES ) x(j) = 1.;
            }
          }
        }
      }
    }

    for(i=0; i < n; ++i) {
      //notice("%d %lf",j,x(j));
      nx += x(i);
    }
    //notice("nx=%lf",nx);
    //abort();
    if ( ( nPass > 0 ) && ( nx > 0 ) && ( nx < n ) ) {
      xt = emx.T * x;
      //xy = xsum = xsq = 0;
      sx = sxx = sxy = 0;
      for(j=0; j < c; ++j) {
        //xy += xt(j)*yt(j);
        //xsum += xt(j);
        //xsq += xt(j)*xt(j);
        sxy += (xt(j)*yt(j));
        sx += xt(j);
        sxx += (xt(j)*xt(j));
      }

      beta = ((c+1.)*sxy-sx*sy)/((c+1.)*sxx-sx*sx);
      varE = 1/(c+1.)/(c-1.)*((c+1.)*syy-sy*sy-beta*beta*((c+1.)*sxx-sx*sx));
      //notice("varE=%lf",varE);
      sebeta = sqrt((c+1.)*varE/((c+1.)*sxx-sx*sx));
      r = ((c+1)*sxy-sx*sy)/sqrt(((c+1)*sxx-sx*sx)*((c+1)*syy-sy*sy));
      t = r * sqrt((c-1)/(1-r*r+pEmmaxHelper::ZEPS));
      pval = pEmmaxHelper::tcdf(t, c-1);

      //sdx = sqrt(xsq/c-xsum*xsum/c/c);
      //sdy = sqrt(ysq/c-ysum*ysum/c/c);

      //r = (xy/c - xsum*ysum/c/c)/(sdx*sdy);
      //beta = r * sdy / sdx;
      //sebeta = sdy * sqrt( (1-r*r+pEmmaxHelper::ZEPS)*c/(c-1) );
      //t = r * sqrt((c-1)/(1-r*r+pEmmaxHelper::ZEPS));
      //pval = pEmmaxHelper::tcdf(t, c-1);

      wf.printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.5lf\t%.4lf\t%.4lg\t%.4lg\t%.4lg\t%.4lg\n",emx.tvcf.groupChrom.c_str(),emx.tvcf.groupBeg,emx.tvcf.groupEnd,emx.tvcf.groupID.c_str(), m, nPass, (int)nx, nx < n ? (double)nx/n : (double)(n-nx)/n,t,pval,beta,sebeta,r*r);
      //wf.printf("%s\t%d\t%d\t%d\t%.5lf\t%.4lf\t%.4lg\t%.4lf\t%.4lf\t%.5lf\n",tvcf.markerGroupName.c_str(), m, nPass, (int)nx, nx < n ? (double)nx/n : (double)(n-nx)/n,t,pval,beta,sebeta,r*r);
    }
    else {
      wf.printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n",emx.tvcf.groupChrom.c_str(),emx.tvcf.groupBeg,emx.tvcf.groupEnd,emx.tvcf.groupID.c_str(), m, nPass, (int)nx);
      //wf.printf("%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n",tvcf.markerGroupName.c_str(), m, nPass, (int)nx);
    }
    //notice("bar");
  }
  return 0;

}

// y(t) ~ Xt + Zt + et
int runEmmaxVT(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("groupf",&arg.groupf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
    LONG_INTPARAMETER("minMAC",&arg.minMAC)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)
    LONG_PARAMETER("sepchr",&arg.sepchr)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("eigf",&arg.ineigf)
    LONG_STRINGPARAMETER("remlf",&arg.remlf)
    LONG_STRINGPARAMETER("scoref",&arg.scoref)
    LONG_INTPARAMETER("seed",&arg.seed)
    LONG_INTPARAMETER("maxperm",&arg.maxperm)
    LONG_INTPARAMETER("minperm",&arg.minperm)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out-assocf",&arg.assocf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.vcf.empty() || arg.phenof.empty() || arg.groupf.empty() || arg.ineigf.empty() || arg.remlf.empty() ) {
    error("--vcf, --phenof, --groupf, --out-remlf and --kinf are required parameters");
  }

  pEmmax emx;
  emx.loadFiles( arg.phenof.c_str(), arg.covf.c_str(), arg.indf.c_str(), NULL, arg.ineigf.c_str(), arg.vcf.c_str(), NULL, NULL, arg.field.c_str(), !arg.ignoreFilter );

  pFile treml(arg.remlf.c_str());
  double delta = 0;
  const char* line = NULL;
  std::vector<std::string> tokens;
  while( (line = treml.getLine()) != NULL ) {
    pFile::tokenizeLine(line," \t\r\n", tokens);
    if ( tokens.size() != 2 )
      error("Phenotype file must have two columns");
    delta = atof(tokens[1].c_str());
    notice("Reading delta = %lf",delta);
    break;
  }

  int c = emx.evalR.rows();
  int n = (int)emx.inds.size();
  emx.T.resize(c,n);
  for(int i=0; i < c; ++i) {
    for(int j=0; j < n; ++j) {
      emx.T(i,j) = emx.evecR(j,i)/sqrt(emx.evalR(i)+delta);
    }
  }

  // transform the coordinates
  VectorXd yt = emx.T * emx.y; // (n-p) * 1 matrix
  double sumyt = yt.sum();
  double sqyt = yt.squaredNorm();
  double nvaryt = sqyt - sumyt*sumyt/c;
  MatrixXd Ytr(arg.minperm+1, c);
  MatrixXd Ytp(arg.minperm,c);
  //MatrixXd Yr(arg.minperm+1, n);
  //MatrixXd Yp(arg.minperm,n);

  std::vector<int> irand;
  int i, j, k;
  /*
  for(i=0; i < n; ++i) {
    irand.push_back(i);
    Yr(0,i) = y(i);
  }
  for(i=1; i < arg.minperm+1; ++i) {
    std::random_shuffle(irand.begin(), irand.end());
    for(j=0; j < n; ++j) {
      Yr(i,j) = y(irand[j]);
    }
  }
  Ytr = Yr * emx.T.transpose();

  VectorXd sumYtr = Ytr.rowwise.sum();
  VectorXd sqYtr = Ytr.rowwise.squaredNorm();
  VectorXd nVarYtr = sqYtr - sumYtr.array().square()/c;
  */

  for(i=0; i < c; ++i) {
    irand.push_back(i);
    Ytr(0,i) = yt(i);
  }
  for(i=1; i < arg.minperm+1; ++i) {
    std::random_shuffle(irand.begin(), irand.end());
    for(j=0; j < c; ++j) {
      Ytr(i,j) = yt(irand[j]);
    }
  }

  // create matrix of collapsing variables
  wFile wf(arg.assocf.c_str());
  double mu = 0;
  int m;

  wf.printf("#CHROM\tBEG\tEND\tMARKER_ID\tTOT_MARKERS\tPASS_MARKERS\tBURDEN_CNT\tFRAC_WITH_RARE\tSTAT\tPVALUE\tR2\tDIRECTION\tOPT_THRES_RAC_OR_SCORE\tOPT_FRAC_WITH_RARE\n");

  // here comes genotype reading part
  // start reading groupf, line by, line
  pFile tgrp(arg.groupf.c_str());
  genomeScore gScore;
  if ( !arg.scoref.empty() ) {
    gScore.setDir(arg.scoref.c_str());
  }

  while( (line = tgrp.getLine()) != NULL ) {
    // Create collapsing variables

    std::vector<std::string> tokens;
    if ( line[0] == '#' ) continue;
    pFile::tokenizeLine(line," \t\r\n", tokens);

    //notice("Reading group %s",tokens[0].c_str());

    emx.tvcf.clear();
    m = emx.tvcf.readMarkerGroup(tokens,1,arg.sepchr);

    std::vector<int> g;
    std::vector<int> macs;

    for(i=0; i < m; ++i) {
      if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) &&
           ( emx.tvcf.MAF(i) <= arg.maxMAF ) &&
           ( emx.tvcf.MAC(i) >= arg.minMAC ) &&
           ( emx.tvcf.callRate(i) >= arg.minCallRate ) &&
           ( ( arg.minRSQ == 0 ) || ( emx.tvcf.RSQ(i) >= arg.minRSQ ) ) ) {

        int offset = i * emx.tvcf.nInds;
        int mac = 0;
        mu = emx.tvcf.alleleFreq(i)*2.;
        for(j=0; j < emx.tvcf.nInds; ++j) {
          if ( !std::isnan(emx.tvcf.genos[j + offset]) ) {
            if ( mu < 1 ) { // AF < .5
              g.push_back(emx.tvcf.genos[j + offset] > DS_THRES ? 1 : 0 );
            }
            else {
              g.push_back(emx.tvcf.genos[j + offset] < 2-DS_THRES ? 1 : 0 );
            }
            mac += g.back();
          }
          else {
            g.push_back(0);
          }
        }
        if ( arg.scoref.empty() ) {
          macs.push_back(mac);
        }
        else {
          // use GERP/PhyloP score * 1000 as threshold
          macs.push_back((int)floor(gScore.baseScore(emx.tvcf.markers[i].c_str())*(-1000)+.5));
        }
      }
    }

    int mv = (int)macs.size(); // number of valid markers

    if ( mv > 0 ) { // at least one marker
      std::set<int> umacSet;
      for(i=0; i < mv; ++i) {
        umacSet.insert(macs[i]);
      }

      std::vector<int> umacs;
      std::map<int,int> iumacs;
      for(std::set<int>::iterator it = umacSet.begin();
          it != umacSet.end(); ++it) {
        iumacs[*it] = umacs.size();
        umacs.push_back(*it);
      }
      int numacs = (int)umacs.size(); // number of unique MACs

      //notice("n = %d, mv = %d, numacs = %d",n, mv, numacs);
      //notice("g.size() = %d",(int)g.size());

      // generate ind * mac incidence matrix
      MatrixXd C = MatrixXd::Zero(n,numacs); // matrix

      // this is 'collapsing-style' variable threshold test
      for(i=0, k=0; i < mv; ++i) {
        for(j=0; j < n; ++j, ++k) {
          if ( g[k] > 0 ) C(j,iumacs[macs[i]]) = 1.;
        }
      }

      // cumulate the incidence matrix
      for(i=1; i < numacs; ++i) {
        for(j=0; j < n; ++j) {
          if ( C(j,i-1) > 0 ) C(j,i) = 1.;
        }
      }

      //notice("Transforming collpasing vector (%d x %d) x (%d x %d) x (%d x %d)",Ytr.rows(), Ytr.cols(), emx.T.rows(), emx.T.cols(), C.rows(), C.cols());

      // C contains collapsing variables
      // Now transform matrix C
      MatrixXd Ct = emx.T * C; // Ct is (n-p) * m matrix
      VectorXd sumCt = Ct.colwise().sum();
      VectorXd sqCt = Ct.colwise().squaredNorm();
      VectorXd sumC = C.colwise().sum();

      // calculate correlation coefficient
      MatrixXd R = Ytr * Ct; // (nperm+1) * numacs matrix
      int optMAC = 0;
      double optBurden = 0;
      double optR2 = 0;
      int optSgn = 0;
      int sgn;

      //notice("Calculating VT statistics");

      for(i=0; i < numacs; ++i) {
        for(j=0; j < arg.minperm+1; ++j) {
          sgn = R(j,i) > 0 ? 1 : -1;
          R(j,i) -= (sumyt * sumCt(i) / c);
          R(j,i) = (R(j,i) * R(j,i)) / ((sqCt(i)-(double)sumCt(i)*sumCt(i)/c)*nvaryt+pEmmaxHelper::ZEPS);
          //R(j,i) = (R(j,i) * R(j,i) - sumyt * sumCt(i)) / ((sqCt(i)-(double)sumCt(i)*sumCt(i)/c)*nvaryt+pEmmaxHelper::ZEPS);
          if ( j == 0 ) {
            if ( R(j,i) > optR2 ) {
              optR2 = R(j,i);
              optBurden = sumC(i);
              optMAC = umacs[i];
              optSgn = sgn;
            }
          }
        }
      }

      // calculate max chi-squared statistics
      VectorXd maxR2 = R.rowwise().maxCoeff();

      // calculate the rank of p-value
      int rank = 0;
      for(i=1; i < arg.minperm+1; ++i) {
        if ( optR2 <= maxR2(i) ) ++rank;
      }

      int burdenMAC = (int)sumC(numacs-1);
      int perm = arg.minperm;

      if ( rank * 10 < arg.minperm ) {
        for(perm=arg.minperm*10; perm <= arg.maxperm; perm *= 10) {
          int rperm = perm*9/10;
          notice("Performing adaptive permutation for additional %d permutations for group %s -- maxR2 = %lg, rank = %d",rperm,tokens[0].c_str(),optR2,rank);
          for(int p=0; p < rperm; p += arg.minperm) {
            for(i=0; i < arg.minperm; ++i) {
              std::random_shuffle(irand.begin(), irand.end());
              for(j=0; j < c; ++j) {
                Ytp(i,j) = yt(irand[j]);
              }
            }
            //Ytp = Yp * emx.T.transpose();
            //VectorXd sumYtp = Ytp.rowwise.sum();
            //VectorXd sqYtp = Ytp.rowwise.squaredNorm();
            //VectorXd nVarYtp = sqYtp - sumYtp.array().square()/c;

            R = Ytp * Ct; // (nperm+1) * numacs matrix
            for(i=0; i < numacs; ++i) {
              for(j=0; j < arg.minperm; ++j) {
                R(j,i) -= (sumyt * sumCt(i)/c);
                R(j,i) = (R(j,i) * R(j,i)) / ((sqCt(i)-(double)sumCt(i)*sumCt(i)/c)*nvaryt+pEmmaxHelper::ZEPS);
                //R(j,i) = (R(j,i) * R(j,i) - sumyt * sumCt(i)) / ((sqCt(i)-(double)sumCt(i)*sumCt(i)/c)*nvaryt+pEmmaxHelper::ZEPS);
              }
            }
            maxR2 = R.rowwise().maxCoeff();
            for(i=0; i < arg.minperm; ++i) {
              if ( optR2 <= maxR2(i) ) ++rank;
            }
          }
          if ( rank * 10 > arg.minperm ) break;
        }
        if ( perm > arg.maxperm ) perm = arg.maxperm;
      }
      wf.printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.5lf\t%.4lf\t%.4lg\t%.4lg\t%s\t%d\t%.5lf\n",emx.tvcf.groupChrom.c_str(),emx.tvcf.groupBeg,emx.tvcf.groupEnd,emx.tvcf.groupID.c_str(), m, mv, burdenMAC, (double)burdenMAC/n, optR2*(c-1), (rank+1>perm ? 1. : (double)(rank+1.)/perm),optR2,optSgn > 0 ? "+" : "-",optMAC,(double)optBurden/n);
    }
    else {
      wf.printf("%s\t%d\t%d\t%s\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n",emx.tvcf.groupChrom.c_str(),emx.tvcf.groupBeg,emx.tvcf.groupEnd,emx.tvcf.groupID.c_str(), m, mv);
    }
  }
  return 0;
}

double chisq1cdf(double x) {
  if ( x < 0 ) x = 0;
  return erfc(sqrt(x/2.));
}

// y ~ X + Z + e
// y(t) ~ Xt + Zt + et
int runGLRT(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  arg.field = "PL";
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_INTPARAMETER("unit",&arg.unit)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
    LONG_INTPARAMETER("minMAC",&arg.minMAC)
    LONG_INTPARAMETER("maxMAC",&arg.maxMAC)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("indf",&arg.indf)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out",&arg.outf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.phenof.empty() || arg.vcf.empty() || arg.outf.empty() ) {
    error("--phenof, --vcf -and --outf are required parameters");
  }

  pEmmax emx;
  emx.loadFiles(arg.phenof.c_str(), arg.covf.c_str(), arg.indf.c_str(), NULL, NULL, arg.vcf.c_str(), arg.region.c_str(), arg.rule.c_str(), arg.field.c_str(), !arg.ignoreFilter);

  wFile wf(arg.outf.c_str());
  int i,m;

  // identify the indices of individuals based on VCF
  std::vector<int> g1, g2;
  for(i=0; i < emx.tvcf.nInds; ++i) {
    if ( emx.y(i) == 1 ) {
      g1.push_back(i);
    }
    else if ( emx.y(i) == 2 ) {
      g2.push_back(i);
    }
    else {
      error("Cannot recognize phenotype value %s. Only 1/2 encoded binary phenotypes are acceptable",emx.tvcf.inds[i].c_str());
    }
  }

  wf.printf("#CHROM\tBEG\tEND\tMARKER_ID\tNS\tMAF\tLLR\tPVALUE\tAF\tAF1\tAF2\n");
  double af, maf, af1, af2;
  for(int M=0; emx.tvcf.readMarkers(arg.unit); ) {
    M += emx.tvcf.nMarkers;
    fprintf(stderr,"Processing %d markers across %d individuals...",M, emx.tvcf.nInds);
    for(i=0, m=0; i < emx.tvcf.nMarkers; ++i) {
      af = emx.tvcf.alleleFreq(i);
      maf = af > 0.5 ? 1-af : af;
      if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) &&
           ( emx.tvcf.MAF(i) <= arg.maxMAF ) &&
           ( emx.tvcf.MAC(i) <= arg.maxMAC ) &&
           ( emx.tvcf.MAC(i) >= arg.minMAC ) &&
           ( emx.tvcf.callRate(i) >= arg.minCallRate ) ) {
        double lrt = emx.tvcf.LRT(i,g1,g2,af1,af2);
        double pval = chisq1cdf(lrt);

        // MARKER_ID NS MAF STAT PVALUE AF AF1 AF2
        wf.printf("%s\t%d\t%d\t%s\t%d\t%.5lf\t%.4lf\t%.4lg\t%.5lf\t%.5lf\t%.5lf\n",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(),emx.tvcf.numAlleles[i]/2,maf,lrt,pval,af,af1,af2);
        ++m;
      }
      else {
        wf.printf("%s\t%d\t%d\t%s\t%d\t%.5lf\tNA\tNA\t%.5lf\tNA\tNA\n",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(),emx.tvcf.numAlleles[i]/2,maf,af);
      }
    }
    fprintf(stderr,"%d markers included\n",m);
  }
  return 0;
}

// y ~ X + Z + e
// y(t) ~ Xt + Zt + et
int runSimul(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Options")
    LONG_STRINGPARAMETER("kinf",&arg.kinf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_PARAMETER("binary",&arg.binary)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_INTPARAMETER("nsimul",&arg.nsimul)
    LONG_DOUBLEPARAMETER("h2",&arg.h2)
    LONG_STRINGPARAMETER("out-phenof",&arg.phenof)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.phenof.empty() || arg.kinf.empty() ) {
    error("--out-phenof, --kinf or --vcf are required parameters");
  }

  pEmmax emx;
  if ( arg.vcf.empty() ) {
    emx.loadFiles(NULL, NULL, arg.indf.c_str(), arg.kinf.c_str(), NULL, NULL, NULL, NULL, NULL, true);
  }

  // set a random seed based on timestamp
  srand(std::time(0));

  notice("Calculating Cholesky Decomposition of the Kinship matrix");
  int n = (int)emx.inds.size();
  //LDLT<MatrixXd> lt = Ks.ldlt();
  LLT<MatrixXd> lt = emx.Ks.llt();
  MatrixXd Y(n,arg.nsimul);
  //MatrixXd L = lt.matrixL() * lt.vectorD().array().sqrt().asDiagonal();

  for(int j=0; j < arg.nsimul; ++j) {
    // generate random numbers following normal distribution
    // using Box-Muller Transformation
    VectorXd r1(n), r2(n);
    std::vector<double> v;
    for(int i=0; i < n; ++i) {
      double x1 = (rand()+.5)/(RAND_MAX+1.);
      double x2 = (rand()+.5)/(RAND_MAX+1.);
      r1[i] = sqrt(-2.0*log(x1)) * cos(M_PI * 2 * x2);
      r2[i] = sqrt(-2.0*log(x1)) * sin(M_PI * 2 * x2);
    }
    //Y.col(j) = sqrt(arg.h2) * (lt.matrixL() * r1) + sqrt(1.-arg.h2) * r2;
    Y.col(j) = (lt.matrixL() * r1) * sqrt(arg.h2) + sqrt(1.-arg.h2) * r2;
    if ( arg.binary  ) {
      for(int i=0; i < n; ++i) {
	v.push_back(Y(i,j));
      }
      std::sort(v.begin(),v.end());
      // currently supports 50:50
      double median = v[n/2];
      for(int i=0; i < n; ++i) {
	if ( Y(i,j) <= median ) {
	  Y(i,j) = 1;
	}
	else {
	  Y(i,j) = 2;
	}
      }
    }
  }

  notice("Writing simulated phenotypes in PED format");
  // write output files in PED format
  wFile wf(arg.phenof.c_str());
  wf.printf("#FAM_ID\tIND_ID\tFAT_ID\tMOT_ID\tSEX");
  for(int j=0; j < arg.nsimul; ++j) {
    wf.printf("\tSIMUL%d",j+1);
  }
  wf.printf("\n");
  for(int i=0; i < n; ++i) {
    wf.printf("%s\t%s\t0\t0\t0",emx.inds[i].c_str(), emx.inds[i].c_str());
    for(int j=0; j < arg.nsimul; ++j) {
      wf.printf("\t%lf",Y(i,j));
    }
    wf.printf("\n");
  }
  return 0;
}

/*
int runGLST(int argc, char** argv) {
  int orgArgc = argc;
  char** orgArgv = argv;

  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_INTPARAMETER("unit",&arg.unit)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("covf",&arg.covf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_PARAMETER("no-intercept",&arg.noIntercept)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out",&arg.outf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.phenof.empty() || arg.vcf.empty() || arg.outf.empty()  ) {
    error("--phenof, --vcf, --out are required parameters");
  }

  pEmmax emx;
  emx.loadFiles(arg.phenof.c_str(), arg.covf.c_str(), arg.indf.c_str(), NULL, NULL, arg.vcf.c_str(), arg.region.c_str(), arg.rule.c_str(), arg.field.c_str(), !arg.ignoreFilter);

  frequency_estimator freqest(&emx.X);

  wFile wf(arg.outf.c_str());
  int i, m;

  for(int M=0; emx.tvcf.readMarkers(arg.unit); ) {
    M += emx.tvcf.nMarkers;

    fprintf(stderr,"Processing %d markers across %d individuals...",M, emx.tvcf.nInds);
    for(i=0, m=0; i < emx.tvcf.nMarkers; ++i) {
      // 
      
      af = emx.tvcf.alleleFreq(i);
      maf = af > 0.5 ? 1-af : af;
      if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) && 
	   ( emx.tvcf.MAF(i) <= arg.maxMAF ) && 
	   ( emx.tvcf.MAC(i) <= arg.maxMAC ) && 
	   ( emx.tvcf.MAC(i) >= arg.minMAC ) && 
	   ( emx.tvcf.callRate(i) >= arg.minCallRate ) ) {
	double lrt = emx.tvcf.LRT(i,g1,g2,af1,af2);
	double pval = chisq1cdf(lrt);

	// MARKER_ID NS MAF STAT PVALUE AF AF1 AF2
	wf.printf("%s\t%d\t%d\t%s\t%d\t%.5lf\t%.4lf\t%.4lg\t%.5lf\t%.5lf\t%.5lf\n",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(),emx.tvcf.numAlleles[i]/2,maf,lrt,pval,af,af1,af2);
	++m;
      }
      else {
	wf.printf("%s\t%d\t%d\t%s\t%d\t%.5lf\tNA\tNA\t%.5lf\tNA\tNA\n",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(),emx.tvcf.numAlleles[i]/2,maf,af);
      }
    }
    fprintf(stderr,"%d markers included\n",m);
    

  // fit the null model first

  std::vector<int> g1, g2;

  // Add additional parmameters : 
  // --out-kinf, --out-eigR, --out-remlf, --out-assocf
  std::string outKinf = arg.outf + ".kin";
  std::string outEigf = arg.outf + ".eigR";
  std::string outRemlf = arg.outf + ".reml";
  std::string outAssocf = arg.outf + ".assoc";

  char** newArgv = new char* [orgArgc + 100];
  for(int i=0; i < orgArgc; ++i) newArgv[i] = orgArgv[i];

  notice("---------------------------------------------------------");
  notice("Running %s gen-kin command",orgArgv[0]);
  notice("---------------------------------------------------------");
  newArgv[orgArgc] = (char*)"--out-kinf"; 
  newArgv[orgArgc+1] = (char*)outKinf.c_str();
  runGenKin(orgArgc+1,newArgv+1);

  notice("---------------------------------------------------------");
  notice("Running %s reml command",orgArgv[0]);
  notice("---------------------------------------------------------");
  newArgv[orgArgc] = (char*)"--kinf"; 
  newArgv[orgArgc+1] = (char*)outKinf.c_str();
  newArgv[orgArgc+2] = (char*)"--out-eigf"; 
  newArgv[orgArgc+3] = (char*)outEigf.c_str();
  newArgv[orgArgc+4] = (char*)"--out-remlf"; 
  newArgv[orgArgc+5] = (char*)outRemlf.c_str();
  runReml(orgArgc+5,newArgv+1);

  notice("---------------------------------------------------------");
  notice("Running %s assoc command",orgArgv[0]);
  notice("---------------------------------------------------------");
  newArgv[orgArgc] = (char*)"--eigf"; 
  newArgv[orgArgc+1] = (char*)outEigf.c_str();
  newArgv[orgArgc+2] = (char*)"--remlf"; 
  newArgv[orgArgc+3] = (char*)outRemlf.c_str();
  newArgv[orgArgc+4] = (char*)"--out-assocf"; 
  newArgv[orgArgc+5] = (char*)outAssocf.c_str();
  runAssoc(orgArgc+5,newArgv+1);

  delete [] newArgv;
  return 0;
}
*/

int runAll(int argc, char** argv) {
  int orgArgc = argc;
  char** orgArgv = argv;

  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_INTPARAMETER("unit",&arg.unit)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("covf",&arg.covf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_PARAMETER("no-intercept",&arg.noIntercept)

    LONG_PARAMETER_GROUP("Empircial Kinship Method")
    EXCLUSIVE_PARAMETER("BN",&arg.BN)
    EXCLUSIVE_PARAMETER("IBS",&arg.IBS)
    EXCLUSIVE_PARAMETER("YV",&arg.YV)

    LONG_PARAMETER_GROUP("Kinship Normalization")
    LONG_PARAMETER("raw",&arg.raw)
    LONG_PARAMETER("cov",&arg.cov)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out",&arg.outf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.phenof.empty() || arg.vcf.empty() || arg.outf.empty()  ) {
    error("--phenof, --vcf, --out are required parameters");
  }

  // Add additional parmameters : 
  // --out-kinf, --out-eigR, --out-remlf, --out-assocf
  std::string outKinf = arg.outf + ".kin";
  std::string outEigf = arg.outf + ".eigR";
  std::string outRemlf = arg.outf + ".reml";
  std::string outAssocf = arg.outf + ".assoc";

  char** newArgv = new char* [orgArgc + 100];
  for(int i=0; i < orgArgc; ++i) newArgv[i] = orgArgv[i];

  notice("---------------------------------------------------------");
  notice("Running %s gen-kin command",orgArgv[0]);
  notice("---------------------------------------------------------");
  newArgv[orgArgc] = (char*)"--out-kinf"; 
  newArgv[orgArgc+1] = (char*)outKinf.c_str();
  runGenKin(orgArgc+1,newArgv+1);

  notice("---------------------------------------------------------");
  notice("Running %s reml command",orgArgv[0]);
  notice("---------------------------------------------------------");
  newArgv[orgArgc] = (char*)"--kinf"; 
  newArgv[orgArgc+1] = (char*)outKinf.c_str();
  newArgv[orgArgc+2] = (char*)"--out-eigf"; 
  newArgv[orgArgc+3] = (char*)outEigf.c_str();
  newArgv[orgArgc+4] = (char*)"--out-remlf"; 
  newArgv[orgArgc+5] = (char*)outRemlf.c_str();
  runReml(orgArgc+5,newArgv+1);

  notice("---------------------------------------------------------");
  notice("Running %s assoc command",orgArgv[0]);
  notice("---------------------------------------------------------");
  newArgv[orgArgc] = (char*)"--eigf"; 
  newArgv[orgArgc+1] = (char*)outEigf.c_str();
  newArgv[orgArgc+2] = (char*)"--remlf"; 
  newArgv[orgArgc+3] = (char*)outRemlf.c_str();
  newArgv[orgArgc+4] = (char*)"--out-assocf"; 
  newArgv[orgArgc+5] = (char*)outAssocf.c_str();
  runAssoc(orgArgc+5,newArgv+1);

  delete [] newArgv;
  return 0;
}

int runRemlAssoc(int argc, char** argv) {
  int orgArgc = argc;
  char** orgArgv = argv;

  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_INTPARAMETER("unit",&arg.unit)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("covf",&arg.covf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("kinf",&arg.kinf)
    LONG_PARAMETER("no-intercept",&arg.noIntercept)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out",&arg.outf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.phenof.empty() || arg.vcf.empty() || arg.outf.empty()  ) {
    error("--phenof, --vcf, --kinf, --out are required parameters (--covf and --indf are also recommended)");
  }

  // Add additional parmameters : 
  // --out-kinf, --out-eigR, --out-remlf, --out-assocf
  std::string outEigf = arg.outf + ".eigR";
  std::string outRemlf = arg.outf + ".reml";
  std::string outAssocf = arg.outf + ".assoc";

  char** newArgv = new char* [orgArgc + 100];
  for(int i=0; i < orgArgc; ++i) newArgv[i] = orgArgv[i];

  notice("---------------------------------------------------------");
  notice("Running %s reml command",orgArgv[0]);
  notice("---------------------------------------------------------");
  newArgv[orgArgc+0] = (char*)"--out-eigf"; 
  newArgv[orgArgc+1] = (char*)outEigf.c_str();
  newArgv[orgArgc+2] = (char*)"--out-remlf"; 
  newArgv[orgArgc+3] = (char*)outRemlf.c_str();
  runReml(orgArgc+3,newArgv+1);

  notice("---------------------------------------------------------");
  notice("Running %s assoc command",orgArgv[0]);
  notice("---------------------------------------------------------");
  newArgv[orgArgc] = (char*)"--eigf"; 
  newArgv[orgArgc+1] = (char*)outEigf.c_str();
  newArgv[orgArgc+2] = (char*)"--remlf"; 
  newArgv[orgArgc+3] = (char*)outRemlf.c_str();
  newArgv[orgArgc+4] = (char*)"--out-assocf"; 
  newArgv[orgArgc+5] = (char*)outAssocf.c_str();
  runAssoc(orgArgc+5,newArgv+1);

  delete [] newArgv;
  return 0;
}

// variable threshold test
int runVT(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  std::string vntf;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("groupf",&arg.groupf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
    LONG_INTPARAMETER("minMAC",&arg.minMAC)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)
    LONG_PARAMETER("sepchr",&arg.sepchr)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("covf",&arg.covf)
    LONG_STRINGPARAMETER("scoref",&arg.scoref)
    LONG_INTPARAMETER("seed",&arg.seed)
    LONG_INTPARAMETER("maxperm",&arg.maxperm)
    LONG_INTPARAMETER("minperm",&arg.minperm)
    LONG_PARAMETER("recessive",&arg.recessive)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out-assocf",&arg.assocf)
    LONG_STRINGPARAMETER("out-variantf",&vntf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.vcf.empty() || arg.phenof.empty() || arg.groupf.empty() || arg.assocf.empty() || vntf.empty() ) {
    error("--vcf, --phenof, --groupf, --out-assocf, --out-variantf are required parameters");
  }

  if ( arg.seed == 0 ) {
    arg.seed = (unsigned int)std::time(NULL);
  }
  srand(arg.seed);

  // p : # permutation
  // m : # of markers (or distinct allele count)
  // maximum number of repetition is p * m
  // first, regression y by covariates

  // ID matching strategy
  // respect orders in the vcf file
  // identify overlapping markers between phe, indf
  // assuming covf will be matching perfectly (can be verified by --in-eig)

  pEmmax emx;
  emx.loadFiles(arg.phenof.c_str(), arg.covf.c_str(), arg.indf.c_str(), NULL, NULL, arg.vcf.c_str(), arg.region.c_str(), NULL, arg.field.c_str(), !arg.ignoreFilter);

  // regress out the covariates
  int n = emx.inds.size();
  JacobiSVD<MatrixXd> svd(emx.X, ComputeThinU | ComputeThinV);
  MatrixXd betasSvd = svd.solve(emx.y);
  VectorXd yr = emx.y - emx.X * betasSvd; // yr is residual phenotypes
  MatrixXd Yr(arg.minperm+1, n); // (nperm+1) * n matrix
  MatrixXd Yp(arg.minperm, n); // (nperm+1) * n matrix

  //std::cout << betasSvd << std::endl;
  //std::cout << X << std::endl;
  //std::cout << yr << std::endl;

  // fill in the minperm matrix
  std::vector<int> irand;
  double sumyr = 0, sqyr = 0;
  int i, j, k;
  for(i=0; i < n; ++i) {
    irand.push_back(i);
    Yr(0,i) = yr(i);
    sumyr += yr(i);
    sqyr += yr(i)*yr(i);
  }
  for(i=1; i < arg.minperm+1; ++i) {
    std::random_shuffle(irand.begin(), irand.end());
    for(j=0; j < n; ++j) {
      Yr(i,j) = yr(irand[j]);
    }
  }
  if ( fabs(sumyr) > 1e-6 ) {
    error("Residual y has cumulative sum %lf, sqsum %lf",sumyr,sqyr);
  }

  wFile wf(arg.assocf.c_str());
  wFile vf(vntf.c_str());
  double mu = 0;
  int m;

  wf.printf("#CHROM\tBEG\tEND\tMARKER_ID\tTOT_MARKERS\tPASS_MARKERS\tBURDEN_CNT\tFRAC_WITH_RARE\tSTAT\tPVALUE\tR2\tDIRECTION\tOPT_THRES_%s\tOPT_FRAC_WITH_RARE\n",arg.scoref.empty() ? "RAC" : "SCORE");
  vf.printf("#GROUP\tMARKER\tAF\tMAC\tOPT_MAC\tIND\tGENOTYPE\n");

  // here comes genotype reading part
  // start reading groupf, line by, line
  pFile tgrp(arg.groupf.c_str());
  genomeScore gScore;
  if ( !arg.scoref.empty() ) {
    gScore.setDir(arg.scoref.c_str());
  }

  const char* line = NULL;
  while( (line = tgrp.getLine()) != NULL ) {
    std::vector<std::string> tokens;
    if ( line[0] == '#' ) continue;
    pFile::tokenizeLine(line," \t\r\n", tokens);

    //notice("Reading group %s",tokens[0].c_str());

    emx.tvcf.clear();
    m = emx.tvcf.readMarkerGroup(tokens,1,arg.sepchr);

    std::vector<int> g;     // genotypes of (n * m) matrix
    std::vector<int> rg;     // genotypes of (n * m) matrix
    std::vector<int> macs;  // MAC or score for m markers
    std::vector<int> ipass;
    const char* genoLabels[4] = {"./.","0/0","0/1","1/1"};

    for(i=0; i < m; ++i) {
      if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) &&
           ( emx.tvcf.MAF(i) <= arg.maxMAF ) &&
           ( emx.tvcf.MAC(i) >= arg.minMAC ) &&
           ( emx.tvcf.callRate(i) >= arg.minCallRate ) &&
           ( ( arg.minRSQ == 0 ) || ( emx.tvcf.RSQ(i) >= arg.minRSQ ) ) &&
           ( (!arg.recessive) || ( emx.tvcf.HOMMINC(i) > 0 ) )
        ) {

        int offset = i * emx.tvcf.nInds;
        int mac = emx.tvcf.MAC(i);
        mu = emx.tvcf.alleleFreq(i)*2.;
        for(j=0; j < emx.tvcf.nInds; ++j) {
          if ( !std::isnan(emx.tvcf.genos[j + offset]) ) {
            rg.push_back(1+(int)emx.tvcf.genos[j + offset]);
            if ( mu < 1 ) { // AF < .5
              if ( arg.recessive ) {
                g.push_back(emx.tvcf.genos[j + offset] > 1 ? 1 : 0 );
              }
              else {
                g.push_back(emx.tvcf.genos[j + offset] > 0 ? 1 : 0 );
              }
            }
            else {
              if ( arg.recessive ) {
                g.push_back(emx.tvcf.genos[j + offset] < 1 ? 1 : 0 );
              }
              else {
                g.push_back(emx.tvcf.genos[j + offset] < 2 ? 1 : 0 );
              }
            }
            //mac += g.back();
          }
          else {
            g.push_back(0);
            rg.push_back(0);
          }
        }

        ipass.push_back(i);
        if ( arg.scoref.empty() ) {
          macs.push_back(mac);
        }
        else {
          // use GERP/PhyloP score * 1000 as threshold
          macs.push_back((int)floor(gScore.baseScore(emx.tvcf.markers[i].c_str())*(-1000)+.5));
        }
      }
    }

    int mv = (int)macs.size(); // number of valid markers

    if ( mv > 0 ) { // at least one marker
      std::set<int> umacSet;          // contains the list of all MACs/SCOREs
      for(i=0; i < mv; ++i) {
        umacSet.insert(macs[i]);
      }

      std::vector<int> umacs;    // list of MAC/Score thresholds
      std::map<int,int> iumacs;  // iumacs[MAC] gives index
      for(std::set<int>::iterator it = umacSet.begin();
          it != umacSet.end(); ++it) {
        iumacs[*it] = umacs.size();
        umacs.push_back(*it);
      }
      int numacs = (int)umacs.size(); // number of unique MACs

      //notice("n = %d, mv = %d, numacs = %d",n, mv, numacs);
      //notice("g.size() = %d",(int)g.size());

      // generate ind * mac incidence matrix
      MatrixXd C = MatrixXd::Zero(n,numacs); // matrix of n * m'

      // this is 'collapsing-style' variable threshold test
      for(i=0, k=0; i < mv; ++i) {
        for(j=0; j < n; ++j, ++k) {
          if ( g[k] > 0 ) C(j,iumacs[macs[i]]) = 1.;
        }
      }

      // cumulate the incidence matrix
      for(i=1; i < numacs; ++i) {
        for(j=0; j < n; ++j) {
          if ( C(j,i-1) > 0 ) C(j,i) = 1.;
        }
      }

      // so C(j,i) = 1 if individual j carries RV with MAC <= umacs[i]
      VectorXd sumC = C.colwise().sum();
      VectorXd sqC = C.colwise().squaredNorm();

      // calculate correlation coefficient
      MatrixXd R = Yr * C; // (nperm+1) * numacs matrix
      int optMAC = 0;
      double optBurden = 0;
      double optR2 = 0;
      int optSgn = 0;
      int sgn;

      for(i=0; i < numacs; ++i) {
        for(j=0; j < arg.minperm+1; ++j) {
          sgn = R(j,i) > 0 ? 1 : -1;
          R(j,i) = R(j,i) * R(j,i) / ((sqC(i)-(double)sumC(i)*sumC(i)/n)*sqyr+pEmmaxHelper::ZEPS);
          if ( j == 0 ) {
            if ( R(j,i) > optR2 ) {
              optR2 = R(j,i);
              optBurden = sumC(i);
              optMAC = umacs[i];
              optSgn = sgn;
            }
          }
        }
      }

      // calculate max chi-squared statistics
      VectorXd maxR2 = R.rowwise().maxCoeff();

      // calculate the rank of p-value
      int rank = 0;
      for(i=1; i < arg.minperm+1; ++i) {
        if ( optR2 <= maxR2(i) ) ++rank;
      }

      int burdenMAC = sumC(numacs-1);
      int perm = arg.minperm;

      if ( rank * 10 < arg.minperm ) {
        for(perm=arg.minperm*10; perm <= arg.maxperm; perm *= 10) {
          int rperm = perm*9/10;
          notice("Performing adaptive permutation for additional %d permutations for group %s -- maxR2 = %lg, rank = %d",rperm,tokens[0].c_str(),optR2,rank);
          for(int p=0; p < rperm; p += arg.minperm) {
            for(i=0; i < arg.minperm; ++i) {
              std::random_shuffle(irand.begin(), irand.end());
              for(j=0; j < n; ++j) {
                Yp(i,j) = yr(irand[j]);
              }
            }
            R = Yp * C; // (nperm+1) * numacs matrix
            for(i=0; i < numacs; ++i) {
              for(j=0; j < arg.minperm; ++j) {
                R(j,i) = R(j,i) * R(j,i) / ((sqC(i)-(double)sumC(i)*sumC(i)/n)*sqyr+pEmmaxHelper::ZEPS);
              }
            }
            maxR2 = R.rowwise().maxCoeff();
            for(i=0; i < arg.minperm; ++i) {
              if ( optR2 <= maxR2(i) ) ++rank;
            }
          }
          if ( rank * 10 > arg.minperm ) break;
        }
        if ( perm > arg.maxperm ) perm = arg.maxperm;
      }
      wf.printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.5lf\t%.4lf\t%.4lg\t%.4lg\t%s\t%d\t%.5lf\n",emx.tvcf.groupChrom.c_str(),emx.tvcf.groupBeg,emx.tvcf.groupEnd,emx.tvcf.groupID.c_str(), m, mv, burdenMAC, (double)burdenMAC/n, optR2*n, (rank+1>perm ? 1. : (double)(rank+1.)/perm),optR2,optSgn > 0 ? "+" : "-",optMAC,(double)optBurden/n);
      //fprintf(stderr,"%s\t%d\t%d\t%d\t%.5lf\t%.4lf\t%.4lg\t%.4lf\t%.4lf\t%.5lf\n",tokens[0].c_str(), m, mv, 0, 0., 0., (rank+1.)/(arg.minperm),0.,0.,maxChisq(0));
      //wf.printf("%s\t%d\t%d\t%d\t%.5lf\t%.4lf\t%.4lg\t%.4lf\t%.4lf\t%.5lf\n",tokens[0].c_str(), m, mv, (int)nx, nx < n ? (double)nx/n : (double)(n-nx)/n,t,pval,beta,sebeta,r*r);
      //wf.printf("%s\t%d\t%d\t%d\t%.5lf\t%.4lf\t%.4lg\t%.4lf\t%.4lf\t%.5lf\n",tvcf.markerGroupName.c_str(), m, nPass, (int)nx, nx < n ? (double)nx/n : (double)(n-nx)/n,t,pval,beta,sebeta,r*r);

      //int iumac = iumacs[optMAC];
      for(i=0; i < (int)macs.size(); ++i) {
        if ( macs[i] <= optMAC ) {
          for(j=0; j < n; ++j) {
            if ( g[j + i*n] > 0 ) {
              vf.printf("%s\t%s\t%.5lf\t%d\t%d\t%s\t%s\n",
                emx.tvcf.groupID.c_str(),
                emx.tvcf.markers[ipass[i]].c_str(),
                emx.tvcf.alleleFreq(ipass[i]),
                macs[i],
                optMAC,
                emx.tvcf.inds[j].c_str(),
                genoLabels[rg[j + i * n]]);
            }
          }
        }
      }
    }
    else {
      wf.printf("%s\t%d\t%d\t%s\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n",emx.tvcf.groupChrom.c_str(),emx.tvcf.groupBeg,emx.tvcf.groupEnd,emx.tvcf.groupID.c_str(), m, mv);
      //wf.printf("%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n",tokens[0].c_str(), m, mv, 0);
      //wf.printf("%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n",tvcf.markerGroupName.c_str(), m, nPass, (int)nx);
    }
    //notice("bar");
  }
  return 0;
}

// variable threshold test
int runDump(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("markerset",&arg.markerset)
    LONG_STRINGPARAMETER("groupf",&arg.groupf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
    LONG_INTPARAMETER("minMAC",&arg.minMAC)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)
    LONG_PARAMETER("sepchr",&arg.sepchr)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("covf",&arg.covf)
    LONG_STRINGPARAMETER("scoref",&arg.scoref)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("outf",&arg.outf)
    LONG_PARAMETER("verbose",&arg.verbose)
    LONG_PARAMETER("summary",&arg.summarize)
    LONG_INTPARAMETER("bins",&arg.bins)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.vcf.empty() || arg.phenof.empty() || arg.markerset.empty() || arg.outf.empty() ) {
    error("--vcf, --phenof, --markerset, --outf are required parameters");
  }

  pEmmax emx;
  emx.loadFiles(arg.phenof.c_str(), arg.covf.c_str(), arg.indf.c_str(), NULL, NULL, arg.vcf.c_str(), arg.region.c_str(), NULL, arg.field.c_str(), !arg.ignoreFilter);

  // categorize phenotypes into categories
  std::map<double,int> phenoCnts;
  int n = (int)emx.inds.size();

  std::vector<double> mins, maxs;
  std::vector<int> bins, cnts;
  int nbins = 0;

  if ( arg.summarize ) {
    for(int i=0; i < n; ++i) {
      ++phenoCnts[emx.y(i)];
    }

    if ( phenoCnts.size() < 10 ) {
      notice("Found %d distinct phenotype values. Considering them categorical values",(int)phenoCnts.size());
      std::map<double,int> pheno2bin;
      for(std::map<double,int>::iterator it = phenoCnts.begin();
          it != phenoCnts.end(); ++it) {
        pheno2bin[it->first] = (int)cnts.size();
        mins.push_back(it->first);
        maxs.push_back(it->first);
        cnts.push_back(it->second);
      }
      for(int i=0; i < n; ++i) {
        bins.push_back(pheno2bin[emx.y(i)]);
      }
      nbins = (int)mins.size();
    }
    else {
      notice("Found >10 distinct phenotype values. Categorizing into %d bins",arg.bins);

      nbins = arg.bins;
      int nsize = emx.y.size();
      int cumCnts = 0;
      int curBin = 0;
      for(std::map<double,int>::iterator it = phenoCnts.begin();
          it != phenoCnts.end(); ++it) {
        if ( curBin == (int)mins.size() ) {
          mins.push_back(it->first);
          maxs.push_back(it->first);
          cnts.push_back(it->second);
        }
        else {
          maxs.back() = it->first;
          cnts.back() += it->second;
        }
        cumCnts += it->second;

        if ( cumCnts > (double)(curBin + 1.)/nbins*nsize ) { ++curBin; }
      }

      if ( nbins != (int)mins.size() ) {
        error("Total number of bins is %d when expecting %d",(int)mins.size(),nbins);
      }

      //notice("foo");

      int i,j;
      bins.resize(n,0);
      for(i=0; i < n; ++i) {
        for(j=0; j < nbins; ++j) {
          if ( ( emx.y(i) >= mins[j] ) && ( emx.y(i) <= maxs[j] ) ) {
            bins[i] = j;
            break;
          }
        }
        if ( j == nbins ) {
          error("Cannot categorize phenotype value %lf",emx.y(i));
        }
      }

      if ( n != (int)bins.size() ) {
        error("Total number of instanced is %d when expecting %d",(int)bins.size(),n);
      }
    }
  }

  std::vector<std::string> tokens;
  std::vector<float> genos;
  std::vector<float> afs;
  if ( arg.groupf.empty() ) {
    pFile::tokenizeLine(arg.markerset.c_str()," \t\r\n",tokens);
  }
  else {
    pFile tgrp(arg.groupf.c_str());
    const char* line = NULL;
    while( (line = tgrp.getLine()) != NULL ) {
      std::vector<std::string> tokens2;
      if ( line[0] == '#' ) continue;
      pFile::tokenizeLine(line," \t\r\n", tokens2);
      if ( tokens2[0] == arg.markerset ) {
        for(int i=1; i < (int)tokens2.size(); ++i) {
          tokens.push_back(tokens2[i]);
        }
        break;
      }
    }
    if ( tokens.size() == 0 ) {
      error("Cannot find the markerset ID %s from file %s",arg.markerset.c_str(),arg.groupf.c_str());
    }
  }

  emx.tvcf.clear();
  int m = emx.tvcf.readMarkerGroup(tokens,0,arg.sepchr);

  int i, j, nPass = 0;
  std::vector<std::string> passMarkers;
  for(i=0; i < m; ++i) {
    if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) &&
         ( emx.tvcf.MAF(i) <= arg.maxMAF ) &&
         ( emx.tvcf.MAC(i) >= arg.minMAC ) &&
         ( emx.tvcf.callRate(i) >= arg.minCallRate ) &&
         ( ( arg.minRSQ == 0 ) || ( emx.tvcf.RSQ(i) >= arg.minRSQ ) ) ) {
      for(j=0; j < n; ++j) {
        genos.push_back(emx.tvcf.genos[j + n*i]);
      }
      afs.push_back(emx.tvcf.alleleFreq(i));
      passMarkers.push_back(tokens[i]);
      ++nPass;
      //notice("i = %d, m= %d",i,m);
    }
  }

  notice("Found %d of %d markers matching the criteria acroess %d samples - %d genotypes",nPass,m,n,(int)genos.size());

  //notice("Writing to %s -- %d markers, %d covs",arg.outf.c_str(),(int)passMarkers.size(),emx.ncovs);
  wFile wf(arg.outf.c_str());
  if ( arg.summarize ) {
    genomeScore gScore;
    if ( !arg.scoref.empty() ) {
      gScore.setDir(arg.scoref.c_str());
    }

    wf.printf("---------------------------------------------------------------------------------\n");
    wf.printf("#MARKER_ID\tMAF");
    if ( !arg.scoref.empty() ) { wf.printf("\tSCORE"); }
    for(i=0; i < nbins; ++i) {
      if ( mins[i] == maxs[i] ) {
        wf.printf("\t%lg (%d)",mins[i],cnts[i]);
      }
      else {
        wf.printf("\t[%lg,%lg] (%d)",mins[i],maxs[i],cnts[i]);
      }
    }
    wf.printf("\n");
    wf.printf("---------------------------------------------------------------------------------\n");

    std::vector<int> collapse(n,0);
    for(i=0; i < nPass; ++i) {
      wf.printf("%s",passMarkers[i].c_str());
      wf.printf("\t%.5lf",afs[i]);
      if ( !arg.scoref.empty() ) {
        wf.printf("\t%lg",gScore.baseScore(passMarkers[i].c_str()));
      }

      std::vector<int> cnts(nbins*3,0);
      for(j=0; j < n; ++j) {
        float g = genos[j+n*i];
        if ( !std::isnan(g) ) {
          ++cnts[g + bins[j]*3];
          if (g != (( afs[i] < 0.5 ) ? 0. : 2.)) collapse[j] = 1;
        }
      }

      for(j=0; j < nbins; ++j) {
        wf.printf("\t%d,%d,%d",cnts[j*3],cnts[j*3+1],cnts[j*3+2]);
      }
      wf.printf("\n");
    }

    int ncollapse = 0;
    std::vector<int> bincollapse(nbins,0);
    for(i=0; i < n; ++i) {
      if ( collapse[i] ) {
        ++ncollapse;
        ++bincollapse[bins[i]];
      }
    }

    wf.printf("---------------------------------------------------------------------------------\n");
    wf.printf("FRAC_WITH_RARE:\t%.5lf",(double)ncollapse/n);
    if ( !arg.scoref.empty() ) { wf.printf("\tNA"); }
    for(i=0; i < nbins; ++i) {
      wf.printf("\t%d,%d",cnts[i]-bincollapse[i],bincollapse[i]);
    }
    wf.printf("\n");
    wf.printf("---------------------------------------------------------------------------------\n");
  }
  else {
    wf.printf("#IND_ID\tPHENO");
    for(i=0; i < emx.ncovs; ++i) {
      wf.printf("\tCOV%d",i+1);
    }
    for(i=0; i < nPass; ++i) {
      wf.printf("\t%s",passMarkers[i].c_str());
    }
    wf.printf("\n");

    for(i=0; i < n; ++i) {
      wf.printf("%s",emx.inds[i].c_str());
      wf.printf("\t%lg",emx.y(i));
      for(j=0; j < emx.ncovs; ++j) {
        wf.printf("\t%lg",emx.X(i,j));
      }
      for(j=0; j < (int)nPass; ++j) {
        float g = genos[i + j*n];
        if ( std::isnan(g) ) {
          wf.printf("\tNA");
        }
        else {
          wf.printf("\t%lf",g);
        }
      }
      wf.printf("\n");
    }
  }
  wf.close();
  return 0;
}

// variable threshold test
int runVcfInfo(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  arg.field = "PL";
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_INTPARAMETER("unit",&arg.unit)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("scoref",&arg.scoref)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("outf",&arg.outf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.vcf.empty() || arg.outf.empty() ) {
    error("--vcf, --outf are required parameters");
  }

  pEmmax emx;
  emx.loadFiles(NULL, NULL, arg.indf.c_str(), NULL, NULL, arg.vcf.c_str(), arg.region.c_str(), NULL, arg.field.c_str(), !arg.ignoreFilter);

  wFile wf(arg.outf.c_str());
  wf.printf("#CHROM\tBEG\tEND\tMARKER_ID\tNS\tAF\tAC\tCALLRATE\tRSQ\tSCORE\n");

  genomeScore gScore;
  if ( !arg.scoref.empty() ) {
    gScore.setDir(arg.scoref.c_str());
  }

  for(int M=0; emx.tvcf.readMarkers(arg.unit); ) {
    M += emx.tvcf.nMarkers;
    fprintf(stderr,"Processing %d markers across %d individuals...",M, emx.tvcf.nInds);
    for(int i=0; i < emx.tvcf.nMarkers; ++i) {
      wf.printf("%s\t%d\t%d\t%s\t%d\t%.5lf\t%.2lf\t%.5lf\t%.5lf",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(),emx.tvcf.numAlleles[i]/2,emx.tvcf.alleleFreq(i),emx.tvcf.sumAlleles[i],emx.tvcf.callRate(i),emx.tvcf.RSQ(i));
      if ( arg.scoref.empty() ) {
        wf.printf("\tNA\n");
      }
      else {
        wf.printf("\t%.4lf\n",gScore.baseScore(emx.tvcf.markers[i].c_str()));
      }
    }
  }
  wf.close();
  return 0;
}

int runKinUtil(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  bool txtkin = false;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Options")
    LONG_STRINGPARAMETER("kinf",&arg.kinf)
    LONG_PARAMETER("text",&txtkin)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("outf",&arg.outf)
    LONG_INTPARAMETER("pca",&arg.pca)
    LONG_DOUBLEPARAMETER("relpair",&arg.relpair)
    LONG_PARAMETER("maxkin",&arg.maxkin)
    LONG_PARAMETER("dump",&arg.dump)
    LONG_INTPARAMETER("digits",&arg.digits)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  //int argstart = pl.ReadWithTrailer(argc, argv) + 1;
  pl.Read(argc, argv);
  pl.Status();

  if ( arg.outf.empty() || arg.kinf.empty() ) {
    error("ERROR: --outf and --kinf parameters are required");
  }

  if ( ( arg.pca == 0 ) && ( arg.relpair == 1 ) && ( !arg.dump ) && ( arg.maxkin == 0 ) ) {
    error("ERROR: one of --pca , --relpair, --dumpm, --maxkin are required");
  }
  
  MatrixXd K;
  int flag;
  double wsum;
  std::vector<std::string> ids;

  // read kinship matrix
  notice("Reading kinship matrix file %s",arg.kinf.c_str());
  if ( txtkin ) {
    pEmmaxHelper::readTxtKinWithIDs(arg.kinf.c_str(), K, flag, wsum, ids);
  }
  else {
    pEmmaxHelper::readKinWithIDs(arg.kinf.c_str(), K, flag, wsum, ids);
  }

  notice("Successfully loaded kinship matrix for %d individuals", (int)ids.size());

  // first, compute PCA if needed
  if ( arg.pca > 0 ) {
    notice("Computing %d principal components", arg.pca);

    int n = K.rows();
    int p = arg.pca;
    SelfAdjointEigenSolver<MatrixXd> eig(K);

    notice("Writing eigenvectors to %s.evecs and eigenvalues to %s.evals", arg.outf.c_str(), arg.outf.c_str());

    VectorXd v = eig.eigenvalues();
    VectorXd evals(p);
    for(int i=n-p; i < n; ++i) {
      evals(i-n+p) = v(i);
    }

    MatrixXd evecs = eig.eigenvectors().block(0,n-p,n,p);

    std::string tmp = arg.outf + ".evecs";
    wFile wf(tmp.c_str());
    wf.printf("#IND_ID");
    for(int i=0; i < p; ++i) {
      wf.printf("\tPC%d",i+1);
    }
    wf.printf("\n");

    for(int i=0; i < n; ++i) {
      wf.printf("%s",ids[i].c_str());
      for(int j=0; j < p; ++j) {
	wf.printf("\t%.*lf",arg.digits,evecs(i,p-j-1));
      }
      wf.printf("\n");
    }
    wf.close();

    tmp = arg.outf + ".evals";
    wFile wf2(tmp.c_str());
    wf2.printf("PC\tEIGENVALUE\n");
    for(int i=0; i < p; ++i) {
      wf2.printf("PC%d\t%.*lf\n",i+1,arg.digits,evals(p-i-1));
    }
    wf2.close();
  }

  if ( arg.relpair < 1 ) {
    notice("Finding individuals with kinship coefficient > %lf", arg.relpair);

    int n = K.rows();
    std::string tmp = arg.outf + ".relpair";
    wFile wf(tmp.c_str());
    for(int i=0; i < n; ++i) {
      for(int j=0; j < i; ++j) {
	if ( K(j,i) > arg.relpair ) {
	  wf.printf("%s\t%s\t%.*lf\n",ids[j].c_str(),ids[i].c_str(),arg.digits,K(j,i));
	}
      }
    }
    wf.close();
  }

  if ( arg.dump ) {
    if ( txtkin ) {
      std::string tmp = arg.outf + ".kinf";
      notice("Dumping the entire kinship matrix into a binary file %s",tmp.c_str());
      pEmmaxHelper::writeKinWithIDs(tmp.c_str(), K, wsum, ids);
    }
    else {
      std::string tmp = arg.outf + ".kindump";
      notice("Dumping the entire kinship matrix into a text file %s",tmp.c_str());
      
      int n = K.rows();
      wFile wf(tmp.c_str());
      for(int i=0; i < n; ++i) {
	if ( i > 0 ) wf.printf("\t");
	wf.printf("%s",ids[i].c_str());
      }
      wf.printf("\n");
      
      for(int i=0; i < n; ++i) {
	wf.printf("%s",ids[i].c_str());
	for(int j=0; j < n; ++j) {
	  wf.printf("\t%.*lf",arg.digits,K(j,i));
	}
	wf.printf("\n");
      }
      wf.close();
    }
  }

  if ( arg.maxkin ) {
    std::string tmp = arg.outf + ".maxkin";
    notice("Finding the closest individual for each individual");

    int n = K.rows();
    std::vector<double> maxkin(n,0);
    std::vector<int> imaxkin(n,0);
    for(int i=0; i < n; ++i) {
      for(int j=0; j < i; ++j) {
	if ( K(j,i) > maxkin[i] ) {
	  maxkin[i] = K(j,i);
	  imaxkin[i] = j;
	}
	if ( K(i,j) > maxkin[j] ) {
	  maxkin[j] = K(j,i);
	  imaxkin[j] = i;
	}
      }
    }

    wFile wf(tmp.c_str());
    wf.printf("#IND_ID\tCLOSEST_IND\tKINCOEFF\n");
    for(int i=0; i < n; ++i) {
      wf.printf("%s\t%s\t%.*lf\n",ids[i].c_str(),ids[imaxkin[i]].c_str(),arg.digits,maxkin[i]);
    }
    wf.close();
  }
  return 0;
}

int runMultiAssoc(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_INTPARAMETER("unit",&arg.unit)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
    LONG_INTPARAMETER("minMAC",&arg.minMAC)
    LONG_INTPARAMETER("maxMAC",&arg.maxMAC)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("tryf",&arg.tryf)
    LONG_STRINGPARAMETER("eigf",&arg.ineigf)
    LONG_STRINGPARAMETER("remlf",&arg.remlf)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_PARAMETER("compact",&arg.compact)
    LONG_DOUBLEPARAMETER("maxP",&arg.maxP)
    LONG_STRINGPARAMETER("out-assocf",&arg.assocf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.tryf.empty() || arg.ineigf.empty() || arg.remlf.empty() ) {
    error("--tryf, --remlf and --kinf are required parameters");
  }

  pEmmaxMulti emx;
  emx.normalize = arg.normalize;
  emx.loadFiles(arg.tryf.c_str(), NULL, arg.indf.c_str(), NULL, arg.ineigf.c_str(), arg.vcf.c_str(), arg.region.c_str(), arg.rule.c_str(), arg.field.c_str(), !arg.ignoreFilter);

  pFile treml(arg.remlf.c_str());
  std::vector<double> deltas;
  std::vector<std::string> tokens;
  const char* line = NULL;
  while( (line = treml.getLine()) != NULL ) {
    if ( line[0] == '#' ) continue;
    pFile::tokenizeLine(line," \t\r\n", tokens);
    if ( tokens.size() != 7 )
      error("Cannot parse REML line %s. Must have seven columns",line);
    if ( emx.pedcols[(int)deltas.size()] != tokens[0] ) {
      error("REML file and phenotype file does not match. %s observed, %s expected",tokens[0].c_str(), emx.pedcols[(int)deltas.size()].c_str());
    }
    deltas.push_back(atof(tokens[1].c_str()));
  }

  int c = emx.evalR.rows();
  int g = emx.Y.cols();
  int n = emx.inds.size();
  MatrixXd Ur = emx.evecR.transpose();

  wFile wf(arg.assocf.c_str());
  double mu = 0;
  int i,m,j,k;
  double sx = 0, sxx = 0, sxy = 0, r = 0, t = 0, beta = 0, pval = 0; //, sebeta = 0, varE = 0;

  // calculate sufficient statistics for var(y)
  std::vector<double> sys, syys;
  for(int i=0; i < g; ++i) {
    double sy = 0, syy = 0, tmp;
    for(int j=0; j < c; ++j) {
      tmp = emx.Y(j,i) / sqrt(deltas[i] + emx.evalR(j));
      sy += tmp;
      syy += (tmp * tmp);
    }
    sys.push_back(sy);
    syys.push_back(syy);
  }

  int genocnts[3] = {0,0,0};
  VectorXd x(n), xt(n);
  wf.printf("#CHROM\tBEG\tEND\tMARKER_ID\tNS\tAC\tCALLRATE\tGENOCNT\tMAF");
  if ( !arg.compact ) {
    for(i=0; i < g; ++i) {
      wf.printf("\t%s.P\t%s.B",emx.pedcols[i].c_str(),emx.pedcols[i].c_str());
    }
  }
  wf.printf("\n");

  for(int M=0; emx.tvcf.readMarkers(arg.unit); ) {
    M += emx.tvcf.nMarkers;
    notice("Processing %d markers across %d individuals...", M, emx.tvcf.nInds);
    for(i=0, m=0; i < emx.tvcf.nMarkers; ++i) {
      if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) &&
           ( emx.tvcf.MAF(i) <= arg.maxMAF ) &&
           ( emx.tvcf.MAC(i) >= arg.minMAC ) &&
           ( emx.tvcf.MAC(i) <= arg.maxMAC ) &&
           ( emx.tvcf.callRate(i) >= arg.minCallRate ) &&
           ( ( arg.minRSQ == 0 ) || ( emx.tvcf.RSQ(i) >= arg.minRSQ ) ) ) {

        //fprintf(stderr,"%d\n",i);
        int offset = i * emx.tvcf.nInds;
        mu = emx.tvcf.alleleFreq(i)*2.;
        for(j=0; j < emx.tvcf.nInds; ++j) {
          if ( std::isnan(emx.tvcf.genos[j + offset]) ) {
            x(j) = mu;
          }
          else {
            x(j) = emx.tvcf.genos[j + offset];
          }
        }

        emx.tvcf.GENOCNT(i,genocnts);

        // CHROM BEG END MARKER_ID NS AC CALLRATE GENOCNT MAF
        wf.printf("%s\t%d\t%d\t%s\t%d\t%.2lf\t%.5lf\t%d/%d/%d\t%.5lf",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(), emx.tvcf.numAlleles[i]/2, emx.tvcf.sumAlleles[i], emx.tvcf.callRate(i), genocnts[0], genocnts[1], genocnts[2], emx.tvcf.MAF(i));

        xt = Ur * x; // size q
        for(k=0; k < g; ++k) {
          sx = sxx = sxy = 0;
          //xy = xsum = xsq = 0;
          for(j=0; j < c; ++j) {
            sxy += (xt(j)*emx.Y(j,k)/(deltas[k] + emx.evalR(j)));
            sx += (xt(j)/sqrt(deltas[k] + emx.evalR(j)));
            sxx += ((xt(j)*xt(j))/(deltas[k] + emx.evalR(j)));
          }
          // this is the estimates including the intercept always
          beta = ((c+1)*sxy-sx*sys[k])/((c+1)*sxx-sx*sx);
          //varE = 1/(c+1.)/(c-1.)*((c+1)*syys[k]-sys[k]*sys[k]-beta*beta*((c+1)*sxx-sx*sx));
          //sebeta = sqrt((c+1)*varE/((c+1)*sxx-sx*sx));
          r = ((c+1)*sxy-sx*sys[k])/sqrt(((c+1)*sxx-sx*sx)*((c+1)*syys[k]-sys[k]*sys[k]));
          t = r * sqrt((c-1)/(1-r*r+pEmmaxHelper::ZEPS));
          pval = pEmmaxHelper::tcdf(t, c-1);
          if ( pval <= arg.maxP ) {
            if ( arg.compact ) {
              wf.printf("\t%s,%.3lg,%.3lg",emx.pedcols[k].c_str(),pval,beta);
            }
            else {
              wf.printf("\t%.3lg\t%.3lg",pval,beta);
            }
          }
          else if ( !arg.compact ) {
            wf.printf("\tNA\tNA");
          }
        }
        wf.printf("\n");
        ++m;
      }
      else if ( ! arg.compact ) {
        mu = emx.tvcf.alleleFreq(i)*2.;
        emx.tvcf.GENOCNT(i,genocnts);
        wf.printf("%s\t%d\t%d\t%s\t%d\t%.2lf\t%.5lf\t%d/%d/%d\t%.5lf",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(), emx.tvcf.numAlleles[i]/2, emx.tvcf.sumAlleles[i], emx.tvcf.callRate(i), genocnts[0], genocnts[1], genocnts[2], emx.tvcf.MAF(i));
        for(k=0; k < g; ++k) {
          wf.printf("\tNA\tNA");
        }
        wf.printf("\n");
      }
    }
    fprintf(stderr,"%d markers included\n",m);

  }
  return 0;
}

// Perform REML fitting
int runMultiReml(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Options")
    LONG_STRINGPARAMETER("phenof",&arg.phenof)
    LONG_STRINGPARAMETER("covf",&arg.covf)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_PARAMETER("no-intercept",&arg.noIntercept)
    LONG_STRINGPARAMETER("kinf",&arg.kinf)
    LONG_PARAMETER("normalize",&arg.normalize)
    LONG_STRINGPARAMETER("in-eigf",&arg.ineigf)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_STRINGPARAMETER("out-eigf",&arg.outeigf)
    LONG_STRINGPARAMETER("out-remlf",&arg.remlf)
    LONG_STRINGPARAMETER("out-tryf",&arg.tryf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.phenof.empty() || arg.kinf.empty() || arg.remlf.empty() ) {
    error("--phenof, --kinf, and --out-remlf are required parameters");
  }

  pEmmaxMulti emx;
  emx.normalize = arg.normalize;

  // load a set of phenotypes, covariates, and kinship, (and eigendecomposition)
  emx.loadFiles(arg.phenof.c_str(), arg.covf.c_str(), arg.indf.c_str(), arg.kinf.c_str(), arg.ineigf.c_str(), arg.vcf.c_str(), NULL, NULL, "GT", true);

  if ( !arg.ineigf.empty() && pEmmaxHelper::exists(arg.ineigf.c_str()) ) { 
    notice("Loading eigendecomposition matrix %s",arg.ineigf.c_str());
  }
  else {
    notice("Computing eigendecomposition of residual kinship matrix with size %d x %d, predictor matrix %d x %d",emx.K.rows(),emx.K.cols(),emx.X.rows(),emx.X.cols());
    pEmmaxHelper::computeEigenRestricted(emx.K, emx.X, emx.evecR, emx.evalR, emx.trK); // compute eigL
    if ( ! arg.outeigf.empty() ) {
      notice("Writing eigendecomposition into kinship matrix",arg.outeigf.c_str());
      pEmmaxHelper::writeEigenWithIDs(arg.outeigf.c_str(), emx.evecR, emx.evalR, emx.trK, emx.inds); // save eigL
    }
  }

  if ( !arg.phenof.empty() ) {
    MatrixXd U = emx.evecR.transpose() * emx.Y; // (n-p) * g matrix
    int r = U.rows();
    int c = U.cols();
    wFile yf(arg.tryf.c_str());

    // print out header information of transposed phenotypes
    yf.printf("#IND_ID");
    for(int i=0; i < (int)emx.pedcols.size(); ++i) {
      yf.printf("\t%s",emx.pedcols[i].c_str());
    }
    yf.printf("\n");

    for(int i=0; i < (int)emx.inds.size(); ++i) {  // although the transformed space do not involve individual IDs, write them down!
      yf.printf("%s",emx.inds[i].c_str());
      if ( i < r ) {
	for(int j=0; j < c; ++j) {
	  yf.printf("\t%lg",U(i,j));
	}
      }
      else {
	for(int j=0; j < c; ++j) {
	  yf.printf("\t0");
	}
      }
      yf.printf("\n");
    }
    yf.close();

    // calculate optimal deltas for each traits
    wFile wf(arg.remlf.c_str());
    wf.printf("#TRAIT\tDELTA\tV_G\tV_E\tLLK0\tLLK1\tH2\n");
    for(int k=0; k < c; ++k) {
      REML reml = pEmmaxHelper::computeREMLU(U.col(k), emx.evecR, emx.evalR, emx.trK);
      wf.printf("%s",emx.pedcols[k].c_str());
      wf.printf("\t%lg",reml.delta);
      wf.printf("\t%lg",reml.vg);
      wf.printf("\t%lg",reml.ve);
      wf.printf("\t%lg",reml.LLK0);
      wf.printf("\t%lg",reml.LLK1);
      wf.printf("\t%lg\n",reml.h2);
    }
    wf.close();
  }
  return 0;
}

int runMultiAssocPlain(int argc, char** argv) {
  // Parse the input arguments
  pEmmaxArgs arg;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("VCF Input Options")
    LONG_STRINGPARAMETER("vcf",&arg.vcf)
    LONG_STRINGPARAMETER("region",&arg.region)
    LONG_INTPARAMETER("unit",&arg.unit)
    LONG_STRINGPARAMETER("indf",&arg.indf)
    LONG_STRINGPARAMETER("field",&arg.field)
    LONG_STRINGPARAMETER("rule",&arg.rule)
    LONG_DOUBLEPARAMETER("minMAF",&arg.minMAF)
    LONG_DOUBLEPARAMETER("maxMAF",&arg.maxMAF)
    LONG_INTPARAMETER("minMAC",&arg.minMAC)
    LONG_INTPARAMETER("maxMAC",&arg.maxMAC)
    LONG_DOUBLEPARAMETER("minCallRate",&arg.minCallRate)
    LONG_DOUBLEPARAMETER("minRSQ",&arg.minRSQ)
    LONG_PARAMETER("ignoreFilter",&arg.ignoreFilter)

    LONG_PARAMETER_GROUP("Other Input Options")
    LONG_STRINGPARAMETER("pheno",&arg.phenof)

    LONG_PARAMETER_GROUP("Output Options")
    LONG_PARAMETER("compact",&arg.compact)
    LONG_DOUBLEPARAMETER("maxP",&arg.maxP)
    LONG_STRINGPARAMETER("out-assocf",&arg.assocf)
    LONG_PARAMETER("verbose",&arg.verbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options", longParameters));
  pl.Read(argc,argv);
  pl.Status();

  // sanity check of input arguments
  if ( arg.phenof.empty() && arg.vcf.empty() ) {
    error("--pheno, --vcf are required parameters");
  }

  pEmmaxMulti emx;
  emx.normalize = arg.normalize;
  emx.loadFiles(arg.phenof.c_str(), NULL, arg.indf.c_str(), NULL, NULL, arg.vcf.c_str(), arg.region.c_str(), arg.rule.c_str(), arg.field.c_str(), !arg.ignoreFilter);

  int g = emx.Y.cols();
  int n = emx.inds.size();

  wFile wf(arg.assocf.c_str());
  double mu = 0;
  int i,m,j,k;
  double sx = 0, sxx = 0, sxy = 0, r = 0, t = 0, beta = 0, pval = 0; //, sebeta = 0, varE = 0;

  // calculate sufficient statistics for var(y)
  std::vector<double> sys, syys;
  for(int i=0; i < g; ++i) {
    double sy = 0, syy = 0, sum = 0, tmp;
    for(int j=0; j < n; ++j) {
      sum += emx.Y(j,i);
    }
    sum /= (double)n;
    for(int j=0; j < n; ++j) {
      emx.Y(j,i) -= sum;
      tmp = emx.Y(j,i);
      sy += tmp;
      syy += (tmp * tmp);
    }
    sys.push_back(sy);
    syys.push_back(syy);
  }

  int genocnts[3] = {0,0,0};
  VectorXd x(n), xt(n);
  wf.printf("#CHROM\tBEG\tEND\tMARKER_ID\tNS\tAC\tCALLRATE\tGENOCNT\tMAF");
  if ( !arg.compact ) {
    for(i=0; i < g; ++i) {
      wf.printf("\t%s.P\t%s.B",emx.pedcols[i].c_str(),emx.pedcols[i].c_str());
    }
  }
  wf.printf("\n");

  for(int M=0; emx.tvcf.readMarkers(arg.unit); ) {
    M += emx.tvcf.nMarkers;
    notice("Processing %d markers across %d individuals...", M, emx.tvcf.nInds);
    for(i=0, m=0; i < emx.tvcf.nMarkers; ++i) {
      if ( ( emx.tvcf.MAF(i) >= arg.minMAF ) &&
           ( emx.tvcf.MAF(i) <= arg.maxMAF ) &&
           ( emx.tvcf.MAC(i) >= arg.minMAC ) &&
           ( emx.tvcf.MAC(i) <= arg.maxMAC ) &&
           ( emx.tvcf.callRate(i) >= arg.minCallRate ) &&
           ( ( arg.minRSQ == 0 ) || ( emx.tvcf.RSQ(i) >= arg.minRSQ ) ) ) {

        //fprintf(stderr,"%d\n",i);
        int offset = i * emx.tvcf.nInds;
        mu = emx.tvcf.alleleFreq(i)*2.;
        for(j=0; j < emx.tvcf.nInds; ++j) {
          if ( std::isnan(emx.tvcf.genos[j + offset]) ) {
            x(j) = 0;
          }
          else {
            x(j) = emx.tvcf.genos[j + offset] - mu;
          }
        }

        emx.tvcf.GENOCNT(i,genocnts);

        // CHROM BEG END MARKER_ID NS AC CALLRATE GENOCNT MAF
        wf.printf("%s\t%d\t%d\t%s\t%d\t%.2lf\t%.5lf\t%d/%d/%d\t%.5lf",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(), emx.tvcf.numAlleles[i]/2, emx.tvcf.sumAlleles[i], emx.tvcf.callRate(i), genocnts[0], genocnts[1], genocnts[2], emx.tvcf.MAF(i));

        for(k=0; k < g; ++k) {
          sx = sxx = sxy = 0;
          //xy = xsum = xsq = 0;
          for(j=0; j < n; ++j) {
            sxy += x(j) * emx.Y(j,k);
            sx += x(j);
            sxx += (x(j)*x(j));
          }
          // this is the estimates including the intercept always
          beta = (n*sxy-sx*sys[k])/(n*sxx-sx*sx);
          //varE = 1/(c+1.)/(c-1.)*((c+1)*syys[k]-sys[k]*sys[k]-beta*beta*((c+1)*sxx-sx*sx));
          //sebeta = sqrt((c+1)*varE/((c+1)*sxx-sx*sx));
          r = (n*sxy-sx*sys[k])/sqrt((n*sxx-sx*sx)*(n*syys[k]-sys[k]*sys[k]));
          t = r * sqrt((n-2)/(1-r*r+pEmmaxHelper::ZEPS));
          pval = pEmmaxHelper::tcdf(t, n-2);
          if ( pval <= arg.maxP ) {
            if ( arg.compact ) {
              wf.printf("\t%s,%.3lg,%.3lg",emx.pedcols[k].c_str(),pval,beta);
            }
            else {
              wf.printf("\t%.3lg\t%.3lg",pval,beta);
            }
          }
          else if ( !arg.compact ) {
            wf.printf("\tNA\tNA");
          }
        }
        wf.printf("\n");
        ++m;
      }
      else if ( ! arg.compact ) {
        mu = emx.tvcf.alleleFreq(i)*2.;
        emx.tvcf.GENOCNT(i,genocnts);
        wf.printf("%s\t%d\t%d\t%s\t%d\t%.2lf\t%.5lf\t%d/%d/%d\t%.5lf",emx.tvcf.chroms[i].c_str(),emx.tvcf.pos1s[i],emx.tvcf.pos1s[i]+emx.tvcf.refs[i].size()-1,emx.tvcf.markers[i].c_str(), emx.tvcf.numAlleles[i]/2, emx.tvcf.sumAlleles[i], emx.tvcf.callRate(i), genocnts[0], genocnts[1], genocnts[2], emx.tvcf.MAF(i));
        for(k=0; k < g; ++k) {
          wf.printf("\tNA\tNA");
        }
        wf.printf("\n");
      }
    }
    fprintf(stderr,"%d markers included\n",m);

  }
  return 0;
}


int main(int argc, char** argv) {
  if ( argc < 2 ) {
    printf("EMMAX for VCF v1.0.0\n");
    printf("Usage : %s [gen-kin | merge-kin | reml | assoc | all] [options]\n",argv[0]);
    printf("\tType one of the following commands below to get detailed usage\n");
    printf("\t%s gen-kin     [options] : create a kinship matrix\n",argv[0]);
    printf("\t%s merge-kin   [options] : merge multiple kinship matrices\n",argv[0]);
    printf("\t%s reml        [options] : fit the phenotypes to mixed model via REML\n",argv[0]);
    printf("\t%s assoc       [options] : run EMMAX association analysis\n",argv[0]);
    printf("\t%s burden-assoc [options] : run EMMAX association analysis\n",argv[0]);
    printf("\t%s remlassoc   [options] : fit REML & run EMMAX association analysis\n",argv[0]);
    printf("\t%s GLRT        [options] : simple genotype likelihood ratio test for binary phenotypes\n",argv[0]);
    printf("\t%s simul       [options] : simulate phenotypes with polygenic background\n",argv[0]);
    printf("\t%s all         [options] : perform all-in-one association analysis sequentially\n",argv[0]);
    return -1;
  }

  std::string cmd(argv[1]);
  if ( cmd == "gen-kin" ) {
    return runGenKin(argc-1,argv+1);
  }
  else if ( cmd == "join" ) {
    return runGenKin(argc-1,argv+1);
  }
  else if ( cmd == "merge-kin" ) {
    return runMergeKin(argc-1,argv+1);
  }
  else if ( cmd == "reml" ) {
    return runReml(argc-1,argv+1);
  }
  else if ( cmd == "assoc" ) {
    return runAssoc(argc-1,argv+1);
  }
  else if ( cmd == "burden-assoc" ) {
    return runBurdenAssoc(argc-1,argv+1);
  }
  else if ( cmd == "simul" ) {
    return runSimul(argc-1,argv+1);
  }
  else if ( cmd == "GLRT" ) {
    return runGLRT(argc-1,argv+1);
  }
  else if ( cmd == "VT" ) {
    return runVT(argc-1,argv+1);
  }
  else if ( cmd == "MMVT" ) {
    return runEmmaxVT(argc-1,argv+1);
  }
  else if ( cmd == "dump") {
    return runDump(argc-1,argv+1);
  }
  else if ( cmd == "vcfinfo") {
    return runVcfInfo(argc-1,argv+1);
  }
  else if ( cmd == "kin-util") {
    return runKinUtil(argc-1,argv+1);
  }
  else if ( cmd == "all" ) {
    return runAll(argc,argv);
  }
  else if ( cmd == "multi-reml") {
    return runMultiReml(argc-1,argv+1);
  }
  else if ( cmd == "multi-assoc") {
    return runMultiAssoc(argc-1,argv+1);
  }
  else if ( cmd == "multi-assoc-plain") {
    return runMultiAssocPlain(argc-1,argv+1);
  }
  else {
    error("Unrecognized command %s\n",argv[1]);
  }
}
