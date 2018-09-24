#include <iostream>
#include <vector>
#include <map>
#include "fVcf.h"
#include "pFile.h"
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <cmath>

#define N_MAGIC 8
#define MAGIC_KIN "EMMA_KIN"
#define MAGIC_EIG "EMMA_EIG"

// perform tabixed vTest
extern "C" {
  char* strAllocCopy(SEXP s, int idx = 0) {
    PROTECT(s = AS_CHARACTER(s));
    char* p = R_alloc(strlen(CHAR(STRING_ELT(s,idx))), sizeof(char));
    strcpy(p, CHAR(STRING_ELT(s, idx)));    
    UNPROTECT(1);
    return(p);
  }

  // write matrix to a file. same to write.table, but write a header in a different mannaer
  SEXP writeMatrix(SEXP matrix, SEXP outname) {
    //fprintf(stderr,"writeMatrix() called\n");
    const char* soutf = strAllocCopy(outname);
    FILE* fp = fopen(soutf,"w");
    if ( fp == NULL ) {
      std::cerr << "Cannot write to file " << soutf << std::endl;
      return (R_NilValue);
    }
    SEXP dname, rname, cname;
    int *dimX = INTEGER(coerceVector(getAttrib(matrix,R_DimSymbol),INTSXP));
    int nr = dimX[0], nc = dimX[1];
    int i,j;
    double* p_m;
    double e;
    PROTECT(dname = GET_DIMNAMES(matrix));
    PROTECT(rname = GET_ROWNAMES(dname));
    PROTECT(cname = GET_COLNAMES(dname));
    PROTECT(matrix = AS_NUMERIC(matrix));
    p_m = NUMERIC_POINTER(matrix);
    //fprintf(fp,"#MARKER");  // print # for headers
    fprintf(fp,"#CHROM\tBEGIN\tEND\tMARKER_ID");  // print # for headers
    for(i=0; i < nc; ++i) {
      fprintf(fp,"\t%s",CHAR(STRING_ELT(cname,i)));
    }
    fprintf(fp,"\n");

    std::string chrom, beg, end, name, anno;
    const char* prname;
    for(i=0; i < nr; ++i) {
      //fprintf(fp,"%s",CHAR(STRING_ELT(rname,i)));
      prname = CHAR(STRING_ELT(rname,i));
      if ( fVcf::parseMarkerID(prname, chrom, beg, end, name, anno) >= 2 )
	fprintf(fp,"%s\t%s\t%s\t%s",chrom.c_str(),beg.c_str(),end.c_str(),prname);
      else
	error("Error in parsing marker ID %s",prname);
      for(j=0; j < nc; ++j) {
	e = p_m[i+j*nr];
	if ( std::isnan(e) ) {
	  fprintf(fp,"\tNA");
	}
	else {
	  fprintf(fp,"\t%.5lg",e);
	}
      }
      fprintf(fp,"\n");
    }
    UNPROTECT(4);
    fprintf(stderr,"Sucessfully wrote ( %d * %d ) matrix\n",nr,nc);
    return (R_NilValue);
  }

  
  // write matrix to a file. same to write.table, but write a header in a different mannaer
  SEXP writeNamedMatrix(SEXP rowNames, SEXP rowVals, SEXP outname) {
    // create the output file
    const char* soutf = strAllocCopy(outname);
    FILE* fp = fopen(soutf,"w");
    if ( fp == NULL ) {
      std::cerr << "Cannot write to file " << soutf << std::endl;
      return (R_NilValue);
    }

    // get the dimensions of two input matrices
    SEXP dval, cval, dname, cname;
    int *dimV = INTEGER(coerceVector(getAttrib(rowVals,R_DimSymbol),INTSXP));
    int nvr = dimV[0], nvc = dimV[1];
    int *dimN = INTEGER(coerceVector(getAttrib(rowNames,R_DimSymbol),INTSXP));
    int nnr = dimN[0], nnc = dimN[1];
    if ( nvr != nnr ) {
      std::cerr << "writeMatrix() : input matices are incompatible (" << nnr << "," << nnc << ") vs (" << nvr << "," << nvc << ")" << std::endl;
      return (R_NilValue);
    }

    int i,j;
    double *p_v;
    double e;
    PROTECT(dval = GET_DIMNAMES(rowVals));
    PROTECT(cval = GET_ROWNAMES(dval));
    PROTECT(rowVals = AS_NUMERIC(rowVals));
    p_v = NUMERIC_POINTER(rowVals);

    PROTECT(dname = GET_DIMNAMES(rowNames));
    PROTECT(cname = GET_ROWNAMES(dname));
    PROTECT(rowNames = AS_CHARACTER(rowNames));

    for(i=0; i < nnc; ++i) {
      if ( i == 0 ) fprintf(fp,"#");
      else fprintf(fp,"\t");
      fprintf(fp,"%s",CHAR(STRING_ELT(cname,i)));
    }
    for(i=0; i < nvc; ++i) {
      fprintf(fp,"\t%s",CHAR(STRING_ELT(cval,i)));
    }
    fprintf(fp,"\n");

    for(i=0; i < nvr; ++i) {
      for(j=0; j < nnc; ++j) {
	if ( j > 0 ) fprintf(fp,"\t");
	fprintf(fp,"%s",CHAR(STRING_ELT(rowNames,i+j*nvr)));
      }
      for(j=0; j < nvc; ++j) {
	e = p_v[i+j*nvr];
	if ( std::isnan(e) ) {
	  fprintf(fp,"\tNA");
	}
	else {
	  fprintf(fp,"\t%.5lg",e);
	}
      }
      fprintf(fp,"\n");
    }
    UNPROTECT(6);
    fprintf(stderr,"Sucessfully written (%d * %d) string matrix and ( %d * %d ) double matrix\n",nvr,nnc,nvr,nvc);
    return (R_NilValue);
  }

  // read a VCF file
  // fname : VCF file name
  // region : region to retrieve [chr:start-end] or NULL
  // field : field to retrieve
  // pass : TRUE if only PASS SNPs are retreived, FALSE for all SNPs
  // ind : indices of individuals to retrieve
  // matchRule : substring match rule
  SEXP readVcf(SEXP fname, SEXP region, SEXP field, SEXP pass, SEXP ind, SEXP rule) {
    const char* sfname = strAllocCopy(fname);
    const char* sregion = (region == R_NilValue) ? NULL : strAllocCopy(region);
    const char* sfield = (field == R_NilValue) ? "GT" : strAllocCopy(field);

    bool bpass = LOGICAL_VALUE(pass) ? true : false;
    const char* srule = (rule == R_NilValue) ? NULL : strAllocCopy(rule);
    int *p_ind;
    int nInds = 0;
    if ( ind == R_NilValue ) {
      p_ind = NULL;
    }
    else {
      PROTECT(ind = AS_INTEGER(ind));
      nInds = LENGTH(ind);
      p_ind = INTEGER_POINTER(ind);
    }

    //fprintf(stderr,"readVcf(%s, %s, %s, %d) called\n",sfname, sregion, sfield, nInds);

    //fprintf(stderr,"field : %s\n",sfield);
    //fprintf(stderr,"RULE : %s\n",srule);
    //fVcf vcf(sfname, sregion, sfield, srule, bpass, false, p_ind, nInds); 
    fVcf vcf;
    std::vector<int> subcols(nInds);
    for(int i=0; i < nInds; ++i) { subcols[i] = p_ind[i]; }

    //fprintf(stderr,"*** before vcf.load()\n");
    
    vcf.load(sfname, sregion, sfield, srule, bpass, subcols);

    //fprintf(stderr,"*** before vcf.readMarkers()\n");

    vcf.readMarkers();

    //fprintf(stderr,"*** after vcf.readMarkers()\n");
    vcf.close();

    SEXP genos, rownames, dimnames;    
    SEXP *p_r; //, *p_c;
    PROTECT( genos = allocMatrix(REALSXP, vcf.nMarkers, vcf.nInds) );
    PROTECT( rownames = NEW_CHARACTER(vcf.nMarkers) );
    //PROTECT( colnames = NEW_CHARACTER(vcf.nInds) );
    PROTECT( dimnames = allocVector(VECSXP, 2));

    //fprintf(stderr,"%d markers, %d samples in VCF, %d with phenotypes and covariates\n",vcf.nMarkers,vcf.nInds,nInds);

    p_r = CHARACTER_POINTER(rownames);    
    //p_c = CHARACTER_POINTER(colnames);
    //for(int i=0; i < vcf.nInds; ++i) { 
    //  p_c[i] = mkChar(vcf.inds[p_ind[i]].c_str());
    //}

    float f;
    for(int i=0; i < vcf.nMarkers; ++i) {
      p_r[i] = mkChar(vcf.markers[i].c_str());
      for(int j=0; j < vcf.nInds; ++j) {
	f = vcf.genos[i*vcf.nInds+j];
	if ( std::isnan(f) ) REAL(genos)[i + vcf.nMarkers*j] = NA_REAL;
	else REAL(genos)[i + vcf.nMarkers*j] = (double)f;
      }
    }

    SET_VECTOR_ELT(dimnames, 0, rownames);
    //SET_VECTOR_ELT(dimnames, 1, colnames);
    setAttrib(genos, R_DimNamesSymbol, dimnames);

    if ( p_ind == NULL ) UNPROTECT(3);
    else UNPROTECT(4);

    //fprintf(stderr,"Read %d markers across %d individuals\n",vcf.nMarkers,vcf.nInds);
    return( genos );
  }

  // read a VCF file
  // fname : VCF file name
  // region : region to retrieve [chr:start-end] or NULL
  // field : field to retrieve
  // ind : indices of individuals to retrieve
  // markers : marker IDs to match
  SEXP readVcfGroup(SEXP fname, SEXP field, SEXP ind, SEXP markers, SEXP sepchr) {
    const char* sfname = strAllocCopy(fname);
    const char* sfield = (field == R_NilValue) ? "GT" : strAllocCopy(field);
    int *p_ind;
    int nInds = 0;
    char ** p_markers;
    int nMarkers = 0;
    bool bSepchr = LOGICAL_VALUE(sepchr) ? true : false;

    std::vector<int> subcols;
    if ( ind == R_NilValue ) {
      p_ind = NULL;
    }
    else {
      PROTECT(ind = AS_INTEGER(ind));
      nInds = LENGTH(ind);
      p_ind = INTEGER_POINTER(ind);
      for(int i=0; i < nInds; ++i) 
	subcols.push_back(p_ind[i]);
    }

    if ( markers == R_NilValue ) {
      p_markers = NULL;
    }
    else {
      PROTECT(markers = AS_CHARACTER(markers));
      nMarkers = LENGTH(markers);
      p_markers = (char**)malloc(sizeof(char*) * nMarkers);
      for(int i=0; i < nMarkers; ++i) {
	p_markers[i] = R_alloc(strlen(CHAR(STRING_ELT(markers,i))),sizeof(char));
	strcpy(p_markers[i], CHAR(STRING_ELT(markers,i)));
      }
    }

    fVcf vcf;
    vcf.load(sfname, NULL, sfield, NULL, false, subcols);
    vcf.readMarkerGroup((const char**)p_markers, nMarkers, 0, bSepchr);
    vcf.close();

    SEXP genos, rownames, dimnames;    
    SEXP *p_r; //, *p_c;
    PROTECT( genos = allocMatrix(REALSXP, vcf.nMarkers, vcf.nInds) );
    PROTECT( rownames = NEW_CHARACTER(vcf.nMarkers) );
    //PROTECT( colnames = NEW_CHARACTER(vcf.nInds) );
    PROTECT( dimnames = allocVector(VECSXP, 2));

    //fprintf(stderr,"%d markers, %d samples in VCF, %d with phenotypes and covariates\n",vcf.nMarkers,vcf.nInds,nInds);

    p_r = CHARACTER_POINTER(rownames);    
    //p_c = CHARACTER_POINTER(colnames);
    //for(int i=0; i < vcf.nInds; ++i) { 
    //  p_c[i] = mkChar(vcf.inds[p_ind[i]].c_str());
    //}

    float f;
    for(int i=0; i < vcf.nMarkers; ++i) {
      p_r[i] = mkChar(vcf.markers[i].c_str());
      for(int j=0; j < vcf.nInds; ++j) {
	f = vcf.genos[i*vcf.nInds+j];
	if ( std::isnan(f) ) REAL(genos)[i + vcf.nMarkers*j] = NA_REAL;
	else REAL(genos)[i + vcf.nMarkers*j] = (double)f;
      }
    }

    SET_VECTOR_ELT(dimnames, 0, rownames);
    //SET_VECTOR_ELT(dimnames, 1, colnames);
    setAttrib(genos, R_DimNamesSymbol, dimnames);

    if ( p_ind == NULL ) {
      if ( markers != R_NilValue ) { 
	UNPROTECT(4);
	free(p_markers); 
      }
      else {
	UNPROTECT(3);
      }
    }
    else {
      if ( markers != R_NilValue ) { 
	UNPROTECT(5);
	free(p_markers); 
      }
      else {
	UNPROTECT(4);
      }
    }

    //fprintf(stderr,"Read %d markers across %d individuals\n",vcf.nMarkers,vcf.nInds);
    return( genos );
  }

  SEXP readKinWithIDs(SEXP fname, SEXP indids) {
    //, SEXP field, SEXP ind, SEXP markers, SEXP sepchr) {
    //static bool readKinWithIDs(const char* inname, MatrixXd& K, int& flag, double& wsum, std::vector<std::string>& ids) {
    const char* sfname = strAllocCopy(fname);

    // read kinship values from the file
    pFile tf(sfname);
    char buf[N_MAGIC];
    tf.read(buf,N_MAGIC);
    int i, n, flag;
    double wsum;
    if ( strncmp(buf,MAGIC_KIN,N_MAGIC) != 0 ) error("MAGIC number for kinship matrix does not match");
    tf.read(&n,sizeof(int));
    tf.read(&flag,sizeof(int)); 
    tf.read(&wsum,sizeof(double)); 
    int size = n * ( n + 1 ) / 2;
    double* p = new double[size];
    int nr = tf.read(p, sizeof(double)*size);
    if ( nr != (int)sizeof(double)*size ) error("Kinship matrix size is truncated");
    fprintf(stderr,"n=%d\n",n);

    // read IDs anyways
    std::vector<std::string> ids;
    char s[65535];
    int sz;
    for(int i=0; i < n; ++i) {
      tf.read(&sz,sizeof(int));
      tf.read(s,sz);
      s[sz] = '\0';
      ids.push_back(s);
    }
    tf.close();
    
    // read individual IDs provided and produce a map
    std::vector<int> icols;
    int nv;
    if ( indids != R_NilValue ) {
      std::map<std::string,int> mInds;
      int im;

      PROTECT(indids = AS_CHARACTER(indids));
      int nInds = LENGTH(indids);
      for(i=0; i < nInds; ++i) {
	mInds[CHAR(STRING_ELT(indids,i))] = i;
      }
      UNPROTECT(1);

      // need to read the individual IDs and find the matching one
      nv = nInds;
      icols.resize(nv,-1);
      
      for(i=0; i < n; ++i) {
	if ( mInds.find(ids[i]) != mInds.end() ) { // found match
	  im = mInds[ids[i]];
	  if ( icols[im] >= 0 ) error("Duplicated ID %s was provided",ids[i].c_str());
	  else icols[im] = i;
	}
      }

      // make sure that every individual ID was included in the kinship
      for(i=0; i < nv; ++i) {
	if ( icols[i] < 0 ) {
	  for(std::map<std::string,int>::iterator it = mInds.begin();
	      it != mInds.end(); ++it) {
	    if ( it->second == i ) error("Individual ID %s cannot be found from the kinship file",it->first.c_str());
	  }
	  error("Something went wrong. Missing individual at index %d, but could not find in the dictionary",i);
	}
      }
    }
    else {
      icols.resize(n);
      for(i=0; i < n; ++i) icols[i] = i;
      nv = n;
    }

    SEXP genos, rownames, dimnames;
    SEXP *p_r;
    PROTECT( genos = allocMatrix(REALSXP, nv, nv) );
    PROTECT( rownames = NEW_CHARACTER(nv) );
    p_r = CHARACTER_POINTER(rownames);
    PROTECT( dimnames = allocVector(VECSXP, 2));

    double* p_genos = REAL(genos);
    int k;
    //double x;
    for(int i=0; i < nv; ++i) {
      for(int j=0; j <= i; ++j) {
	//tf.read(&x, sizeof(double));
	k = (icols[i] > icols[j]) ? (icols[i]*(icols[i]+1)/2+icols[j]) : (icols[j]*(icols[j]+1)/2+icols[i]);
	if ( i == j ) {
	  p_genos[i + nv*j] = p[k];
	}
	else {
	  p_genos[i + nv*j] = p_genos[j + nv*i] = p[k];
	}
      }

      p_r[i] = mkChar(ids[icols[i]].c_str());
    }

    delete [] p;

    SET_VECTOR_ELT(dimnames, 0, rownames);
    setAttrib(genos, R_DimNamesSymbol, dimnames);

    UNPROTECT(3);

    return(genos);
  }

  SEXP writeKinWithIDs(SEXP fname, SEXP kin) {
    const char* soutf = strAllocCopy(fname);

    FILE* fp = fopen(soutf,"wb");
    if ( fp == NULL ) {
      std::cerr << "Cannot write to file " << soutf << std::endl;
      return (R_NilValue);
    }

    SEXP dname, rname, cname;
    int *dimX = INTEGER(coerceVector(getAttrib(kin,R_DimSymbol),INTSXP));
    int nr = dimX[0], nc = dimX[1];
    int i,l;

    if ( nr != nc ) error("Row number %d is not identical to %d",nr,nc);

    PROTECT(dname = GET_DIMNAMES(kin));
    PROTECT(rname = GET_ROWNAMES(dname));
    PROTECT(cname = GET_COLNAMES(dname));
    PROTECT(kin = AS_NUMERIC(kin));
    double* p_kin = NUMERIC_POINTER(kin);

    double wsum = 1;
    int nzero = 0;
    fwrite(MAGIC_KIN, sizeof(char), N_MAGIC, fp);
    fwrite(&nr, sizeof(int), 1, fp);
    fwrite(&nzero, sizeof(int), 1, fp);
    fwrite(&wsum, sizeof(double), 1, fp);
    for(i=0; i < nr; ++i) {
      fwrite(&p_kin[i*nr], sizeof(double), (i+1), fp);
      //for(j=0; j <=i; ++j, k) {
      //fwrite(&p_kin[j+i*nr],sizeof(double));
      //}
    }
    const char* prname;
    for(i=0; i < nr; ++i) {
      prname = CHAR(STRING_ELT(rname,i));
      l = strlen(prname);
      fwrite(&l, sizeof(int), 1, fp);
      fwrite(prname, sizeof(char), l, fp);
    }

    UNPROTECT(4);
    fprintf(stderr,"Sucessfully wrote ( %d * %d ) matrix into EPACTS kinship format\n",nr,nc);
    return (R_NilValue);
  }

  SEXP readEigenWithIDs(SEXP fname) {
    const char* sfname = strAllocCopy(fname);

    // read kinship values from the file
    pFile tf(sfname);
    char buf[N_MAGIC];
    tf.read(buf,N_MAGIC);
    int i, nr, nc;
    double trK;
    if ( strncmp(buf,MAGIC_EIG,N_MAGIC) != 0 ) error("MAGIC number for eigen matrix does not match");

    tf.read(&nr,sizeof(int));
    tf.read(&nc,sizeof(int));
    tf.read(&trK,sizeof(double));

    SEXP evecs, evals, rownames, dimnames;

    PROTECT( evecs = allocMatrix(REALSXP, nr, nc) );
    PROTECT( evals = allocVector(REALSXP, nc) );

    double* p_evecs = REAL(evecs);
    double* p_evals = REAL(evals);
    
    tf.read(p_evals, sizeof(double) * nc);
    tf.read(p_evecs, sizeof(double) * nc * nr);

    // read IDs anyways
    std::vector<std::string> ids;
    char s[65535];
    int sz;
    for(int i=0; i < nr; ++i) {
      tf.read(&sz,sizeof(int));
      tf.read(s,sz);
      s[sz] = '\0';
      ids.push_back(s);
    }
    tf.close();

    PROTECT( rownames = NEW_CHARACTER(nr) );
    PROTECT( dimnames = allocVector(VECSXP,2) );

    SEXP* p_r = CHARACTER_POINTER(rownames);
    
    for(int i=0; i < nr; ++i) {
      p_r[i] = mkChar(ids[i].c_str());
    }

    SET_VECTOR_ELT(dimnames, 0, rownames );
    setAttrib( evecs, R_DimNamesSymbol, dimnames );

    SEXP result;
    PROTECT(result = allocVector(VECSXP, 3));


    SEXP tr;
    PROTECT(tr = NEW_NUMERIC(1));
    double* p_tr = NUMERIC_POINTER(tr);
    p_tr[0] = trK;

    SET_VECTOR_ELT(result, 0, evals);
    SET_VECTOR_ELT(result, 1, evecs);    
    SET_VECTOR_ELT(result, 2, tr);    

    SEXP header;
    PROTECT(header = allocVector(STRSXP,3));
    SET_STRING_ELT(header, 0, mkChar("evals"));
    SET_STRING_ELT(header, 1, mkChar("evecs"));
    SET_STRING_ELT(header, 2, mkChar("tr"));

    SET_NAMES(result,header);

    //setAttrib(result, 

    UNPROTECT(7);

    return(result);
  }

  SEXP writeEigenWithIDs(SEXP fname, SEXP evecs, SEXP evals, SEXP trK) {
    const char* soutf = strAllocCopy(fname);

    FILE* fp = fopen(soutf,"wb");
    if ( fp == NULL ) {
      std::cerr << "Cannot write to file " << soutf << std::endl;
      return (R_NilValue);
    }

    SEXP dname, rname, cname;
    int *dimX = INTEGER(coerceVector(getAttrib(evecs,R_DimSymbol),INTSXP));
    int nr = dimX[0], nc = dimX[1];

    PROTECT(dname = GET_DIMNAMES(evecs));
    PROTECT(rname = GET_ROWNAMES(dname));
    PROTECT(evecs = AS_NUMERIC(evecs));
    double* p_evecs = NUMERIC_POINTER(evecs);
    
    double tr = (trK == R_NilValue) ? 1. : NUMERIC_VALUE(trK);

    PROTECT( evals = AS_NUMERIC(evals) );
    int nv = LENGTH(evals);
    double* p_evals = NUMERIC_POINTER(evals);

    int i,l;

    if ( nc != nv ) error("Eigenvector matrix (%d x %d) is not compatiblw with eigenvalue size of %d", nr, nc, nv);

    fwrite(MAGIC_EIG, sizeof(char), N_MAGIC, fp);
    fwrite(&nr, sizeof(int), 1, fp);
    fwrite(&nc, sizeof(int), 1, fp);
    fwrite(&tr, sizeof(double), 1, fp);
    fwrite(p_evals, sizeof(double), nc, fp);
    fwrite(p_evecs, sizeof(double), nc * nr, fp);

    const char* prname;
    for(i=0; i < nr; ++i) {
      prname = CHAR(STRING_ELT(rname,i));
      l = strlen(prname);
      fwrite(&l, sizeof(int), 1, fp);
      fwrite(prname, sizeof(char), l, fp);
    }

    UNPROTECT(4);
    fprintf(stderr,"Sucessfully wrote ( %d * %d ) matrix into EPACTS Eigenvalue format\n",nr,nc);
    return (R_NilValue);
  }

//   SEXP readKinWithIDs(SEXP fname) {
//     //, SEXP field, SEXP ind, SEXP markers, SEXP sepchr) {
//     //static bool readKinWithIDs(const char* inname, MatrixXd& K, int& flag, double& wsum, std::vector<std::string>& ids) {
//     const char* sfname = strAllocCopy(fname);
//     pFile tf(sfname);
//     char buf[N_MAGIC];
//     tf.read(buf,N_MAGIC);
//     int n, flag;
//     double wsum;
//     if ( strncmp(buf,MAGIC_KIN,N_MAGIC) != 0 ) error("MAGIC number for kinship matrix does not match");
//     tf.read(&n,sizeof(int));
//     tf.read(&flag,sizeof(int)); 
//     tf.read(&wsum,sizeof(double)); 
//     int size = n * ( n + 1 ) / 2;
//     double* p = new double[size];
//     int nr = tf.read(p, sizeof(double)*size);
//     if ( nr != (int)sizeof(double)*size ) error("Kinship matrix size is truncated");

//     fprintf(stderr,"n=%d\n",n);

//     SEXP genos, rownames, dimnames;
//     SEXP *p_r;
//     PROTECT( genos = allocMatrix(REALSXP, n, n) );
//     PROTECT( rownames = NEW_CHARACTER(n) );
//     p_r = CHARACTER_POINTER(rownames);
//     PROTECT( dimnames = allocVector(VECSXP, 2));

//     double* p_genos = REAL(genos);
//     //double x;
//     for(int i=0, k=0; i < n; ++i) {
//       for(int j=0; j <= i; ++j, ++k) {
// 	//tf.read(&x, sizeof(double));
// 	if ( i == j ) {
// 	  p_genos[i + n*j] = p[k];
// 	}
// 	else {
// 	  p_genos[i + n*j] = p_genos[j + n*i] = p[k];
// 	}
//       }
//     }

//     delete [] p;

//     std::vector<std::string> ids;
//     char s[65535];
//     int sz;
//     for(int i=0; i < n; ++i) {
//       tf.read(&sz,sizeof(int));
//       tf.read(s,sz);
//       s[sz] = '\0';
//       p_r[i] = mkChar(s);
//     }
//     tf.close();

//     SET_VECTOR_ELT(dimnames, 0, rownames);
//     setAttrib(genos, R_DimNamesSymbol, dimnames);

//     UNPROTECT(3);

//     return(genos);
//   }
}
