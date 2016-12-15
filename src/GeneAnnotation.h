#ifndef _GENEANNOTATION_H_
#define _GENEANNOTATION_H_

#include "AnnotationResultCollection.h"
#include "AnnotationOutput.h"

struct GeneAnnotationParam{
  GeneAnnotationParam():
      upstreamRange(50),
      downstreamRange(50),
      spliceIntoExon(3),
      spliceIntoIntron(8) {};
  int upstreamRange;      // upstream def
  int downstreamRange;    // downstream def
  int spliceIntoExon;     // essential splice site def
  int spliceIntoIntron;   // essentail splice site def
};


class GeneAnnotation{
 public:
  GeneAnnotation():allowMixedVariation(false)
  {};
  virtual ~GeneAnnotation() {
  };
  /// @return 0: success
  int openGeneFile(const char* geneFileName, const char* geneFormat) {
    // set foramt
    GeneFormat format;
  
    std::string s = toLower(geneFormat);
    if (s == "refflat") {
      format.setRefFlatFormat();
      LOG << "Input gene format: refFlat\n";
    } else if (s == "knowngene") {
      format.setUCSCKnownGeneFormat();
      LOG << "Input gene format: knownGene\n";      
    } else if (s == "refgene") {
      format.setRefGeneFormat();
      LOG << "Input gene format: refGene\n";      
    } else {
      fprintf(stderr, "Unknown format (other than refFlat, knownGene, refGene)!\nNow quitting...\n");
      LOG << "Input gene format is wrong!\n";      
      abort();
    }

    // read gene file
    fprintf(stderr, "Load gene file %s...\n", geneFileName);
    std::string line;
    std::vector<std::string> fields;
    LineReader lr(geneFileName);
    int totalGene = 0;
    while (lr.readLine(&line) > 0) {
      stringStrip(&line);
      if ( (line.size()>0 && line[0] == '#') || line.size() == 0) continue; // skip headers and empty lines
      Gene g;
      g.readLine(line.c_str(), format);
      this->geneList[g.chr].push_back(g);
      totalGene ++;
    }
    // make sure genes are ordered
    this->sortGene();
    fprintf(stderr, "DONE: %d gene loaded.\n", totalGene);
    LOG << "Gene file " << geneFileName << " loads succeed!\n";
    return 0;
  };
  void openCodonFile(const char* codonFileName) {
    fprintf(stderr, "Load codon file %s...\n", codonFileName);
    this->codon.open(codonFileName);
    fprintf(stderr, "DONE: codon file loaded.\n");
    LOG << "Codon file " << codonFileName << " loads succeed!\n";
    return;
  };
  void openReferenceGenome(const char* referenceGenomeFileName) {
    fprintf(stderr, "Load reference genome %s...\n", referenceGenomeFileName);
    this->gs.open(referenceGenomeFileName);
    fprintf(stderr, "DONE: %d chromosomes and %ld bases are loaded.\n", this->gs.size(), this->gs.getGenomeLength());
    LOG << "Reference genome file " << referenceGenomeFileName << " loads succeed!\n";
    return;
  };
  void openPriorityFile(const char* fileName) {
    fprintf(stderr, "Load priority file %s...\n", fileName);
    int ret = this->priority.open(fileName);
    fprintf(stderr, "DONE: %d priority annotation types loaded.\n", ret);
    LOG << "Priority file " << fileName << " load succeed!\n";

    this->outputter.setPriority(this->priority);
    return;
  };

  void outputAnnotationStats(const char* outputFileName) {
    // output frequency files
    std::string fn = outputFileName;

    // output annotation frequency (all types of annotation)
    std::string ofs = fn+".anno.frq";
    this->printAnnotationFrequency(ofs.c_str());
    fprintf(stderr, "DONE: Generated frequency of each annotype type in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of each annotation type in " << ofs << " succeed!\n";

    // output annotation frequency
    ofs = fn+".top.anno.frq";
    this->printTopPriorityAnnotationFrequency(ofs.c_str());
    fprintf(stderr, "DONE: Generated frequency of each highest priority annotation type in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of high priority for highest priority annotation type in " << ofs << " succeed!\n";

    // output Ts/Tv ratio
    this->printTsTvRatio();

    // output base change frequency
    ofs = fn+".base.frq";
    this->printBaseChangeFrequency(ofs.c_str());
    fprintf(stderr, "DONE: Generated frequency of each base change in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of each base change in " << ofs << " succeed!\n";

    // output codon change frequency
    ofs = fn+".codon.frq";
    this->printCodonChangeFrequency(ofs.c_str());
    fprintf(stderr, "DONE: Generated frequency of each codon change in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of each codon change in " << ofs << " succeed!\n";

    // output indel length frequency
    ofs = fn+".indel.frq";
    this->printIndelLengthFrequency(ofs.c_str());
    fprintf(stderr, "DONE: Generated frequency of indel length in [ %s ].\n", ofs.c_str());
    LOG << "Generate frequency of indel length in " << ofs << " succeed!\n";
  };
public:
  // annotation will find all overlapping gene, and call annotateByGene for each gene
  // results will be stored in @param annotationResult
  // priority will also be calculated in annotationResult
  void annotate(const std::string& chrom,
                const int pos,
                const std::string& ref,
                const std::string& alt,
		const std::string& markerId ) {
    this->annotate(chrom, pos, ref, alt, markerId, &this->annotationResults);
    this->outputter.setAnnotationResult(this->annotationResults);
  };
  std::string getTopPriorityAnnotation() const{
    return this->outputter.getTopPriorityAnnotation();
  };
  std::string getFullAnnotation() const{
    return this->outputter.getFullAnnotation();
  };
private:
  // annotation will find all overlapping gene, and call annotateByGene for each gene
  // results will be stored in @param annotationResult
  // priority will also be calculated in annotationResult
  void annotate(const std::string& chrom,
                const int pos,
                const std::string& ref,
                const std::string& altParam,
		const std::string& markerId,
                AnnotationResultCollection* annotationResult) {
    // check VARATION_TYPE
    std::string alt = altParam;
    VARIATION_TYPE type = determineVariationType(ref, alt);
    AnnotationResult annotationPerGene;
    this->annotationResults.clear();
    
    if (type == MIXED) {
      // only annotate the first variation
      int commaPos = alt.find(',');
      alt = altParam.substr(0, commaPos);
      type = determineVariationType(ref, alt);
    }
    
    if (type == NO_VARIATION) {
      annotationPerGene.add(MONOMORPHIC);
      this->annotationResults.push_back(annotationPerGene);
    } else {
      // find near target genes
      std::vector<unsigned> potentialGeneIdx;
      this->findInRangeGene(chrom, pos, &potentialGeneIdx);
      // if Intergenic,  we will have (potentialGeneIdx.size() == 0)


      // annotate for each gene
      size_t i;
      for ( i = 0; i < potentialGeneIdx.size(); i++) {
        annotationPerGene.clear();
        this->annotateByGene(potentialGeneIdx[i], chrom, pos, ref, alt, type, markerId, &annotationPerGene);
        // annotationPerGene.dump();
        this->annotationResults.push_back(annotationPerGene);
      }
      if ( this->annotationResults.empty() ) {
        annotationPerGene.add(INTERGENIC);
        this->annotationResults.push_back(annotationPerGene);
      }
    }
    annotationResults.sortByPriority(this->priority);
    // record frquency info
    updateTypeFrequency(type, ref, alt);
    updateAnnotationFrequency(annotationResults);
    return;
  };
public:
  void updateTypeFrequency(const VARIATION_TYPE& type, const std::string& ref, const std::string& alt) {
    switch (type) {
      case SNP:
        this->baseFreq.add(ref + "->" +alt);
        break;
      case INS:
      case DEL:
        this->indelLengthFreq.add(calculateIndelLength(ref, alt));
      default:
        break;
    }
  }
  void updateAnnotationFrequency(const AnnotationResultCollection& result) {
    assert( !result.empty());
    if (result.size() == 1 ) {
      if (result.getTopAnnotation()[0].getType()[0] == INTERGENIC) {
        this->annotationTypeFreq.add(INTERGENIC);
        this->topPriorityAnnotationTypeFreq.add(INTERGENIC);
        return;
      }
      if (result.getTopAnnotation()[0].getType()[0] == MONOMORPHIC) {
        this->annotationTypeFreq.add(INTERGENIC);
        this->topPriorityAnnotationTypeFreq.add(INTERGENIC);
        return;
      }
    }
    
      const std::vector<AnnotationResult>& top = result.getTopAnnotation();
      this->topPriorityAnnotationTypeFreq.add(top[0].getType()[0]);
      
      const std::vector<AnnotationResult>& all = result.getAllAnnotation();
      for (size_t i = 0; i < all.size(); ++i ){
        size_t n = all[i].getType().size();
        for (size_t j = 0; j < n; j++ ){
          this->annotationTypeFreq.add( all[i].getType() [j] );
        }
      }

  };
  void setAnnotationParameter(GeneAnnotationParam& param) {
    this->param = param;
  };
  void printAnnotationFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    unsigned int n = this->annotationTypeFreq.size();
    for (unsigned int i = 0; i < n; i++){
      AnnotationType t;
      int freq;
      this->annotationTypeFreq.at(i, &t, &freq);
      fprintf(fp, "%s\t%d\n", AnnotationString[t], freq);
    }
    fclose(fp);
  };
  void printTopPriorityAnnotationFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    unsigned int n = this->topPriorityAnnotationTypeFreq.size();
    for (unsigned int i = 0; i < n; i++){
      AnnotationType t;
      int freq;
      this->topPriorityAnnotationTypeFreq.at(i, &t, &freq);
      fprintf(fp, "%s\t%d\n", AnnotationString[t], freq);
    }
    fclose(fp);
  };

  void printBaseChangeFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    const unsigned int n = this->baseFreq.size();
    for (unsigned int i = 0; i < n; i++){
      std::string k;
      int freq;
      this->baseFreq.at(i, &k, &freq);
      fprintf(fp, "%s\t%d\n", k.c_str(), freq);
    }
    fclose(fp);
  };
  
  void printCodonChangeFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    const unsigned int n = this->codonFreq.size();
    for (unsigned int i = 0; i < n; i++){
      std::string k;
      int freq;
      this->codonFreq.at(i, &k, &freq);
      fprintf(fp, "%s\t%d\n", k.c_str(), freq);
    }
    fclose(fp);
  };
  void printIndelLengthFrequency(const char* fileName){
    FILE* fp = fopen(fileName, "wt");
    assert(fp);
    const unsigned int n = this->indelLengthFreq.size();
    for (unsigned int i = 0; i < n; i++){
      int key;
      int freq;
      this->indelLengthFreq.at(i, &key, &freq);
      fprintf(fp, "%d\t%d\n", key, freq);
    }
    fclose(fp);
  };
  void printTsTvRatio(){
    const unsigned int n = this->baseFreq.size();
    unsigned int ts = 0;
    unsigned int tv = 0;
    for (unsigned int i = 0; i < n; i++){
      std::string k;
      int freq;
      this->baseFreq.at(i, &k, &freq);

      // summarize ts/tv
      if (k.size() == 4 ||  // only process A->T type of mutation
          k[1] == '-' ||
          k[2] == '>') {
        char c1 = k[0];
        char c2 = k[3];
        // flip G to A and T to C
        // then compare if flipped alleles are the same.
        if (c1 == 'G') c1 = 'A';
        if (c1 == 'T') c1 = 'C';
        if (c2 == 'G') c2 = 'A';
        if (c2 == 'T') c2 = 'C';
        if (c1 != 'A' && c1 != 'C') continue;
        if (c2 != 'A' && c2 != 'C') continue;
        if (c1 == c2) {
          ts += freq;
        } else {
          tv += freq;
        }
      }
    }
    // output ts/tv statistics
    if (tv != 0) {
      double tstv = 1.0 * ts / tv;
      LOG << "Ts/Tv ratio: " << tstv << "\n";
      fprintf(stderr, "Ts/Tv ratio: %g\n", tstv);
    } else {
      LOG << "Ts/Tv ratio: NA\n";
      fprintf(stderr, "Ts/Tv ratio: NA\n");
    }
    LOG << "Ts observed: " << ts << " times; Tv observed: " << tv << " times.\n";
    fprintf(stderr, "Ts observed: %d  times; Tv observed: %d times.\n", ts, tv);
  };
  
 private:
  // make sure genes are ordered
  void sortGene() {
    std::map<std::string, std::vector<Gene> >:: iterator it;
    for (it = this->geneList.begin(); it != this->geneList.end(); it ++){
      std::sort( it->second.begin(), it->second.end(), GeneCompareLess);
    }
  };
  // store results in @param potentialGeneIdx
  // find gene whose range plus downstream/upstream overlaps chr:pos
  void findInRangeGene(const std::string& chr, const int pos, std::vector<unsigned int>* potentialGeneIdx) {
    assert(potentialGeneIdx);
    potentialGeneIdx->clear();

    std::vector<Gene>& g = this->geneList[chr];
    unsigned int gLen = g.size();
    if (gLen == 0) {
      return;
    }
    int maxDist = (param.upstreamRange > param.downstreamRange) ? param.upstreamRange : param.downstreamRange;
    Range r ((pos - maxDist), (pos + maxDist));
    for (unsigned int i = 0; i < gLen; i++ ){
      if (g[i].tx.start <= r.start) {
        if (g[i].tx.end < r.start){
          continue;
        } else
          potentialGeneIdx->push_back(i);
      } else if (r.isInRange(g[i].tx.start)) {
        potentialGeneIdx->push_back(i);
      } else {
        break;
      }
    }
#if 0
    for (unsigned int i = 0 ; i < potentialGeneIdx->size() ; i++){
      printf("%d, ", (*potentialGeneIdx)[i]);
    }
    printf("\n");
#endif
    return;
  };
  /**
   * fill the actual base in @param refTriplet and @param altTriplet
   * we consider @param forwardStrand, so for forward strand, we copy from this->reference,
   * or, we copy the reverse complement from this->reference
   */
  void fillTriplet(const std::string& chr, const int variantPos, const int codonPos[3], bool forwardStrand,
                   const std::string& ref, const std::string& alt,
                   char refTriplet[3], char altTriplet[3]) {
    assert(ref.size() == 1 && alt.size() == 1);
    // deal with no-existing chromosomes (e.g. chrM is in input file, but not in reference file)
    if (!this->gs.exists(chr)) { 
      refTriplet[0] = refTriplet[1] = refTriplet[2] = 'N';
      altTriplet[0] = altTriplet[1] = altTriplet[2] = 'N';
      return;
    }
    const Chromosome& seq = this->gs[chr];
    if (codonPos[0] < 0 || codonPos[2] > seq.size()) {
      refTriplet[0] = refTriplet[1] = refTriplet[2] = 'N';
      altTriplet[0] = altTriplet[1] = altTriplet[2] = 'N';
    } else {
      refTriplet[0] = seq[codonPos[0] - 1];
      refTriplet[1] = seq[codonPos[1] - 1];
      refTriplet[2] = seq[codonPos[2] - 1];
      altTriplet[0] = (variantPos != codonPos[0]) ? seq[codonPos[0] - 1] : alt[0];
      altTriplet[1] = (variantPos != codonPos[1]) ? seq[codonPos[1] - 1] : alt[0];
      altTriplet[2] = (variantPos != codonPos[2]) ? seq[codonPos[2] - 1] : alt[0];
    }
  };
  AnnotationType determineSNVType(const std::string& refAAName, const std::string& altAAName, const int codonNum){
    if (refAAName == Codon::unknownAA || altAAName == Codon::unknownAA) {
      return SNV;
    } else if (Codon::isStopCodon(refAAName) && !Codon::isStopCodon(altAAName)) {
      return STOP_LOSS;
    } else if (!Codon::isStopCodon(refAAName) && Codon::isStopCodon(altAAName)) {
      return STOP_GAIN;
    } else if (refAAName == "Met" && altAAName != "Met" && codonNum <= 3) {
      return START_LOSS;
    } else if (refAAName != "Met" && altAAName == "Met" && codonNum <= 3) {
      return START_GAIN;
    } else if (refAAName == altAAName) {
      return SYNONYMOUS;
    } else {
      return NONSYNONYMOUS;
    }
  };
  /**
   * To annotate for insertion is very similar to annotate SNP, and the only difference is that
   * insertion in the exon could cause frameshift/codon_insertion/codon_deletion
   */
  void annotateIns(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt, const std::string& markerId, AnnotationResult* result) {
    Gene& g = this->geneList[chr][geneIdx];
    result->add(g);

    // might useful vars.
    int dist2Gene;
    int utrPos, utrLen;
    int exonNum; // which exon
    // int codonNum; // which codon
    // int codonPos[3] = {0, 0, 0}; // the codon position
    // AnnotationType type; // could be one of
    int intronNum; // which intron
    bool isEssentialSpliceSite;

    result->add(INSERTION);
    if (g.isUpstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(UPSTREAM);
    } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(DOWNSTREAM);
    } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
      result->add(EXON);
      if (g.isCoding()) {
        if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR5);
        } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR3);
        } else { // cds part has base change
          int insertSize = alt.size() - ref.size();
          if (insertSize % 3 == 0) {
            result->add(CODON_GAIN);
            std::string s;
            char triplet[3];
            std::string aaName;
            s+= WITHIN_GENE_LEFT_DELIM;
            for (unsigned int i = ref.size(); i < alt.size();){
              triplet[0 ] = alt[i++];
              triplet[1] = alt[i++];
              triplet[2] = alt[i++];
              if (!g.forwardStrand)
                reverseComplementTriplet(triplet);
              aaName = this->codon.toAA(triplet);
              s+= aaName;
            }
            s+= WITHIN_GENE_RIGHT_DELIM;
            result->addDetail(CODON_GAIN, s);
          } else {
            result->add(FRAME_SHIFT);
          }
        }
      } else {
        result->add(NONCODING);
      }
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else if (g.isIntron(variantPos, &intronNum)) {
      result->add(INTRON);
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else {
      //annotation.push_back("Intergenic");
    }
  } // end annotateIns(...)
  /**
   * Deletion may across various regions
   * we use std::set to store all regions it came across
   *
   */
  void annotateDel(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt, const std::string& markerId, AnnotationResult* result) {
    Gene& g = this->geneList[chr][geneIdx];
    std::set<AnnotationType> annotationSet;

    // preprocessing alt.
    std::string cleanedAlt;
    if (alt == ".") {
      cleanedAlt = "";
    } else {
      cleanedAlt = alt;
    }

    // might useful vars.
    int dist2Gene;
    int utrPos, utrLen;
    int exonNum; // which exon
    /* int codonNum; // which codon */
    /* int codonPos[3] = {0, 0, 0}; // the codon position */
    /* AnnotationType type; // could be one of */
    int intronNum; // which intron
    bool isEssentialSpliceSite;

    // calculate range of the deletion
    // some cases:
    // e.g. ref = AG cleanedAlt = A
    // e.g. ref = AG cleanedAlt = ""
    int delBeg = variantPos + cleanedAlt.size();  // delBeg: inclusive
    int delEnd = variantPos + ref.size(); // delEnd: exclusive

    std::string overlappedCdsBase; // bases in the cds (if any)
    annotationSet.insert(DELETION);
    for (int pos = delBeg; pos < delEnd; pos++) {
      if (g.isUpstream(pos, param.upstreamRange, &dist2Gene)) {
        annotationSet.insert(UPSTREAM);
      } else if (g.isDownstream(pos, param.upstreamRange, &dist2Gene)) {
        annotationSet.insert(DOWNSTREAM);
      } else if (g.isExon(pos, &exonNum)){//, &codonNum, codonPos)) {
        annotationSet.insert(EXON);
        if (!g.isNonCoding()) {
          if (g.is5PrimeUtr(pos, &utrPos, &utrLen)) {
            annotationSet.insert(UTR5);
          } else if (g.is3PrimeUtr(pos, &utrPos, &utrLen)) {
            annotationSet.insert(UTR3);
          } else { // cds part has base change
            overlappedCdsBase.push_back(ref[ pos - delBeg ]);
          }
        } else {
          annotationSet.insert(NONCODING);
        }
        // check splice site
        if (g.isSpliceSite(pos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
          if (isEssentialSpliceSite)
            annotationSet.insert(ESSENTIAL_SPLICE_SITE);
          else
            annotationSet.insert(NORMAL_SPLICE_SITE);
        }
      } else if (g.isIntron(pos, &intronNum)) {
        annotationSet.insert(INTRON);
        // check splice site
        if (g.isSpliceSite(pos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
          if (isEssentialSpliceSite)
            annotationSet.insert(ESSENTIAL_SPLICE_SITE);
          else
            annotationSet.insert(NORMAL_SPLICE_SITE);
        }
      } else {
        //annotation.push_back("Intergenic");
      }
    } // end for
    // check how many codon in cds are delete
    if (overlappedCdsBase.size() > 0) {
      if (overlappedCdsBase.size() % 3 == 0) {
        annotationSet.insert(CODON_LOSS);
      } else {
        annotationSet.insert(FRAME_SHIFT);
      }
    }
    // store all existing annotation
    result->add(g);
    for (std::set<AnnotationType>::const_iterator it = annotationSet.begin();
         it != annotationSet.end();
         it++) {
      result->add(*it);
    };
  }; // end annotateDel
  /**
   * SV is the most complex scenario. fully support this is an ongoing work.
   * We will just annotation the region in the rough scale.
   */
  void annotateSV(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt, const std::string& markerId, AnnotationResult* result) {
    Gene& g = this->geneList[chr][geneIdx];
    result->add(g);

    // might useful vars.
    int dist2Gene;
    int utrPos, utrLen;
    int exonNum; // which exon
    /* int codonNum; // which codon */
    /* int codonPos[3] = {0, 0, 0}; // the codon position */
    /* AnnotationType type; // could be one of */
    int intronNum; // which intron
    bool isEssentialSpliceSite;
    result->add(STRUCTURE_VARIATION);
    if (g.isUpstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(UPSTREAM);
    } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(DOWNSTREAM);
    } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
      result->add(EXON);
      if (!g.isNonCoding()) {
        if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR5);
        } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR3);
        } else { // cds part has base change
          result->add(CODON_REGION);
        }
      } else{
        result->add(NONCODING);
      }
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else if (g.isIntron(variantPos, &intronNum)) {
      result->add(INTRON);
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else {
      //annotation.push_back("Intergenic");
    }
  } // end annotateSV
  void annotateSNP(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt, const std::string& markerId, AnnotationResult* result) {
    Gene& g = this->geneList[chr][geneIdx];
    result->add(g);

    // might useful vars.
    int dist2Gene;
    int utrPos, utrLen;
    int exonNum; // which exon
    int codonNum; // which codon
    int codonPos[3] = {0, 0, 0}; // the codon position
    /* AnnotationType type; // could be one of */
    int intronNum; // which intron
    bool isEssentialSpliceSite;

    if (g.isUpstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(UPSTREAM);
    } else if (g.isDownstream(variantPos, param.upstreamRange, &dist2Gene)) {
      result->add(DOWNSTREAM);
    } else if (g.isExon(variantPos, &exonNum)){//, &codonNum, codonPos)) {
      result->add(EXON);
      if (g.isCoding()) {
        if (g.is5PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR5);
        } else if (g.is3PrimeUtr(variantPos, &utrPos, &utrLen)) {
          result->add(UTR3);
        } else { // cds part has base change
          if (g.calculateCodonPosition(variantPos, &codonNum, codonPos)) {
            char refTriplet[3];
            char altTriplet[3];
            std::string refAAName;
            std::string altAAName;
            std::string refLetterName;
            std::string altLetterName;
            AnnotationType annotationType;
            // when reference genome is provided
            if (this->gs.size() >0){
              this->fillTriplet(chr, variantPos, codonPos, g.forwardStrand, ref, alt, refTriplet, altTriplet);
              if (!g.forwardStrand){
                complementTriplet(refTriplet);
                complementTriplet(altTriplet);
              };
              refAAName = this->codon.toAA(refTriplet);
              altAAName = this->codon.toAA(altTriplet);
              refLetterName = this->codon.toLetter(refTriplet);
              altLetterName = this->codon.toLetter(altTriplet);
              annotationType = this->determineSNVType(refAAName, altAAName, codonNum);

              result->add(annotationType);
              std::string s;
              s += WITHIN_GENE_LEFT_DELIM;
              s += refTriplet[0];
              s += refTriplet[1];
              s += refTriplet[2];
              s += "/";
              s += refAAName;
              s += "/";
              s += refLetterName;
              s += "->";
              s += altTriplet[0];
              s += altTriplet[1];
              s += altTriplet[2];
              s += "/";
              s += altAAName;
              s += "/";
              s += altLetterName;
              s += WITHIN_GENE_SEPARATOR;
              // quick patch about codon number
              char buf[128];
              sprintf(buf, "Base%d/%d",
                      codonNum + 1, g.getCDSLength());
              s += buf;
              s += WITHIN_GENE_SEPARATOR;
              s += "Codon";
              s += toStr( codonNum / 3 + 1 );
              s += "/";
              s += toStr(g.getCDSLength() / 3);
              s += WITHIN_GENE_SEPARATOR;
              s += "Exon";
              s += toStr(exonNum + 1); // convert 0 indexed to 1 indexed
              s += "/";
              s += toStr( (int)( g.exon.size()));
              s += WITHIN_GENE_RIGHT_DELIM;
              result->addDetail(annotationType, s);
              // record frequency
              this->codonFreq.add(refAAName+"->"+altAAName);
            } else {
              result->add(SNV);
            }
          } else {
            result->add(SNV);
          }
        }
      } else{
        result->add(NONCODING);
      }
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else if (g.isIntron(variantPos, &intronNum)) {
      result->add(INTRON);
      // check splice site
      if (g.isSpliceSite(variantPos, param.spliceIntoExon, param.spliceIntoIntron, &isEssentialSpliceSite)){
        if (isEssentialSpliceSite)
          result->add(ESSENTIAL_SPLICE_SITE);
        else
          result->add(NORMAL_SPLICE_SITE);
      }
    } else {
      //annotation.push_back("Intergenic");
    }
  } // end annotateSNP
  /**
   * annotation results will be store in @param result
   */
  void annotateByGene(int geneIdx, const std::string& chr, const int& variantPos, const std::string& ref, const std::string& alt, 
                      const VARIATION_TYPE& type, const std::string& markerId,
                      AnnotationResult* result){
    result->clear();
    switch(type) {
      case SNP:
        this->annotateSNP(geneIdx, chr, variantPos, ref, alt, markerId, result);
        break;
      case INS:
        this->annotateIns(geneIdx, chr, variantPos, ref, alt, markerId, result);
        break;
      case DEL:
        this->annotateDel(geneIdx, chr, variantPos, ref, alt, markerId, result);
        break;
      case SV:
        this->annotateSV(geneIdx, chr, variantPos, ref, alt, markerId, result);
        break;
      case MIXED:
      case UNKNOWN:
      default:
        LOG << "Currently we don't support this variation type: " << type << "\n";
        break;
    };
    
    return;
  };
  /**
   * @return indel length, for insertion, return positive number; or return negative number
   */
  int calculateIndelLength(const std::string& ref, const std::string& alt){
    int refLen = ref.size();
    int altLen = alt.size();
    if (alt == "." || alt == "<DEL>") {
      altLen = 0 ;
    }
    return (altLen - refLen);
  };
  /**
   * @return the variation type depending on the first entry in the alt field
   */
  VARIATION_TYPE determineVariationType(const std::string& ref, const std::string& alt) {
    if (alt == ".") {
      return NO_VARIATION;
    }
    if (alt.find(',') != std::string::npos) {
      return MIXED;
    }
    // NOTE: a single "." is used for deletion; but "G." can represent uncertain breakpoint
    //       we will check the former case first.
    if (alt == ".") {
      return DEL;
    }

    const char* ALLOWED_BASE = "ACGT";
    unsigned int refLen = ref.size();
    unsigned int altLen = alt.size();
    if (alt.find_first_not_of(ALLOWED_BASE) != std::string::npos) {
      // NOTE: SV usually contain "[" or "]" for rearrangment
      //       ">", "<" for haplotypes or large deletion/insertion.
      return SV;
    }
    if (refLen == altLen) {
      if (refLen == 1){
        return SNP;
      } else {
        return UNKNOWN;
      }
    } else if (refLen > altLen) {
      return DEL;
    } else if (refLen < altLen) {
      return INS;
    }
    return UNKNOWN;
  };
 private:
  // internal data
  std::map <std::string, std::vector<Gene> > geneList;   // chrom -> genes
  std::string annotation;

  
  // parameters 
  GeneAnnotationParam param;
  GenomeSequence gs;
  Codon codon;
  Priority priority;

  bool allowMixedVariation;       // VCF ALT field may have more than one variation e..g A,C

  // output related variables
  AnnotationResultCollection annotationResults;
  AnnotationOutput outputter;               // control output format

  // frequency related variables
  FreqTable<AnnotationType> annotationTypeFreq;                 // base change frequency
  FreqTable<AnnotationType> topPriorityAnnotationTypeFreq;      // base change frequency of top priority
  FreqTable<std::string> baseFreq;                              // base change frequency
  FreqTable<std::string> codonFreq;                             // codon change frequency
  FreqTable<int> indelLengthFreq;                               // for insertion, the value is positive; for deletion, positive
}; // end class GeneAnnotation



#endif /* _GENEANNOTATION_H_ */
