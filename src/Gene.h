#ifndef _GENE_H_
#define _GENE_H_

#include "GeneFormat.h"
#include "Range.h"
#include "StringUtil.h"

class Gene{
public:
  // read refFlat file into internal data structure
  // chr1 will be chopped to 1.
  void readLine(const char* line, const GeneFormat& format) {
    std::vector< std::string > field;
    std::vector< std::string > exon_beg;
    std::vector< std::string > exon_end;
    int nf = stringTokenize(line, "\t", &field);
    if (nf < format.getMinimumExpectedColumn()) {
      static int nTimeError = 0;
      fprintf(stderr, "Unable to read this gene from: %s\n", line);
      if (nTimeError++ > 10) {
        fprintf(stderr, "Too many errors, now quiting...\n");
        exit(1);
      }
      return;
    }
    // read name
    this->geneName.clear();
    for (unsigned int i = 0; i < format.geneNameCol.size(); i++) {
      if (i) {
        this->geneName += "/";
      }
      this->geneName += field[format.geneNameCol[i]];
    }
    this->transcriptName.clear();
    for (unsigned int i = 0; i < format.transcriptNameCol.size(); i++) {
      if (i) {
        this->transcriptName += "/";
      }
      this->transcriptName += field[format.transcriptNameCol[i]];
    }
    // read others
    this->chr = chopChr(field[format.chrCol]);
    this->forwardStrand = (field[format.strandCol] == "+" ? true: false);
    if (!str2int(field[format.txStartCol], &this->tx.start)) {
      fprintf(stderr, "Gene file format error!\n");
      exit(1);
    }
    this->tx.start ++;
    this->tx.end = toInt(field[format.txEndCol]);
    int cdsStart = toInt(field[format.cdsStartCol]) + 1;
    int cdsEnd = toInt(field[format.cdsEndCol]);
    if (this->length(cdsStart, cdsEnd) == 0) {
      this->isNonCodingGene = true;
    } else {
      this->isNonCodingGene = false;
    }
    unsigned int nExon = toInt(field[format.exonNumCol]);
    stringTokenize(field[format.exonStartCol], ',', &exon_beg);
    stringTokenize(field[format.exonEndCol], ',', &exon_end);
    for (unsigned int i = 0; i < nExon; i++ ){
      this->exon.push_back(Range(toInt(exon_beg[i]) + 1, toInt(exon_end[i])));
    }
    if (!this->isNonCodingGene) {
      // we will load utr5 and utr3, at the last step, we will
      // swap them for reverse strand
      // load till left of cdsBegin
      unsigned i = 0;
      for (i = 0; i < nExon; i++ ){
        int beg = exon[i].start;
        int end = exon[i].end;
        if (this->isInRange(cdsStart, beg, end)) {
          if (beg != cdsStart) // avoid add empty range
            this->utr5.push_back(Range(beg, cdsStart - 1));
          break;
        } else {
          this->utr5.push_back(Range(beg, end));
        }
      }
      // load cds region
      for (; i < nExon; i++) {
        int beg = exon[i].start;
        int end = exon[i].end;
        if (this->isInRange(cdsStart, beg, end)) {
          beg = cdsStart;
        }
        if (this->isInRange(cdsEnd, beg, end)) {
          end = cdsEnd;
          this->cds.push_back(Range(beg, cdsEnd));
          break;
        } else {
          this->cds.push_back(Range(beg,end));
        }
      }
      // load from cdsEnd to the end of exon
      for (; i < nExon; i++) {
        int beg = exon[i].start;
        int end = exon[i].end;
        if (this->isInRange(cdsEnd, beg, end)) {
          if (cdsEnd != end) // avoid add empty range
            this->utr3.push_back(Range(cdsEnd + 1, end));
        } else {
          this->utr3.push_back(Range(beg, end));
        }
      }
      if (i != nExon) {
        assert (i == nExon);
      }
      if (!this->forwardStrand) {
        std::swap(this->utr5, this->utr3);
      }
    }

#if 0
    // debug code
    if (name == "DDX53") {
      assert ( 0 == getCDSLength() % 3 );
    }
#endif
#if 0
    // just for my curiosity
    if (!isNonCoding()) {
      if (name == "DDX53") {
        printf("\n");
      }
      printf("%s (%d) with %d exon has 5'utr (%d), 3'utr(%d), cds(%d), cds module 3(%d)\n",
             name.c_str(),
             getGeneLength(),
             (int)exon.size(),
             get5PrimeUTRLength(),
             get3PrimeUTRLength(), getCDSLength(),
             getCDSLength() % 3 );
    }
#endif
  };
  /**
   *@return true if @param pos is in upstream and return how far it is from the beginning of the gene
   */
  bool isUpstream(const int pos, const int upstreamRange, int* dist) {
    if (this->forwardStrand) {
      if (this->tx.start - upstreamRange < pos && pos < this->tx.start) {
        *dist = this->tx.start - pos;
        return true;
      }
    } else {
      if (this->tx.end < pos && pos < this->tx.end + upstreamRange) {
        *dist = this->tx.end - pos;
        return true;
      }
    }
    return false;
  };
  bool isDownstream(const int pos, const int downstreamRange, int* dist) {
    if (this->forwardStrand) {
      if (this->tx.end < pos && pos < this->tx.end + downstreamRange) {
        *dist = pos - this->tx.end;
        return true;
      }
    } else {
      if (this->tx.start - downstreamRange < pos && pos < this->tx.start) {
        *dist = this->tx.start - pos;
        return true;
      }
    }
    return false;
  };
  /**
   * @return true is @param variousPos is i 5'-UTR region,
   * //@param utrPos will store the relative position of @param variousPos to the leftmost position of 5' UTR
   * //@param utrLen will store the length of the 5'-UTR region
   */
  bool is5PrimeUtr(const int variantPos, int* utrPos, int* utrLen) {
    if (this->isNonCoding()) return false;
    if (this->isInRange(variantPos, this->utr5)){
      return true;
    };
    return false;
  };
  bool is3PrimeUtr(const int variantPos, int* utrPos, int* utrLen) {
    if (this->isNonCoding()) return false;
    if (this->isInRange(variantPos, this->utr3)){
      return true;
    };
    return false;
  };
  bool isExon(const int variantPos, int* exonNum){
    bool ret = this->isInRange(variantPos, this->exon, exonNum);
    if (ret && !this->forwardStrand){
      *exonNum = this->exon.size() - 1 - *exonNum; // e.g. exonNum = 0 and exon.size() = 5, then we should get exonNum = 4
    }
    return ret;
    /*
      if (isNonCoding()) {
      if (this->isInRange(variantPos, this->exon, exonNum))
      return true;
      } else {
      if (this->isInRange(variantPos, this->utr5) ||
      this->isInRange(variantPos, this->cds) ||
      this->isInRange(variantPos, this->utr3))
      return true;
      }
      return false;
    */
  };
  bool isCodingRegion(const int variantPos, int* codonNum){
    if (isNonCoding())
      return false;
    if (this->isInRange(variantPos, this->cds)){
      return true;
    }
    return false;
  };
  /**
   * @return the position of the codon right next to @param currentPos,
   * @param cdsIdx indicates which codon (1, 2, 3..) is the @param currentPos
   * @param offset indicate the calculated position is on the left or right
   */
  int nextCodonPos(const int currentPos, int* cdsIdx, const int offset) {
    if (currentPos < 0) return -1;
    assert (offset == 1 || offset == -1);
    int nextPos = -1;
    if (offset == 1) {
      nextPos = currentPos + 1;
      if (!this->isInRange(nextPos, this->cds[*cdsIdx])) {
        (*cdsIdx) ++;
        if (*cdsIdx >= (int)this->cds.size()) {
          return -1;
        }
        nextPos = this->cds[*cdsIdx].start;
      }
    } else {
      nextPos = currentPos - 1;
      if (!this->isInRange(nextPos, this->cds[*cdsIdx])) {
        (*cdsIdx) --;
        if (*cdsIdx < 0) {
          return -1;
        }
        nextPos = this->cds[*cdsIdx].end;
      }
    }
    return nextPos;
  };
  /**
   * @return true: if codonPos[3] are all valid position
   * @param codonNum : which base (inclusive, 1-based) has mutation. possible values: 1, 2, 3 ...
   *        e.g. codonNum = 5 meaning if we concatenate all codon together
   *        the mutation is in the 5th position, the 2rd codon.
   */
  bool calculateCodonPosition(const int variantPos, int* codonNum, int codonPos[3]){
    *codonNum = 0;
    if (this->forwardStrand) {
      unsigned int i;
      for (i = 0; i < this->cds.size() ; i++) {
        if (this->isInRange(variantPos, this->cds[i])){
          *codonNum += variantPos - this->cds[i].start + 1;
          break;
        } else {
          *codonNum += this->cds[i].length();
        }
      }
      int n = (*codonNum) % 3;
      int cdsIdx = i;
      switch(n){
        case 0:
          codonPos[2] = variantPos;
          codonPos[1] = nextCodonPos(codonPos[2], &cdsIdx, -1);
          codonPos[0] = nextCodonPos(codonPos[1], &cdsIdx, -1);
          break;
        case 1:
          codonPos[0] = variantPos;
          codonPos[1] = nextCodonPos(codonPos[0], &cdsIdx, 1);
          codonPos[2] = nextCodonPos(codonPos[1], &cdsIdx, 1);
          break;
        case 2:
          codonPos[1] = variantPos;
          codonPos[2] = nextCodonPos(codonPos[1], &cdsIdx, 1);
          cdsIdx = i;
          codonPos[0] = nextCodonPos(codonPos[1], &cdsIdx, -1);
          break;
      }
    } else { // backward
      int i;
      for (i = this->cds.size() - 1; i >= 0 ; i--) {
        if (this->isInRange(variantPos, this->cds[i])){
          *codonNum += this->cds[i].end - variantPos + 1;
          break;
        } else {
          *codonNum += this->cds[i].length();
        }
      }
      int n = (*codonNum) % 3;
      int cdsIdx = i;
      switch(n){
        case 0:
          codonPos[2] = variantPos;
          codonPos[1] = nextCodonPos(codonPos[2], &cdsIdx, +1);
          codonPos[0] = nextCodonPos(codonPos[1], &cdsIdx, +1);
          break;
        case 1:
          codonPos[0] = variantPos;
          codonPos[1] = nextCodonPos(codonPos[0], &cdsIdx, -1);
          codonPos[2] = nextCodonPos(codonPos[1], &cdsIdx, -1);
          break;
        case 2:
          codonPos[1] = variantPos;
          codonPos[2] = nextCodonPos(codonPos[1], &cdsIdx, -1);
          cdsIdx = i;
          codonPos[0] = nextCodonPos(codonPos[1], &cdsIdx, +1);
          break;
      }
    }
    if (codonPos[0] != -1 && codonPos[1] != -1 && codonPos[2] != -1)
      return true;
    return false;
  };
  bool isIntron(const int variantPos, int* intronNum){
    // strand is not an issue here
    for (unsigned int i = 1; i < this->exon.size(); i++) {
      if (this->exon[i-1].end < variantPos && variantPos < this->exon[i].start) {
        return true;
      }
    }
    return false;
  };
  bool isSpliceSite(const int variantPos, int spliceIntoExon, int spliceIntoIntron, bool* isEssentialSpliceSite){
    *isEssentialSpliceSite = false;
    unsigned int exonNumber = this->exon.size();
    // first check splice into exon
    if (this->isInRange(variantPos, this->exon[0].end - (spliceIntoExon - 1), this->exon[0].end)) {
      return true;
    }
    if (this->isInRange(variantPos, this->exon[exonNumber - 1].start, this->exon[exonNumber - 1].start + (spliceIntoExon - 1))){
      return true;
    }
    for (unsigned int i = 1; i < exonNumber - 1; i ++) {
      if (this->isInRange(variantPos, this->exon[i].start, this->exon[i].start + (spliceIntoExon - 1)))
        return true;
      if (this->isInRange(variantPos, this->exon[i].end - (spliceIntoExon - 1), this->exon[i].end))
        return true;
    }
    // check splice into intron (also mark isEssentialSpliceSite)
    // we define essential splice site is (GU...AG) in the intron, and next to exon
    // and GU, AG both have length 2.
    for (unsigned int i = 0; i < exonNumber - 1; i++ ) {
      if (this->isInRange(variantPos, this->exon[i].end + 1, this->exon[i].end + 2)) {
        *isEssentialSpliceSite = true;
        return true;
      } else if (this->isInRange(variantPos, this->exon[i+1].start - 2, this->exon[i+1].start - 1)) {
        *isEssentialSpliceSite = true;
        return true;
      }
      if (this->isInRange(variantPos, this->exon[i].end+1, this->exon[i].end + 1 + (spliceIntoIntron - 1))) {
        return true;
      } else if (this->isInRange(variantPos, this->exon[i+1].start - 1 - (spliceIntoIntron - 1), this->exon[i+1].start - 1)) {
        return true;
      }
    }
    return false;
  };
  size_t getExonNumber() const {
    return this->exon.size();
  }
  int getTotalLength(const std::vector<Range>& v) {
    int l = 0;
    for (unsigned int i = 0; i < v.size() ; i++ )
      l += v[i].length();
    return l;
  };
  int getExonLength() {
    return this->getTotalLength(this->exon);
  };
  int getCDSLength() {
    return this->getTotalLength(this->cds);
  };
  int get5PrimeUTRLength() {
    return this->getTotalLength(this->utr5);
  };
  int get3PrimeUTRLength() {
    return this->getTotalLength(this->utr3);
  };
  int getGeneLength() {
    return this->tx.length();
  };
  bool isNonCoding() {
    return this->isNonCodingGene;
  };
  bool isCoding() {
    return !this->isNonCodingGene;
  };
  /**
   * @return true if @param pos is in the range [@param beg, @param end] (inclusive on the boundaries).
   */
  bool isInRange(const int pos, const int beg, const int end) {
    if (beg > end) {
      fprintf(stderr, "in isInRange beg(%d) > end(%d).\n", beg, end);
    }
    if (beg <= pos && pos <= end)
      return true;
    return false;
  };
  bool isInRange(const int pos, const Range& r) {
    return (this->isInRange(pos, r.start, r.end));
  };
  bool isInRange(const int pos, const std::vector<Range>& r) {
    int idx;
    return this->isInRange(pos, r, &idx);
  };
  /**
   * Check if @param pos is in the range (@param r), and put the index to @param whichRange
   * so that @param pos is within @param r [ @param whichRange]
   */
  bool isInRange(const int pos, const std::vector<Range>& r, int* whichRange) {
    assert(whichRange);
    for (unsigned int i = 0; i < r.size(); i++) {
      if (this->isInRange(pos, r[i])){
        *whichRange = i;
        return true;
      }
    }
    *whichRange = -1;
    return false;
  };
  /**
   * @return the total length from @param beg to @param end, inclusive on the boundaries
   */
  int length(int beg, int end) {
    if (beg > end+1) {
      fprintf(stderr, "In length beg(%d) > end(%d) + 1; please check gene file format.\n", beg, end);
    }
    return (end - beg + 1);
  };
public:
  std::string geneName;
  std::string transcriptName;
  std::string chr;
  bool forwardStrand;
  Range tx;
  // used for nonCoding gene
  std::vector<Range> exon;
  // used for coding gene
  std::vector<Range> cds;
  std::vector<Range> utr5;
  std::vector<Range> utr3;
  bool isNonCodingGene;
  GeneFormat format;
};
bool GeneCompareLess(const Gene& a, const Gene& b) {
  if (a.tx.start != b.tx.start)
    return a.tx.start < b.tx.start;
  else
    return a.exon.size() < b.exon.size();
};
#endif /* _GENE_H_ */
