#ifndef _GENEFORMAT_H_
#define _GENEFORMAT_H_

/**
 * The following format all use 0-based index
 * 
 */
class GeneFormat{
public:
  /// WASH7P	NR_024540	chr1	-	14361	29370	29370	29370	11	14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,	14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,
  void setRefFlatFormat(){
    this->geneNameCol.clear();
    this->geneNameCol.push_back(0);
    this->transcriptNameCol.clear();
    this->transcriptNameCol.push_back(1);
    this->chrCol = 2;
    this->strandCol = 3;
    this->txStartCol = 4;
    this->txEndCol = 5;
    this->cdsStartCol = 6;
    this->cdsEndCol = 7;
    this->exonNumCol = 8;
    this->exonStartCol = 9;
    this->exonEndCol = 10;

    this->minimumExpectedColumn = 11;
  };
  /**
   * Can be 11 column or 12 column, however, the last column is usually the same as first one, so we will ignore it.
   uc001aaa.3      chr1    +       11873   14409   11873   11873   3       11873,12612,13220,      12227,12721,14409,
           uc001aaa.3
   uc010nxq.1      chr1    +       11873   14409   12189   13639   3       11873,12594,13402,      12227,12721,14409,      B7ZGX9    uc010nxq.1

   *
   * Schema: http://genome.csdb.cn/cgi-bin/hgTables?db=rn4&hgta_group=genes&hgta_track=knownGene&hgta_table=knownGene&hgta_doSchema=describe+table+schema
   */
  void setUCSCKnownGeneFormat(){
    this->geneNameCol.clear();
    this->geneNameCol.push_back(10);
    this->transcriptNameCol.clear();
    this->transcriptNameCol.push_back(0);
    
    this->chrCol = 1;
    this->strandCol = 2;
    this->txStartCol = 3;
    this->txEndCol = 4;
    this->cdsStartCol = 5;
    this->cdsEndCol = 6;
    this->exonNumCol = 7;
    this->exonStartCol = 8;
    this->exonEndCol = 9;

    this->minimumExpectedColumn = 11;
  };

  
  void setRefGeneFormat() {
    this->geneNameCol.clear();
    this->geneNameCol.push_back(12);
    this->transcriptNameCol.clear();
    this->transcriptNameCol.push_back(1);
    
    this->chrCol = 2;
    this->strandCol = 3;
    this->txStartCol = 4;
    this->txEndCol = 5;
    this->cdsStartCol = 6;
    this->cdsEndCol = 7;
    this->exonNumCol = 8;
    this->exonStartCol = 9;
    this->exonEndCol = 10;

    this->minimumExpectedColumn = 16;
  }
  int getMinimumExpectedColumn() const {
    return this->minimumExpectedColumn;
  };

public:
  std::vector<int> geneNameCol;   // gene name(s), join by '/'
  std::vector<int> transcriptNameCol;  // transcripts name(s), join by '/'
  int chrCol;
  int strandCol;
  int txStartCol;
  int txEndCol;
  int cdsStartCol;
  int cdsEndCol;
  int exonNumCol;
  int exonStartCol;
  int exonEndCol;
  int minimumExpectedColumn;
  static const char NAME_SEP = '/';
};

#endif /* _GENEFORMAT_H_ */
