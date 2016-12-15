#ifndef _ANNOTATIONSTRING_H_
#define _ANNOTATIONSTRING_H_

/**
 * A class that provide strings corresponding to above annotation types.
 */
class OutputAnnotationString{
 public:
  OutputAnnotationString(){
    this->setFormat("default");
  }
  // support vaious output format, e.g. default, or epacts
  void setFormat(const char* format) {
    std::string f = format;
    f = toLower(f);
    if ( f == "default" ) {
      this->annotationString = OutputAnnotationString::defaultAnnotationString;
    } else if (f == "epact") {
      this->annotationString = OutputAnnotationString::epactAnnotationString;
    } else {
      fprintf(stderr, "Cannot recoginized format: [ %s ]!\n", format);
    };;
  };
  const char* operator [] (const int idx) {
    return this->annotationString[idx];
  };
 private:
  const char** annotationString;
  static const char* defaultAnnotationString[];
  static const char* epactAnnotationString[];
}; // end class OutputAnnotationString

const char* OutputAnnotationString::defaultAnnotationString[] = {
  "StructuralVariation",
  "Stop_Gain",
  "Stop_Loss",
  "Start_Gain",
  "Start_Loss",
  "Frameshift",
  "CodonGain",
  "CodonLoss",
  "CodonRegion",
  "Insertion",
  "Deletion",
  "Nonsynonymous",
  "Synonymous",
  "Essential_Splice_Site",
  "Normal_Splice_Site",
  "Utr5",
  "Utr3",
  "Exon",
  "Intron",
  "Upstream",
  "Downstream",
  "SNV",
  "Noncoding",
  "Intergenic",
  "Monomorphic"
};

const char* OutputAnnotationString::epactAnnotationString[]= {
  "StructuralVariation",
  "Nonsense",                 // diff
  "Stop_Loss",
  "Start_Gain",
  "Start_Loss",
  "Frameshift",
  "CodonGain",
  "CodonLoss",
  "CodonRegion",
  "Insertion",
  "Deletion",
  "Missense",                 // diff
  "Silent",                   // diff
  "Essential_Splice_Site",
  "Normal_Splice_Site",
  "Utr5",
  "Utr3",
  "Exon",
  "Intron",
  "Upstream",
  "Downstream",
  "SNV",
  "Noncoding",
  "Intergenic"
};


#endif /* _ANNOTATIONSTRING_H_ */
