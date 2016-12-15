#ifndef _ANNOTATIONOUTPUT_H_
#define _ANNOTATIONOUTPUT_H_

#include "Common.h"

/**
 * AnnotationOutput links the class AnnotationResult to a string representation
 * It contain a string templates.
 */
class AnnotationOutput{
 public:
  AnnotationOutput():
      topPriorityTemplate("$(TOP_TYPE):$[$(ALL_TOP_GENE) |]"),
      geneTemplate("$(GENE_NAME):$(GENE_STRAND):$[$(TYPE) :]"),
      fullTemplate("$[$(GENE_ANNOTATION) |]"),
      annotationResult(NULL),
      priority(NULL)
  {
  };
  void setAnnotationResult(const AnnotationResultCollection& r) {
    this->annotationResult = &r;
    this->buildKeywordDict();
  };
  void setPriority(const Priority& p) {
    this->priority = &p;
  };
  void buildKeywordDict() {
    if (this->annotationResult->empty()) {
      // something should not be happen, at least one type of anntation should be here.
      fprintf(stderr, "Internal logic error: no annotation found. \n");
      return;
    }
    if (this->annotationResult->size() == 1) { // intergenic or monomorhpic
      if ((*this->annotationResult)[0].getType()[0] == MONOMORPHIC ||
          (*this->annotationResult)[0].getType()[0] == INTERGENIC) {
        // will handle later in each separate functions.
        return;
      }
    };
    const std::vector<AnnotationResult>& top = this->annotationResult->getTopAnnotation();
    const std::vector<AnnotationResult>& all = this->annotationResult->getAllAnnotation();

    std::vector<std::string> res;
    for (size_t i= 0; i < top.size() ; ++i) {
      res.push_back(top[i].getGeneName());
    }

    // only keep unique gene names
    inplace_make_set(&res);

    topPriorityTemplate.add("TOP_GENE", res.size() ? res[0]: AnnotationString[INTERGENIC]);
    topPriorityTemplate.add("ALL_TOP_GENE", res);
    topPriorityTemplate.add("TOP_TYPE", AnnotationString[top[0].getType()[0]]);

    std::vector<std::string> geneTemplate;
    for (unsigned int i = 0; i != all.size(); i++){
      geneTemplate.push_back(getGeneAnnotation(all[i]));
    }
    fullTemplate.add("GENE_ANNOTATION", geneTemplate);
  };
  // Format:
  //  Most_priority:gene1|gene2
  std::string getTopPriorityAnnotation() const{
    if (this->annotationResult->size() == 1) { // intergenic or monomorhpic
      if ((*this->annotationResult)[0].getType()[0] == MONOMORPHIC) {
        return AnnotationString[MONOMORPHIC];
      }
      if ((*this->annotationResult)[0].getType()[0] == INTERGENIC) {
        return AnnotationString[INTERGENIC];
      }
    }
        
    std::string s;
    if (this->topPriorityTemplate.translate(&s)) {
      fprintf(stderr, "topPriorityTemplate failed translation!\n");
    }
    return s;
  };

  std::string getFullAnnotation() const{
    if (this->annotationResult->size() == 1) { // intergenic or monomorhpic
      if ((*this->annotationResult)[0].getType()[0] == MONOMORPHIC) {
        return AnnotationString[MONOMORPHIC];
      }
      if ((*this->annotationResult)[0].getType()[0] == INTERGENIC) {
        return AnnotationString[INTERGENIC];
      }
    }
    
    std::string s;
    if (this->fullTemplate.translate(&s)) {
      fprintf(stderr, "fullTemplate failed translation!\n");
    }
    return s;
  };

  // get per gene annotation
  std::string getGeneAnnotation(const AnnotationResult& res){
    const std::vector<AnnotationType>& type = res.getType();
    const std::map<AnnotationType, std::string>& detail = res.getDetail();

    std::vector<std::string> args;
    std::string s;
    for (size_t i = 0; i < type.size(); i++) {
      s =  AnnotationString[type[i]];
      // for synonymous or non-synonymous, we may output details
      std::map<AnnotationType, std::string>::const_iterator iter;
      iter = detail.find(type[i]);
      if (iter != detail.end()) {
        s += iter->second;
      }
      args.push_back(s);
    }

    this->geneTemplate.add("GENE_NAME", res.getFullName());
    this->geneTemplate.add("GENE_STRAND", res.hasForwardStrand() ? FORWARD_STRAND_STRING : REVERSE_STRAND_STRING);
    this->geneTemplate.add("TYPE", args);

    this->geneTemplate.translate(&s);

    return s;
  };

  int inplace_make_set(std::vector<std::string>* input) {
    std::vector<std::string>& v = * input;
    std::sort(v.begin(), v.end());
    std::vector<std::string>::iterator it;
    it = std::unique(v.begin(), v.end());
    v.resize( it - v.begin());
    return v.size();
  };

 private:
  StringTemplate topPriorityTemplate;
  StringTemplate geneTemplate;
  StringTemplate fullTemplate;
  const AnnotationResultCollection* annotationResult;
  const Priority* priority;
}; // end AnnotationOutput


#endif /* _ANNOTATIONOUTPUT_H_ */
