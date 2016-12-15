#ifndef _ANNOTATIONRESULT_H_
#define _ANNOTATIONRESULT_H_

#include "Priority.h"

/**
 * For each gene, we use AnnotationResult to store all annotation results.
 */
class AnnotationResult{
 public:
  AnnotationResult() {
    this->clear();
  }
  void clear() {
    this->gene = NULL;
    this->type.clear();
    this->detail.clear();
  };

  void add(const Gene& g){
    if (this->type.size() != 0) {
      // we usually record gene name and its strand first
      fprintf(stderr, "Something weired happen\n");
    }
    this->gene = &g;
  };
  void add(const AnnotationType& t) {
    this->type.push_back(t);
  };
  // add extra details such as "(CCT/Pro->CAT/His)" to the last element
  void addDetail(const AnnotationType& t, const std::string& s) {
    this->detail[t] = s;
  };
  void sortByPriority(const Priority& p) {
    Comparator compareFunction(p);
    // for (size_t i = 0; i < this->type.size(); ++i) {
    //   fprintf(stderr, "%zu -> %s ", i, AnnotationString[type[i]]);
    //   p.getPriority(type[i]).dump();
    // };
    // puts("-----");
    // fprintf(stderr, " sort %d elements: \n", this->type.end() - this->type.begin());
    std::sort(this->type.begin(), this->type.end(), compareFunction);
    // for (size_t i = 0; i < this->type.size(); ++i) {
    //   fprintf(stderr, "%zu -> %s ", i, AnnotationString[type[i]]);
    //   p.getPriority(type[i]).dump();
    // };
    // exit(1);
  };

  //////////////////////////////////////////////////////////////////////
  // getters
  const Gene& getGene() const {
    return (*this->gene);
  };
  const std::string& getGeneName() const{
    return this->gene->geneName;
  };
  const std::string& getTranscriptName() const {
    return this->gene->transcriptName;
  };
  const std::string getFullName() const {
    return this->gene->geneName + "/" + this->gene->transcriptName;
  };
  const size_t getExonNumber() const {
    return this->gene->getExonNumber();
  };
  const std::vector<AnnotationType>& getType() const{
    return this->type;
  };
  const std::map<AnnotationType, std::string>& getDetail() const {
    return this->detail;
  };
  bool hasForwardStrand() const {
    return this->gene->forwardStrand;
  };
  void dump() const { //debugging code
    printf("[ %s ]", gene->geneName.c_str());
    for (size_t i = 0; i < type.size(); ++i) {
      printf(" %s ", AnnotationString[type[i]]);
    };
    puts("");
  };
 private:
  const Gene* gene;
  std::vector<AnnotationType> type;
  std::map<AnnotationType, std::string> detail;

  struct Comparator{
    Comparator(const Priority& p): priority(p) {
      //fprintf(stderr, "create priority comparator\n");
    };
    bool operator() (const AnnotationType& t1,
                     const AnnotationType& t2) const {
      // this->priority.getPriority(t1) .dump();
      // this->priority.getPriority(t2) .dump();
      // fprintf(stderr, "do comparing\n");

      // const Priority::Level& l1 = this->priority.getPriority(t1);
      // const Priority::Level& l2 = this->priority.getPriority(t2);
      // return (l1 < l2 );

      return (this->priority.getPriority(t1)) < (this->priority.getPriority(t2));
    };
   private:
    const Priority& priority;
  };

  // int topPriorityIndex;  // this->type[this->topPriorityIndex] has top priority; <0, means unknown, will remove soon
}; // end AnnotationResult



#endif /* _ANNOTATIONRESULT_H_ */
