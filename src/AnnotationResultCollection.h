#ifndef _ANNOTATIONRESULTCOLLECTION_H_
#define _ANNOTATIONRESULTCOLLECTION_H_

#include "AnnotationResult.h"

class AnnotationResultCollection{
 public:
  AnnotationResultCollection(): sorted(false) {};
  // STL like accessors
  void push_back(const AnnotationResult& r) {
    if (r.getType().size() > 0) {
      this->data.push_back(r);
      this->sorted = false;
    }
  };
  AnnotationResult& operator[] (int idx) {
    return data[idx];
  }
  const AnnotationResult& operator[] (const int idx) const {
    return data[idx];
  }
  bool empty() const {
    return this->data.empty();
  };
  size_t size() const {
    return this->data.size();
  };
  void clear() {
    this->data.clear();
    this->top.clear();
  };
  // sort all annoations and calculate top annotations
  void sortByPriority(const Priority& p) {
    for (size_t i = 0; i < data.size(); ++i) {
      data[i].sortByPriority(p);
    }
    this->sorted = true;

    top.clear();
    if (this->data.size() == 0) return;

    size_t highestIdx = 0;
    Priority::Level highestPriority = p.getPriority(this->data[highestIdx].getType()[0]); // [0] is the highest priority
    this->top.push_back(this->data[0]);
    for (size_t i = 1; i != this->data.size(); i++) {
      Priority::Level newPriority = p.getPriority(this->data[i].getType()[0]);
      if (newPriority < highestPriority) {
        highestIdx = i;
        highestPriority = newPriority;
        this->top.clear();
        this->top.push_back(this->data[i]);
      } else if (newPriority == highestPriority) {
        this->top.push_back(this->data[i]);
      }
    }
    // sort by exon numbers
    std::sort(this->top.begin(), this->top.end(), AnnotationResultCollection::CompareAnnotationResultByExonNumberGreater);
    // for (size_t i = 0; i < this->top.size(); ++i) {
    //   printf("%zu ", i);
    //   this->top[i].dump();
    // }
    // puts("==========-");
  };
  /// remember to call sort by priority before
  const std::vector<AnnotationResult>& getTopAnnotation() const{
    assert(sorted);
    return this->top;
  };
  const std::vector<AnnotationResult>& getAllAnnotation() const {
    return this->data;
  }
 private:
  static bool CompareAnnotationResultByExonNumberGreater(const AnnotationResult& r1,
                                                         const AnnotationResult& r2) {
    return r1.getExonNumber() > r2.getExonNumber();
  };
  std::vector<AnnotationResult> data;
  std::vector<AnnotationResult> top; // top priority result
  bool sorted;
}; // end class AnnotationResultCollection


#endif /* _ANNOTATIONRESULTCOLLECTION_H_ */
