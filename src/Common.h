#ifndef _COMMON_H_
#define _COMMON_H_

// here we define format or the annotation.
#define FORWARD_STRAND_STRING "+"
#define REVERSE_STRAND_STRING "-"
#define WITHIN_GENE_LEFT_DELIM "("
#define WITHIN_GENE_RIGHT_DELIM ")"
#define VCF_ANNOTATION_START_TAG "ANNO="
#define VCF_ANNOTATION_START_TAG_FULL "ANNOFULL="
#define VCF_INFO_SEPARATOR ";"
#define WITHIN_GENE_SEPARATOR ":"
// #define GENE_SEPARATOR "|"

inline bool hasSuffix(const std::string& s, const std::string& suffix) {
  size_t l = s.size();
  size_t l2 = suffix.size();
  if ( l < l2) return false;
  for (size_t i = 1; i <= l2; ++i) {
    if (suffix[l2 -i ] != s[l-i])
      return false;
  }
  return true;
};




#endif /* _COMMON_H_ */
