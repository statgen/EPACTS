#ifndef _PRIORITY_H_
#define _PRIORITY_H_

#include "AnnotationString.h"

extern OutputAnnotationString AnnotationString; // global variable

/**
 * priority are given by lineNo
 * small number -> high priority (more important)
 * contain a relationship between AnnotationType and int(priority)
 */
struct Priority{
 public:
  // Level is just a wrapper for int, we use it only for type checking.
  struct Level{
    Level(int i):level(i){};
    bool operator<(const Priority::Level& l) const{
      // fprintf(stderr, "called <1: %d, %d \n", this->level, l.level);
      return this->level < l.level;
    };
    bool operator==(const Priority::Level& l) const{
      return this->level == l.level;
    };
    void dump() {
      fprintf(stderr, "level = %d\n", this->level);
    };
    int level;
  };     // end struct Level
  int open(const char* fileName) {
    // Load priority.txt
    this->priorityIdx = 0;
    this->priorityInt2Str.clear();
    this->priorityStr2Int.clear();
    LineReader lr(fileName);
    std::vector<std::string> fd;

    while (lr.readLineBySep(&fd, " \t")){
      if (fd.size() == 0) continue;
      if (fd[0][0] == '#') continue;
      if (fd[0].size() == 0) continue;
      priorityIdx ++;
      // fprintf(stderr, "add priority [%s]\n", fd[0].c_str());
      priorityInt2Str[priorityIdx] = fd[0];
      priorityStr2Int[fd[0]] = priorityIdx;
    }
    return priorityIdx;
  };
  Level getPriority(const AnnotationType& t) const{
    std::map<std::string, int>::const_iterator it;
    it = this->priorityStr2Int.find( AnnotationString[t] );
    if (it == this->priorityStr2Int.end()) {
      fprintf(stderr, "Cannot find annotation type [ %s ] from priority files!\n", AnnotationString[t]);
      Level l(-1);
      return l;
    } else {
      Level l(it->second);
      return l;
    }
  };

  std::string getAnnotationString (const int& i) const{
    std::map<int, std::string>::const_iterator it;
    it = this->priorityInt2Str.find( i );
    if (it == this->priorityInt2Str.end()) {
      fprintf(stderr, "Cannot find priority [ %d ] from priority files!\n", i);
      return "";
    } else {
      return it->second;
    }
  };
  /**
   * @return the string representation of @param priority
   */
  std::string toString(const int priority) const {
    std::string s;
    std::map<int, std::string>::const_iterator it;
    it = this->priorityInt2Str.find(priority);
    if (it == this->priorityInt2Str.end()) {
      return s;
    } else{
      return it->second;
    }
  };

  static int getLeastPriority(){
    return 9999;
  }
 private:
  int priorityIdx;
  std::map<int, std::string> priorityInt2Str;
  std::map<std::string, int> priorityStr2Int;
}; // end Priority



#endif /* _PRIORITY_H_ */
