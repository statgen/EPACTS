#ifndef _STRINGUTIL_H_
#define _STRINGUTIL_H_

#include <algorithm>
/** tokenize the string
    @return number of tokens we obtained
    e.g. For empty input string, we will return 1, and result will have
only 1 element (the empty string)
*/
int stringTokenize(const std::string& str, const std::string& delim, std::vector<std::string>* result){
    assert(result);
    result->clear();
    std::string s;
    unsigned int l = str.size();
    unsigned int i = 0;
    while (i < l) {
        if (delim.find(str[i]) != std::string::npos) { // it's a delimeter
            result->push_back(s);
            s.clear();
        } else {
            s.push_back(str[i]);
        }
        ++i;
    };
    result->push_back(s);    
    return result->size();
};
int stringTokenize(const std::string& str, const char delim, std::vector<std::string>* result){
    std::string d;
    d.push_back(delim);
    return (stringTokenize(str, d, result));
};

bool isEmptyString(const std::string& s){
    return s.size() == 0;
};
int stringNaturalTokenize(const std::string& str, const std::string& delim, std::vector<std::string>* result){
    int ret = stringTokenize(str, delim, result);
    ret = remove_if(result->begin(), result->end(), isEmptyString) - result->begin();
    result->resize(ret);
    return ret;
};
int stringNaturalTokenize(const std::string& str, const char delim, std::vector<std::string>* result){
    std::string d;
    return (stringNaturalTokenize(str, d, result));
};

//remove leading and trailing characters
void stringStrip(std::string* input, const char* characters = " ") {
  if (!input || !input->size()) return;
    size_t beg = input->find_first_not_of(characters);
    size_t end = input->find_last_not_of(characters);
    input->assign( input->substr(beg, end - beg + 1) );
};

//extract piece of string from @param beg(inclusive) to @param end(exclusive)
//NOTE: beg and end can be negative meaning count from right hand side. 
void stringSlice(std::string* input, int beg, int end) {
    assert(input);
    unsigned int len = input->size();
    if (beg < 0) beg += len;
    if (end < 0) end += len;
    assert (beg >= 0 && end >= 0);
    input -> assign ( input->substr(beg, end- beg)) ;
};

template <class T>
std::string stringJoin(const std::vector<std::string>& input, const T delim) {
    std::string s;
    if (input.size() == 0) {
        return s;
    }
    s = input[0];
    for (unsigned int i = 1; i < input.size(); i++) {
        s+= delim;
        s+= input[i];
    }
    return s;
};
/**
 * for std::string type, we use reference to save memory.
 */
template <>
std::string stringJoin<const std::string&>(const std::vector<std::string>& input, const std::string& delim) {
    std::string s;
    if (input.size() == 0) {
        return s;
    }
    s = input[0];
    for (unsigned int i = 1; i < input.size(); i++) {
        s+= delim;
        s+= input[i];
    }
    return s;
};

std::string toUpper(const std::string& s) {
    std::string r;
    for (unsigned int i = 0; i < s.size(); i++) {
        r.push_back(toupper(s[i]));
    }
    return r;
};

std::string toLower(const std::string& s) {
    std::string r;
    for (unsigned int i = 0; i < s.size(); i++) {
        r.push_back(tolower(s[i]));
    }
    return r;
};

// remove the leading 'chr' if any
std::string chopChr(const std::string& s) {
    if (s.size() > 3 && 
        (s[0] == 'c' || s[0] == 'C') &&
        (s[1] == 'h' || s[1] == 'H') &&
        (s[2] == 'r' || s[2] == 'R')){
        return s.substr(3);
    }
    return s;
};

char _bufferStr[128];
// convert number to char*
const char* toStr(const int i) {
    sprintf(_bufferStr, "%d", i);
    return _bufferStr;
};
const char* toStr(const double d) {
    sprintf(_bufferStr, "%lf", d);
    return _bufferStr;
};

void tolower(std::string* s) {
  for (std::string::iterator i = s->begin();
       i != s->end();
       ++i)
    (*i) = tolower(*i);
};

std::string tolower(const std::string& s) {
  std::string ret(s);
  tolower(&ret);
  return ret;
};
#endif /* _STRINGUTIL_H_ */
