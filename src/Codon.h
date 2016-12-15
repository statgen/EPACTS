#ifndef _CODON_H_
#define _CODON_H_

class Codon{
public:
  /// @return 0 if success
  bool open(const char* codonFile) {
    LineReader lr(codonFile);
    std::string line;
    std::vector<std::string> f;
    while(lr.readLine(&line)>0) {
      if (line.size() > 0 && line[0] != '#') {
        stringTokenize(line, '\t', &f);
        this->codon2aa[f[0]] = f[1];
        this->codon2letter[f[0]] = f[2];
        this->codon2fullName[f[0]] = f[3];
      };
    }
    return true;
  };
  const std::string& toAA(const char s[3]) {
    std::string key;
    key.push_back(s[0]);
    key.push_back(s[1]);
    key.push_back(s[2]);
    return safeAccess(this->codon2aa, key, Codon::unknownAA);
  };
  const std::string& toLetter(const char s[3]) {
    std::string key;
    key.push_back(s[0]);
    key.push_back(s[1]);
    key.push_back(s[2]);
    return safeAccess(this->codon2letter, key, Codon::unknownLetter);
  };
public:
  static bool isStopCodon(const std::string& a) {
    if (a == "Stp") return true;
    return false;
  };
public:
  static std::string unknownAA;
  static std::string unknownLetter;
  static std::string unknownFullName;
private:
  /**
   * @return if the key exists, return data[key]
   * or @return @param defaultValue
   */
  const std::string& safeAccess( const std::map<std::string, std::string>& data,
                                 const std::string& key,
                                 const std::string& defaultValue) const {
    std::map<std::string, std::string>::const_iterator it = data.find(key);
    if (it == data.end())
      return defaultValue;
    return it->second;
  }

private:
  std::map<std::string, std::string> codon2aa;        // three letter amino acid
  std::map<std::string, std::string> codon2letter;    // amino acio letter
  std::map<std::string, std::string> codon2fullName;  // full amino acid name

};
std::string Codon::unknownAA = "N/A";
std::string Codon::unknownLetter = "*";
std::string Codon::unknownFullName ="UnknownAminoAcid";


#endif /* _CODON_H_ */
