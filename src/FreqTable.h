#ifndef _FREQTABLE_H_
#define _FREQTABLE_H_


template <class T>
class FreqTable{
public:
    void add(const T& t) {
        if (this->data.find(t) == this->data.end()) {
            this->data[t] = 1;
        } else {
            this->data[t] ++; 
        }
        this->isSorted = false;
    };
    void remove(const T& t) {
        if (this->data.find(t) == this->data.end()) {
            return;
        }
        this->data[t] -- ;
        this->isSorted = false;
    };
    size_t size() const{ return this->data.size();}; 
    // return the frequency in ascending order
    void at(const unsigned int idx, T* t, int* v) {
        if (!this->isSorted) 
            this->sortByFrequency();
        *v = this->orderedData[idx].first;
        *t = *(this->orderedData[idx].second);
    };
    void clear() {
        this->data.clear();
        this->orderedData.clear();
        this->isSorted = false;
    };
    void dump() {
      for (size_t i = 0; i < this->orderedData.size(); ++i) {
        std::cout << i << "\t"
                  << orderedData[i].first << "\t"
                  << *(orderedData[i].second) << "\n";
      }
    };
      
private:
    void sortByFrequency() {
        this->sortByKey();
        std::stable_sort(this->orderedData.begin(), this->orderedData.end(), FreqTable::sortFirstInPair);
        this->isSorted = true;
        /* dump();  */
    };
    void sortByKey() {
        this->orderedData.clear();
        typename std::map<T, int >::iterator it;
        for (it = this->data.begin(); 
             it != this->data.end() ; it++) {
            this->orderedData.push_back(std::make_pair( (*it).second, &((*it).first)) );
        }
        /* dump(); */
        /* std::stable_sort(this->orderedData.begin(), this->orderedData.end(), FreqTable::sortSecondInPair); */
        /* dump(); */
        this->isSorted = true;
    };
    static bool sortFirstInPair( const std::pair<int, const T*>& t1,
                                 const std::pair<int, const T*>& t2) {
      return t1.first < t2.first;
    };
    static bool sortSecondInPair( const std::pair<int, const T*>& t1,
                                 const std::pair<int, const T*>& t2) {
      return *(t1.second) < *(t2.second);
    };
    std::map< T, int > data;
    std::vector< std::pair<int, const T*> > orderedData;
    bool isSorted;
 };

//////////////////////////////////////////////////////////////////////
// Test case 1 //
/*
    FreqTable<std::string> codonFreq;
    codonFreq.add("A");
    codonFreq.add("B");
    codonFreq.add("C");
    codonFreq.add("D");
    codonFreq.add("A");
    codonFreq.add("B");
    codonFreq.add("C");
    codonFreq.add("A");
    codonFreq.add("B");
    codonFreq.add("A");
    
    std::string s;
    int f;
    for (unsigned int i = 0 ; i < codonFreq.size(); i++) {
        codonFreq.at(i, &s, &f);
        printf("freq of %s is %d\n", s.c_str(), f);
    };
    return 0;
*/

// Test case 2 //
/*
    FreqTable<char> codonFreq;
    codonFreq.add('A');
    codonFreq.add('B');
    codonFreq.add('C');
    codonFreq.add('D');
    codonFreq.add('A');
    codonFreq.add('B');
    codonFreq.add('C');
    codonFreq.add('A');
    codonFreq.add('B');
    codonFreq.add('A');
    
    char s;
    int f;
    for (unsigned int i = 0 ; i < codonFreq.size(); i++) {
        codonFreq.at(i, &s, &f);
        printf("freq of %c is %d\n", s, f);
    };
    return 0;

*/

#endif /* _FREQTABLE_H_ */
