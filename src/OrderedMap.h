#ifndef _ORDEREDMAP_H_
#define _ORDEREDMAP_H_

template <class KEY, class TYPE>
class OrderedMap{
  public:
    bool find(const KEY& key) {
        if (this->keyTypeMap.find(key) == this->keyTypeMap.end()){
            return false;
        }
        return true;
    }
    TYPE& operator[] (const KEY& key) {
        if (!this->find(key)){
            this->keyVec.push_back(key);
        }
        return this->keyTypeMap[key];
    }
    void front(KEY* k, TYPE* v) {
        *k = this->keyVec.front();
        *v = this->keyTypeMap[(*k)];
    }
    bool at(unsigned int idx, KEY* k, TYPE* v) {
        if (idx >= this->size()) return false;
        *k = this->keyVec[idx];
        *v = this->keyTypeMap[(*k)];
        return true;
    }
    bool at(unsigned int idx, KEY* k, TYPE* v) const{
        if (idx >= this->size()) return false;
        *k = this->keyVec[idx];

        typename std::map < KEY, TYPE >::const_iterator iter;
        iter = this->keyTypeMap.find(*k);
        if (iter != this->keyTypeMap.end())
          *v = iter->second;
        else
          return false;
        
        return true;
    };
    bool value(const KEY& key, TYPE* val) const{
      typename std::map < KEY, TYPE >::const_iterator iter;
      iter = this->keyTypeMap.find(key);
      if (iter != this->keyTypeMap.end()){
        *val = iter->second;
        return true;
      } else {
        return false;
      }
    }
    size_t size() const { return this->keyVec.size();} ;
    void clear() {
      keyVec.clear();
      keyTypeMap.clear();
    };
    const std::vector<KEY>& getKey() const{
      return this->keyVec;
    };
  private:
    std::vector < KEY > keyVec;
    std::map < KEY, TYPE > keyTypeMap;

};

#endif /* _ORDEREDMAP_H_ */
