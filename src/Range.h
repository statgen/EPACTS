#ifndef _RANGE_H_
#define _RANGE_H_

// all are 1-based index, inclusive on boundaries.
// NOTE: UCSC use 0-based start and 1-based end, this is convenient to calculate length
//       but for complex case, this causes confusion.
/**
 * Example, 1-index range 3-4, inclusive on the boundary
 * in UCSC setting, this range is coded as 2-4 (see below)
 *             @
 * 0-base: 0 1 2 3 4 5 
 * 1-base: 1 2 3 4 5 6
 *               ^
 * the tricky thing is 0-length range, in such case,
 * the UCSC coded the start and end as the same value
 */
struct Range{
    int start;
    int end;
    Range(int s, int e): start(s), end(e) {};
    Range(): start(-1), end(-1) {};
    int length() const { 
        int l = end - start + 1;
        if (l < 0) {
            printf("getLength() < 0 for start(%d) and end(%d)\n", start, end);
        }
        return l;
    };
    bool isInRange(int pos) {
        if (this->start <= pos && pos <= this->end) 
            return true;
        return false;
    };
};

#endif /* _RANGE_H_ */
