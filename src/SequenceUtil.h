#ifndef _SEQUENCEUTIL_H_
#define _SEQUENCEUTIL_H_

/**
 * @return the complement of @param c in capital letter
 */
char complementBase(const char c){
    switch (c){
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    case 'a':
        return 'T';
    case 't':
        return 'A';
    case 'g':
        return 'C';
    case 'c':
        return 'G';
    default:
        return 'N';
    }
};

void complementTriplet(char s[3]){
    s[0] = complementBase(s[0]);
    s[1] = complementBase(s[1]);
    s[2] = complementBase(s[2]);
};

void reverseComplementTriplet(char s[3]){
    char t = s[0];
    s[0] = complementBase(s[2]);
    s[2] = complementBase(t);
    s[1] = complementBase(s[1]);
};

#endif /* _SEQUENCEUTIL_H_ */
