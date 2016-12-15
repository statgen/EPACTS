#ifndef _TYPECONVERSION_H_
#define _TYPECONVERSION_H_

#include <errno.h>
#include <limits.h>
#include <math.h> // for HUGE_VALH, HUGE_VALL
#include <sstream>

// convert double/int/byte to string type
template<class T>
std::string toString(T i){
    std::stringstream ss;
    ss << i;
    return ss.str();
}

// convert std::string to integer
// @return true if conversion succeed
bool str2int(const char* input, int* output) {
    char* endptr;
    long val;
    errno = 0;
    val = strtol(input, &endptr, 10);
    
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
        || (errno != 0 && val == 0)) {
        perror("strtol");
        return false;
    }

    if (endptr == input) {
        // no digits found
        return false;
    }
    *output = val;
    return true;
}

bool str2int(const std::string& input, int* output) {
    return str2int(input.c_str(), output);
}

int toInt(const char* input) {
    return (atoi(input));
};
int toInt(const std::string& input) {
    return (atoi(input.c_str()));
};

float toFloat(const char* input) {
    return (atof(input));
};
float toFloat(const std::string& input) {
    return (atof(input.c_str()));
};

double toDouble(const char* input) {
    return (atof(input));
};
double toDouble(const std::string& input) {
    return (atof(input.c_str()));
};

// convert std::string to double
// @return true if conversion succeed
bool str2double(const char* input, double* output) {
    char* endptr;
    double val;

    errno = 0;
    val = strtod(input, &endptr);
    
    if ((errno == ERANGE && (val == HUGE_VALF || val == HUGE_VALL))
        || (errno != 0 && val == 0.)) {
        perror("strtod");
        return false;
    }

    if (endptr == input) {
        // no digits found
        return false;
    }
    *output = val;
    return true;
}

#endif /* _TYPECONVERSION_H_ */
