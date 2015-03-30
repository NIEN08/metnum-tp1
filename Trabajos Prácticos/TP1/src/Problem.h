//
// Created by julian on 3/21/15.
//

#ifndef _TP1_PROBLEM_H_
#define _TP1_PROBLEM_H_ 1

#include <string>
#include <cstdint>
#include "BDouble.h"
#include "Matrix.h"

class Problem {
public:
    enum Method {
        BAND_GAUSSIAN_ELEMINATION,
        LU_FACTORIZATION,
        SIMPLE_ALGORITHM,
        SHERMAN_MORRISON
    };

    Problem(std::string, std::string, enum Method);
    int run();
private:
    std::string input;
    std::string output;
    enum Method method;
    unsigned width;
    unsigned height;
    Matrix temperatures;
    BDouble h;
    unsigned amount;
};


#endif //_TP1_PROBLEM_H_
