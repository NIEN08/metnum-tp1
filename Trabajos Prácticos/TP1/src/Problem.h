//
// Created by julian on 3/21/15.
//

#ifndef _TP1_PROBLEM_H_
#define _TP1_PROBLEM_H_ 1

#include <string>
#include <cstdint>

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
    uint64_t width;
    uint64_t height;
    double h;
    uint64_t amount;
};


#endif //_TP1_PROBLEM_H_
