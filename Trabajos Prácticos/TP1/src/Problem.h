#ifndef _TP1_PROBLEM_H_
#define _TP1_PROBLEM_H_ 1

#include "Matrix.h"

enum Method {
    BAND_GAUSSIAN_ELIMINATION,
    LU_FACTORIZATION,
    SIMPLE_ALGORITHM,
    SHERMAN_MORRISON
};

class Problem {
public:
    Problem(enum Method method, const Matrix &input)
            : temperatures(input), method(method) { }

    Matrix run() {
        switch (method) {
            case BAND_GAUSSIAN_ELIMINATION:
                break;
            case LU_FACTORIZATION:
                break;
            case SIMPLE_ALGORITHM:
                break;
            case SHERMAN_MORRISON:
                break;
        }
    }
private:
    Matrix temperatures;
    enum Method method;
};


#endif //_TP1_PROBLEM_H_
