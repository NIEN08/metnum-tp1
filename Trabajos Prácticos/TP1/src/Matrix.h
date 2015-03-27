//
// Created by julian on 3/21/15.
//

#ifndef _TP1_MATRIX_H_
#define _TP1_MATRIX_H_ 1

#include <cstdint>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include "BDouble.h"
#include <limits>

/*
* Matriz Banda.
*/
class BandMatrix {
public:
    BandMatrix(std::size_t, std::size_t, std::size_t lband = SIZE_MAX, std::size_t uband = SIZE_MAX);
    BandMatrix(const BandMatrix &);
    std::size_t rows() const;;
    std::size_t columns() const;
    std::size_t upper_bandwidth() const;
    std::size_t lower_bandwidth() const;
    BandMatrix & operator=(const BandMatrix &);
    bool operator==(const BandMatrix &) const;
    bool operator!=(const BandMatrix &) const;
    BDouble & operator()(std::size_t, std::size_t);
    const BDouble & operator()(std::size_t, std::size_t) const;
    BandMatrix & operator+=(const BandMatrix &);
    BandMatrix & operator*=(const BDouble &);
    BDouble * triangulate(BDouble *);

    virtual ~BandMatrix();
private:
    // Matrix
    std::size_t N;
    std::size_t M;
    std::size_t uband;
    std::size_t lband;
    BDouble **matrix;
};

extern std::ostream &operator<<(std::ostream &, const BandMatrix &);
extern BandMatrix operator+(const BandMatrix &, const BandMatrix &);
extern BandMatrix operator*(const BandMatrix &, const BDouble &);
extern BandMatrix operator*(const BandMatrix &, const BandMatrix &);

#endif //_TP1_MATRIX_H_
