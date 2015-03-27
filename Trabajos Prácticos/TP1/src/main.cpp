#include "Matrix.h"
#include <iostream>
#include <istream>

int main(int argc, char *argv[]) {
    BandMatrix m(3, 3);
    m(0,0) = 2.0;
    m(0,1) = 4.0;
    m(1,0) = 4.0;
    m(1,1) = 2.0;

    BandMatrix n(m);
    std::cout << m + n;
    n(0, 0) = 30.0;

    std::cout << m + n;

    return 0;
}
