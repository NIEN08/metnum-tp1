#include "Matrix.h"
#include <iostream>
#include <istream>

int main(int argc, char *argv[]) {
    BandMatrix m(3, 3);
    m(0,0) = 1.0;
    m(1,1) = 1.0;
    m(2,2) = 1.0;

    std::cout << m << std::endl;
    BandMatrix g(3, 3);
    g(0,0) = 1.0;
    g(1,1) = 1.0;
    g(2,2) = 1.0;
    g(1,2) = 6.0;
    m += g * BDouble(123.0) + m;
    std::cout << m + g << std::endl;
    return 0;
}
