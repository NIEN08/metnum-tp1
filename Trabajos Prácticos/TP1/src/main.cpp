#include "Matrix.h"
#include <iostream>
#include <istream>

int main(int argc, char *argv[]) {
    BandMatrix m(2, 2);
    m(0,0) = 2.0;
    m(0,1) = 4.0;
    m(1,0) = 4.0;
    m(1,1) = 2.0;

    BandMatrix a(2, 2);
    a(0,0) = 3.0;
    a(0,1) = 2.0;
    a(1,0) = 2.0;
    a(1,1) = 3.0;
    std::cout << m*a << std::endl;

    return 0;
}
