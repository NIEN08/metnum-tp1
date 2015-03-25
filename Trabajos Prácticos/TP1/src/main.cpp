#include "Matrix.h"

int main(int argc, char *argv[]) {
    InefficientMatrix<int, 0> m(3, 3);
    m(0,0) = 1;
    m(1,1) = 1;
    m(2,2) = 1;

    std::cout << m << std::endl;
    InefficientMatrix<int, 0> g(3, 3);
    g(0,0) = 1;
    g(1,1) = 1;
    g(2,2) = 1;
    g(1,2) = 6;
    m += g * 123 + m;
    std::cout << m + g << std::endl;
}