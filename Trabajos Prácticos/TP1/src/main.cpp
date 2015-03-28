#include "Matrix.h"

int main(int argc, char *argv[]) {
    Matrix m(3,3, 0, 1);
    m(0,1) = 1.0;
    Matrix g(3,4);
    g(0, 1) = 1.0;
    g(1,1) = 5.6;
    g(2,3) = 7.8;

    std::cout << m;
    std::cout << g;
    std::cout << (m * g);
    return 0;
}
