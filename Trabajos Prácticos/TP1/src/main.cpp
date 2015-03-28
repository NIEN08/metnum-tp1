#include "Matrix.h"

int main(int argc, char *argv[]) {
    Matrix m(3, 3);
    m(0,0) = 2.0;
    m(0,1) = 4.0;
    m(1,0) = 4.0;
    m(1,1) = 2.0;

    BDouble b[] = {2.0, 5.0, 0.0};
    std::pair<BDouble *, enum Solutions> sol = m.gaussian_elimination(b);

    for (std::size_t k = 0; k < 3; ++k) {
        std::cout << sol.first[k] << " ";
    }
    std::cout << std::endl;

    return 0;
}
