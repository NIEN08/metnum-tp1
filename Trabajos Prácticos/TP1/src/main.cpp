#include "Problem.h"

int main(int argc, char *argv[]) {
    Problem p(std::string(argv[1]), std::string(argv[2]), argv[3]);
    return p.run();
}