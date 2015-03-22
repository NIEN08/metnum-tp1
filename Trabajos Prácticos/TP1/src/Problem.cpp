//
// Created by julian on 3/21/15.
//

#include "Problem.h"
#include <fstream>

Problem::Problem(std::string input, std::string output, enum Method method) :
        input(input), output(output), method(method), width(0), height(0), h(0.0), amount(0) {
    std::ifstream handle(input, std::ifstream::in);
    handle >> width >> height >> h >> amount;

    auto i = 0;

    // Espero que esto tire exceptions...
    while (i < amount) {
        double x = 0, y = 0, r = 0, t = 0;
        handle >> x >> y >> r >> t;
        // Acá hay que meter los puntos en algún lado que sea útil
        ++i;
    }
}

int Problem::run() {

}