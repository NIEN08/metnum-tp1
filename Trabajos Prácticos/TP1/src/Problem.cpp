//
// Created by julian on 3/21/15.
//

#include "Problem.h"
#include <fstream>
#include <iostream>
//using namespace std;

Problem::Problem(enum Method method, std::string input = "entrada.txt", std::string output = "salida.txt") :
        input(input), output(output), method(method), width(0), height(0), h(0.0), amount(0) {
    std::ifstream handle(input, std::ifstream::in);

    handle >> this->width >> this->height >> this->h >> this->amount;

    unsigned xCoordinates = this->width / this->h + 1;
    unsigned yCoordinates = this->height / this->h + 1;

    // Creamos la matriz de temperaturas
    this->temperatures = Matrix(xCoordinates, yCoordinates);

    // Ponemos los bordes en -100.0C
    for (std::size_t i = 0; i < xCoordinates; ++i) {
        if (i == 0 || i == xCoordinates - 1) {
            for (std::size_t j = 0; j < yCoordinates; ++j) {
                this->temperatures(i, j) = -100.0;
            }
        } else {
            this->temperatures(i, 0) = -100.0;
            this->temperatures(i, yCoordinates - 1) = -100.0;
        }
    }

    // Levantamos las posiciones de las sanguijuelas
    for (unsigned i = 0; i < this->amount; ++i) {
        BDouble x, y, r, t;
        handle >> x >> y >> r >> t;

        // Distribuimos las temperaturas de las sanguijuelas
        BDouble ceilY = std::ceil(y), floorY = std::floor(y);
        BDouble ceilX = std::ceil(x), floorX = std::floor(x);

        BDouble ceilYR = std::ceil(y + r), floorY = std::floor(y - r);
        BDouble ceilXR = std::ceil(x + r), floorX = std::floor(x - r);
        std
    }

    switch(method) {
        case BAND_GAUSSIAN_ELEMINATION:
            break;
        case LU_FACTORIZATION:
            break;
        case SIMPLE_ALGORITHM:
            break;
        case SHERMAN_MORRISON:
            break;
    };

    handle.close();
}

int Problem::run() {

}

int main() {
    //Problem p("entrada2.txt", "salida.txt", BAND_GAUSSIAN_ELEMINATION);
}