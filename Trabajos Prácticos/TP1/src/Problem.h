#ifndef _TP1_PROBLEM_H_
#define _TP1_PROBLEM_H_ 1

#include <list>
#include <cmath>
#include "Matrix.h"

enum Method {
    BAND_GAUSSIAN_ELIMINATION,
    LU_FACTORIZATION,
    SIMPLE_ALGORITHM,
    SHERMAN_MORRISON
};

typedef struct _Leech {
public:
    BDouble x;
    BDouble y;
    BDouble radio;
    BDouble temperature;
    BDouble tempPC;
} Leech;

class Problem {
public:
    Problem(enum Method method,
            const BDouble &width,
            const BDouble &height,
            const BDouble &h,
            std::list<Leech> &leeches)
            : width(width), height(height), h(h), leeches(leeches), method(method),
              xCoordinates(static_cast<int>(std::round(height.get() / h.get()))),
              yCoordinates(static_cast<int>(std::round(width.get() / h.get()))),
              temperatures(Matrix(xCoordinates + 1, yCoordinates + 1)) {

        for (std::list<Leech>::iterator b = leeches.begin(); b != leeches.end(); ++b) {
            // Distribuimos las temperaturas de la sanguijuela
            int topX = std::floor(b->x + b->radio);
            int bottomX = std::ceil(b->x - b->radio);
            int topY = std::floor(b->y + b->radio);
            int bottomY = std::ceil(b->y - b->radio);

            // Ponemos las coordenadas en rango
            topX = std::min(std::max(topX, 0), xCoordinates);
            bottomX = std::min(std::max(bottomX, 0), xCoordinates);

            topY = std::min(std::max(topY, 0), yCoordinates);
            bottomY = std::min(std::max(bottomY, 0), yCoordinates);

            // Seteamos las temperaturas en la matriz.
            // Cabe destacar, la temperatura de cada sanguijuela es igual para todos los puntos que cubre.
            for (int i = bottomX; i <= topX; ++i) {
                for (int j = bottomY; j <= topY; ++j) {
                    // Sólo queda la temperatura más alta
                    if (this->temperatures(i, j) < b->temperature) {
                        this->temperatures(i, j) = b->temperature;
                    }
                }
            }
        }

        // Ponemos los bordes en -100.0C
        for(int i = 0; i <= xCoordinates; ++i){
            if (i == 0 || i == xCoordinates) {
                for (int j = 0; j <= yCoordinates; ++j) {
                    temperatures(i, j) = -100.0;
                }
            } else {
                temperatures(i, 0) = -100.0;
                temperatures(i, yCoordinates) = -100.0;
            }
        }
    }

    Matrix run() {
        // Dimensiones de la matriz banda de ecuaciones
        int dims = this->temperatures.rows() * this->temperatures.columns();
        // La matriz en sí
        Matrix system(dims, dims, 2, 2);
        // La solución que buscamos
        BDouble *b = new BDouble[dims];

        // Creamos la matriz inicial del sistema
        // TODO: chequear que esto anda
        for (int d = 0; d < dims; ++d) {
            for (int j = std::min(0, d - 4); d < std::min(d + 4, this->temperatures.columns()); ++j) {
                system(d, j) = 1.0;
            }

            system(d, d) = 0.0;

            int j = d % this->temperatures.columns();
            int k = j % this->temperatures.rows();

            b[d] = this->temperatures(k, j) * 4.0;
        }

        // Hacemos el "pase de variables al otro lado" para conseguir el b que buscamos.

        // Variamos nuestra solución según el método
        switch (method) {
            case BAND_GAUSSIAN_ELIMINATION:
                this->band_gaussian_elimination(system, b);
                break;
            case LU_FACTORIZATION:
                this->lu_factorization(system, b);
                break;
            case SIMPLE_ALGORITHM:
                this->simple_algorithm(system, b);
                break;
            case SHERMAN_MORRISON:
                this->sherman_morrison(system, b);
                break;
        }

        return this->temperatures;
    }
private:
    void band_gaussian_elimination(const Matrix &system, BDouble *b) {
        // Resolver el problema
        std::pair<BDouble *, enum Solutions> solution = gaussian_elimination(system, b);

        // Cargar los datos en la matriz

        // Borrar el espacio extra
        delete[] solution.first;
    }

    void lu_factorization(const Matrix &system, BDouble *b) {
        std::pair<Matrix, Matrix> factors = LU_factorization(system);
        std::pair<BDouble *, enum Solutions> j = gaussian_elimination(factors.first, b);
        std::pair<BDouble *, enum Solutions> x = gaussian_elimination(factors.second, j.first);
        // x tiene la solución definitiva al sistema.

        // Liberamos la memoria que usamos.
        delete[] j.first;
        delete[] x.first;
    }

    void simple_algorithm(const Matrix &system, BDouble *b) {

    }

    void sherman_morrison(const Matrix &system, BDouble *b) {

    }

    BDouble width;
    BDouble height;
    BDouble h;
    std::list<Leech> leeches;
    enum Method method;
    int xCoordinates;
    int yCoordinates;
    Matrix temperatures;
};


#endif //_TP1_PROBLEM_H_
