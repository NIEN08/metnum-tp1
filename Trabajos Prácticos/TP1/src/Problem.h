#ifndef _TP1_PROBLEM_H_
#define _TP1_PROBLEM_H_ 1

#include <list>
#include <algorithm>
#include <map>
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
} Leech;

class Problem {
public:
    // Invariante:
    // h /= 0
    // height /= 0
    // width /= 0
    // h | height
    // h | width
    Problem(enum Method method,
            const BDouble &width,
            const BDouble &height,
            const BDouble &h,
            std::list<Leech> &leeches)
            : width(width), height(height), h(h), leeches(leeches), method(method) {
    }

    Matrix run() {
        int rows = round(height / h) + 1;
        int columns = round(width / h) + 1;
        int dims = rows * columns;

        Matrix system(dims, dims, columns, columns);
        BDouble *b = new BDouble[dims];
        Matrix temperatures(rows, columns);

        switch (method) {
            case BAND_GAUSSIAN_ELIMINATION:
                band_gaussian_elimination(system, b, temperatures);
                break;
            case LU_FACTORIZATION:
                lu_factorization(system, b, temperatures);
                break;
            case SIMPLE_ALGORITHM:
                simple_algorithm(system, b, temperatures);
                break;
            case SHERMAN_MORRISON:
                sherman_morrison(system, b, temperatures);
                break;
        }

        delete[] b;
        return temperatures;
    }
private:
    void band_gaussian_elimination(Matrix &system, BDouble *b, Matrix &temperatures) {
        build_system(system, b, this->leeches);
        std::pair<BDouble *, enum Solutions> solution = gaussian_elimination(system, b);

        // Cargar los datos en la matriz
        for (int i = 0; i < temperatures.rows(); ++i) {
            for (int j = 0; j < temperatures.columns(); ++j) {
                temperatures(i, j) = solution.first[(i * temperatures.columns()) + j];
            }
        }

        // Borrar el espacio extra
        delete[] solution.first;
    }

    void lu_factorization(Matrix &A, BDouble *b, Matrix &temperatures) {
        build_system(A, b, this->leeches);

        // Sea A la matriz del sistema de ecuaciones,
        // factorizamos A = LU con L, U triangulares inferior/superior
        std::pair<Matrix, Matrix> factors = LU_factorization(A);
        Matrix& L = factors.first;
        Matrix& U = factors.second;

        //Resolvemos el sistema Ly = b
        std::pair<BDouble *, enum Solutions> partialSolution = forward_substitution(L, b);
        BDouble* y = partialSolution.first;

        //Resolvemos el sistema Ux = y
        std::pair<BDouble *, enum Solutions> finalSolution = backward_substitution(U, y);
        BDouble* x	= finalSolution.first;

        //Cargamos la solucion en la matriz de temperaturas
        for (int i = 0; i < temperatures.rows(); i++) {
            for (int j = 0; j < temperatures.columns(); j++) {
                temperatures(i, j) = x[(i * temperatures.columns()) + j];
            }
        }

        // Liberamos la memoria que usamos.
        delete[] y;
        delete[] x;
    }

    void simple_algorithm(Matrix &system, BDouble *b, Matrix &temperatures) {

        Matrix minTemperature(system.rows(), system.columns(), 1, 1);
        BDouble minValue = 0.0;

        for (std::list<Leech>::iterator leech = leeches.begin(); leech != leeches.end(); ++leech) {
            // Generamos una copia de la lista de sanguijuelas sin la actual
            std::list<Leech> temporal(leeches);
            temporal.erase(std::distance(leeches.begin(), leech));
            BDouble *b2 = std::copy(b); // TODO:
            Matrix system2(system);

            // Armamos la matriz del sistema
            this->build_system(system, b, temporal);

            // Resolvemos
            std::pair<BDouble *, enum Solutions> solution = gaussian_elimination(system2, b2);

            // Valor del medio
            BDouble n = 0.0;
            BDouble p = 0.0;

            BDouble centerY = this->height/2.0;
            BDouble centerX = this->width/2.0;
            BDouble topJ = (centerX + this->h)/this->h;
            BDouble bottomJ = (centerX - this->h)/this->h;
            BDouble topI = (centerY + this->h)/this->h;
            BDouble bottomI = (centerY - this->h)/this->h;

            for (int i = std::ceil(bottomI); BDouble(double(i)) <= topI; ++i) {
                BDouble iA = BDouble(double(i));

                for (int j = std::ceil(bottomJ); BDouble(double(j)) <= topJ; ++j) {
                    BDouble iJ = BDouble(double(j));
                    BDouble coef = std::pow(iA * this->h - centerY, 2) + std::pow(iJ * this->h - centerX, 2);

                    if (coef <= h*h) {
                        n += 1.0;
                        p += solution.first[i * temperatures.columns() + j];
                    }
                }
            }

            p /= n;

            if (p > minValue || leech == leeches.begin()) {
                minValue = p;

                // Cargar los datos en la matriz
                for (int i = 0; i < temperatures.rows(); ++i) {
                    for (int j = 0; j < temperatures.columns(); ++j) {
                        minTemperature(i, j) = solution.first[i * temperatures.columns() + j];
                    }
                }
            }

            delete[] solution.first;
            delete[] b2;
        }

        // Vamos a devolver la matriz de temperatura que tengamos
        temperatures = minTemperature;
    }

    // Invariante:
    // al menos 1 sanguijuela
    void sherman_morrison(const Matrix &system, BDouble *b, Matrix &temperatures) {
    }

    /**
     * Construimos:
     * - system, la matriz de ecuaciones que representa la relación de las temperaturas.
     * - b, el vector de resultados que representa las condiciones del sistema.
     * */
    void build_system(Matrix &system,  BDouble *b, const std::list<Leech> &leeches) const {
        int columns = system.upper_bandwidth();
        int rows = system.rows() / columns;
        int limit = columns * rows;

        std::map<std::pair<int, int>, BDouble> associations;

        // Cargamos posiciones afectadas por alguna sanguijuela
        for (auto &leech : leeches) {
            // Ponemos el rango que vamos a generar
            BDouble topJ = std::min(leech.x + leech.radio, this->width - this->h)/h;
            BDouble bottomJ = std::max(leech.x - leech.radio, this->h)/h;

            BDouble topI = std::min(leech.y + leech.radio, this->height - this->h)/h;
            BDouble bottomI = std::max(leech.y - leech.radio, this->h)/h;

            // Seteamos las temperaturas en la matriz.
            // Cabe destacar, la temperatura de cada sanguijuela es igual para todos los puntos que cubre.
            for (int i = std::ceil(bottomI); BDouble(double(i)) <= topI; ++i) {
                BDouble iA = BDouble(double(i));

                for (int j = std::ceil(bottomJ); BDouble(double(j)) <= topJ; ++j) {
                    BDouble iJ = BDouble(double(j));
                    BDouble coef = std::pow(iA*this->h - leech.y, 2) + std::pow(iJ*this->h - leech.x, 2);

                    if (coef <= std::pow(leech.radio, 2)) {
                        try {
                            if (associations.at(std::pair<int, int>(i, j)) < leech.temperature) {
                                associations[std::pair<int, int>(i, j)] = leech.temperature;
                            }
                        } catch(...) {
                            associations[std::pair<int, int>(i, j)] = leech.temperature;
                        }
                    }
                }
            }
        }

        for (int ijEq = 0; ijEq < limit; ijEq++) {
            system(ijEq, ijEq) = 1.0;
            int i = ijEq / columns;
            int j = ijEq % columns;

            if (i == 0 || j == 0 || i == rows - 1 || j == columns - 1) {
                //Si esta en el borde el valor esta fijo en -100.0 y no hay que usar
                //la ecuacion de laplace
                b[ijEq] = -100.0;

            } else {
                try {
                    //Si la posicion se encuentra en el radio de una sanguijuela
                    //la temperatura que afecta la posicion es la de la sanguijuela
                    //y no hay que usar la ecuacion de laplace
                    b[ijEq] = associations.at(std::pair<int, int>(i, j));

                } catch(...) {
                    //Finalmente si no es borde ni sanguijuela, hay que usar la
                    //ecuacion de laplace.
                    //Las posiciones de los bordes se ignoran porque figuran con -100.0
                    //y fija el valor.
                    b[ijEq] = 	0.0;

                    // t[i-1][j] + t[i, j-1] - 4*t[i, j] + t[i+1, j] + t[i, j+1] = 0
                    system(ijEq, (i * columns) + j - 1) = -0.25;
                    system(ijEq, (i * columns) + j + 1) = -0.25;
                    system(ijEq, ((i - 1) * columns) + j) = -0.25;
                    system(ijEq, ((i + 1) * columns) + j) = -0.25;
                }
            }
        }
    }

    BDouble width;
    BDouble height;
    BDouble h;
    std::list<Leech> leeches;
    enum Method method;
};


#endif //_TP1_PROBLEM_H_