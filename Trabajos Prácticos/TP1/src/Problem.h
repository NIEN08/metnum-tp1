#ifndef _TP1_PROBLEM_H_
#define _TP1_PROBLEM_H_ 1

#include <list>
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
                temperatures(i, j) = solution.first[j * temperatures.columns() + i];
            }
        }

        // Borrar el espacio extra
        delete[] solution.first;
    }

    void lu_factorization(Matrix &A, BDouble *b, Matrix &temperatures) {
		std::cout << "executing build_system..." << std::endl;
		build_system(A, b, this->leeches);
		std::cout << "finalizing build_system..." << std::endl;
		// std::cout << A;
		
		// Sea A la matriz del sistema de ecuaciones,
		// factorizamos A = LU con L, U triangulares inferior/superior
		std::cout << "executing LU_factorization..." << std::endl;
        std::pair<Matrix, Matrix> factors = LU_factorization(A);
		std::cout << "finalizing LU_factorization..." << std::endl;
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
                temperatures(i, j) = x[ (i * temperatures.columns()) + j];
            }
        }

        // Liberamos la memoria que usamos.
        delete[] y;
        delete[] x;
    }

    // Invariante:
    // al menos 1 sanguijuela
    void simple_algorithm(const Matrix &system, BDouble *b, Matrix &temperatures) {
        /*Matrix minTemperature(xCoordinates + 1, yCoordinates + 1, 1, 1);
        BDouble minValue = 0.0;

        for (std::list<Leech>::iterator leech = leeches.begin(); leech != leeches.end(); ++leech) {
            // Generamos una copia de la lista de sanguijuelas sin la actual
            std::list<Leech> temporal(leeches);
            temporal.erase(distance(leeches.begin(), leech));

            // Armamos la matriz del sistema
            std::pair<Matrix, BDouble *> tmpSystem = this->generate_system(temporal);

            // Resolvemos
            std::pair<BDouble *, enum Solutions> solution = gaussian_elimination(tmpSystem.first, tmpSystem.second);

            // Liberamos memoria
            delete[] tmpSystem.second;

            // Actualizamos el mínimo
            // TODO: cual es el medio?
            if (solution.first[] > minValue || leech = leeches.begin()) {
                minValue = solution.first[];

                // Cargar los datos en la matriz
                for (int i = 0; i < this->temperatures.rows(); ++i) {
                    for (int j = 0; j < this->temperatures.columns(); ++j) {
                        minTemperature(i, j) = solution.first[j * this->temperatures.columns() + i];
                    }
                }
            }

            delete[] solution.first;
        }

        // Vamos a devolver la matriz de temperatura que tengamos
        this->temperatures = minTemperature;*/
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
		int c = 1;
        for (auto &leech : leeches) {
            std::cout << "leech" << c << "(" << leech.x << " ," << leech.y << " ) : " << leech.temperature << std::endl;
			c++;
			// Distribuimos las temperaturas de la sanguijuela
            int topX = std::floor((leech.x + leech.radio)/h);
            int bottomX = std::ceil((leech.x - leech.radio)/h);
            int topY = std::floor((leech.y + leech.radio)/h);
            int bottomY = std::ceil((leech.y - leech.radio)/h);
			
			
			
            // Ponemos las coordenadas en rango
            topX = std::min(std::max(topX, 0), columns - 1);
            bottomX = std::min(std::max(bottomX, 0), columns - 1);

            topY = std::min(std::max(topY, 0), rows - 1);
            bottomY = std::min(std::max(bottomY, 0), rows - 1);
			
			
			std::cout << "topX " << topX << std::endl;
			std::cout << "bottomX " << bottomX << std::endl;
			std::cout << "topY " << topY << std::endl;
			std::cout << "bottomY " << bottomY << std::endl;
			
            // Seteamos las temperaturas en la matriz.
            // Cabe destacar, la temperatura de cada sanguijuela es igual para todos los puntos que cubre.
            for (int i = bottomX; i <= topX; ++i) {
                for (int j = bottomY; j <= topY; ++j) {
                    //std::cout << i <<" "<<j<< " " << leech.temperature << std::endl;
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

        /**
         * Cada fila del sistema representa la ecuacion que rige la temperatura T(i,j)
         * ijEq = (j * yCoordinates) + i
         * */
        for (int ijEq = 0; ijEq < limit; ijEq++) {
            system(ijEq, ijEq) = 1.0;
            int i = ijEq % columns;
            int j = ijEq / columns;

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
                    if (i > 1) {
                        int ijVar = (j * columns) + i - 1;
                        system(ijEq, ijVar) = -0.25;
                    }

                    if (i < rows - 1) {
                        int ijVar = (j * columns) + i + 1;
                        system(ijEq, ijVar) = -0.25;
                    }

                    if (j > 1) {
                        int ijVar = ((j - 1) * columns) + i;
                        system(ijEq, ijVar) = -0.25;
                    }

                    if (j < columns - 1) {
                        int ijVar = ((j + 1) * columns) + i;
                        system(ijEq, ijVar) = -0.25;
                    }
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