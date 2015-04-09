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
            : width(width), height(height), h(h), leeches(leeches), method(method),
              xCoordinates(std::round(height / h)),
              yCoordinates(std::round(width / h)),
              temperatures(Matrix(xCoordinates + 1, yCoordinates + 1)) {
    }

    Matrix run() {
        // Dimensiones de la matriz banda de ecuaciones
        int dims = xCoordinates * yCoordinates;
        // La matriz de ecuaciones tiene banda inferior y superior como la cantidad de columnas de
        // la matriz de temperaturas. 
        Matrix system(dims, dims, this->temperatures.columns(), this->temperatures.columns());
        // La solución que buscamos
        BDouble *b = new BDouble[dims];


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

        delete[] b;

        return this->temperatures;
    }
private:
    void band_gaussian_elimination(const Matrix &system, BDouble *b) {
        // Resolver el problema
        std::pair<BDouble *, enum Solutions> solution = gaussian_elimination(system, b);

        // Cargar los datos en la matriz
        for (int i = 0; i < this->temperatures.rows(); ++i) {
            for (int j = 0; j < this->temperatures.columns(); ++j) {
                this->temperatures(i, j) = solution.first[j * this->temperatures.columns() + i];
            }
        }

        // Borrar el espacio extra
        delete[] solution.first;
    }

    void lu_factorization(const Matrix &system, BDouble *b) {
        std::pair<Matrix, Matrix> factors = LU_factorization(system);
        std::pair<BDouble *, enum Solutions> temporal = gaussian_elimination(factors.first, b);
        std::pair<BDouble *, enum Solutions> solution = gaussian_elimination(factors.second, temporal.first);

        // Cargar los datos en la matriz
        for (int i = 0; i <       this->temperatures.rows(); ++i) {
            for (int j = 0; j < this->temperatures.columns(); ++j) {
                this->temperatures(i, j) = solution.first[j * this->temperatures.columns() + i];
            }
        }

        // Liberamos la memoria que usamos.
        delete[] solution.first;
        delete[] temporal.first;
    }
    
    
	
	bool is_border(int x, int y) {
		return x == 0 || y == 0 || x == xCoordinates || y == yCoordinates;
 	}
	
	/**
	 * Cada ecuacion de temperatura corresponde a una temperatura final
	 * con posiciones en xCoordinate * yCoordinate.
	 * 
	 * Devuelve la sanguijuela de mayor temperatura que afecte la posicion x, y
	 * 
	 * */
	Leech* high_temp_leech(int x, int y) {
		 Leech* maxLeech = NULL;
		 //BDouble maxTemp = 0.0;
		 for (std::list<Leech>::iterator leech = leeches.begin(); leech != leeches.end(); ++leech) {
            // Distribuimos las temperaturas de la sanguijuela
            int topX = std::floor(leech->x + leech->radio);
            int bottomX = std::ceil(leech->x - leech->radio);
            int topY = std::floor(leech->y + leech->radio);
            int bottomY = std::ceil(leech->y - leech->radio);

            // Ponemos las coordenadas en rango
            topX = std::min(std::max(topX, 0), xCoordinates);
            bottomX = std::min(std::max(bottomX, 0), xCoordinates);

            topY = std::min(std::max(topY, 0), yCoordinates);
            bottomY = std::min(std::max(bottomY, 0), yCoordinates);

            // Seteamos las temperaturas en la matriz.
            // Cabe destacar, la temperatura de cada sanguijuela es igual para todos los puntos que cubre.
            for (int i = bottomX; i <= topX; ++i) {
                for (int j = bottomY; j <= topY; ++j) {
					if(i == x && y == j) {
						if (maxLeech == NULL) {
							maxLeech = &*leech;
						} else if (maxLeech->temperature < leech->temperature) {
							maxLeech = &*leech;
						}					
					}
                }
            }
        }
        return maxLeech;
	}
	
	/**
     * Construimos:
     * - system, la matriz de ecuaciones que representa la relación de las temperaturas.
     * - b, el vector de resultados que representa las condiciones del sistema.
     * */
    void build_system(Matrix &system,  BDouble *b) {
		// Creamos la matriz inicial del sistema
        // TODO: chequear que esto anda
        int tempsAmount = xCoordinates * yCoordinates;
               
        /**
         * Cada fila del sistema representa la ecuacion que rige la temperatura T(i,j)
         * ijEq = (j * yCoordinates) + i
         * */
        for (int ijEq = 0; ijEq < tempsAmount; ijEq++) {
            system(ijEq,ijEq) = 1.0;            
            
            //TODO: Esto funca?
            int i = ijEq % xCoordinates;
            int j = ijEq / xCoordinates;
            

            if(is_border(i,j)) {
				//Si esta en el borde el valor esta fijo en -100.0 y no hay que usar
				//la ecuacion de laplace
				b[ijEq] = -100.0;
				
			} else {
				Leech* leech = high_temp_leech(i,j);
				if (leech != NULL) {
					//Si la posicion se encuentra en el radio de una sanguijuela
					//la temperatura que afecta la posicion es la de la sanguijuela
					//y no hay que usar la ecuacion de laplace
					b[ijEq] = leech->temperature;
				} else {
					//Finalmente si no es borde ni sanguijuela, hay que usar la
					//ecuacion de laplace.
					//Las posiciones de los bordes se ignoran porque figuran con -100.0
					//y fija el valor.
					b[ijEq] = 0.0;
					if (i > 1) {
						int ijVar = (j * yCoordinates) + (i-1);
						system(ijEq, ijVar) = -1.0 / 4.0;					
					}
					if (i < yCoordinates-1) {
						int ijVar = (j * yCoordinates) + (i+1);
						system(ijEq, ijVar) = -1.0 / 4.0;
					}
					if (j > 1) {
						int ijVar = ( (j-1) * yCoordinates) + i;
						system(ijEq, ijVar) = -1.0 / 4.0;
					}
					if (j < xCoordinates-1) {
						int ijVar = ( (j+1) * yCoordinates) + i;
						system(ijEq, ijVar) = -1.0 / 4.0;
					}
					
				}
			}
            
            
        }
        
	} 
	
		
	/*
	void insert_leech(Matrix &system, BDouble *b, int x, int y, Leech* leech) {
			leeches.push_back(*leech);
			
			// Distribuimos las temperaturas de la sanguijuela
			int topX = std::floor(leech->x + leech->radio);
			int bottomX = std::ceil(leech->x - leech->radio);
			int topY = std::floor(leech->y + leech->radio);
			int bottomY = std::ceil(leech->y - leech->radio);

			// Ponemos las coordenadas en rango
			topX = std::min(std::max(topX, 0), xCoordinates);
			bottomX = std::min(std::max(bottomX, 0), xCoordinates);

			topY = std::min(std::max(topY, 0), yCoordinates);
			bottomY = std::min(std::max(bottomY, 0), yCoordinates);

			// Seteamos las temperaturas en la matriz.
			// Cabe destacar, la temperatura de cada sanguijuela es igual para todos los puntos que cubre.
			for (int i = bottomX; i <= topX; ++i) {
				for (int j = bottomY; j <= topY; ++j) {
					int ijEq = (j * yCoordinates) + i;
					
					if (b[ijEq] < leech->temperature) {
						b[ijEq] = leech->temperature;			
					}
					
					//El tamaño de la banda en la matriz de ecuaciones es la cantidad de columnas
					//del parabrisas discretizado pensado como matriz
					for (int h = 1; h < xCoordinates; h++) {
						//Limpiamos todos los coeficientes de la ecuacion ijE excepto el de la diagonal
						//el valor sera T(i,j) = b[ijEq]
						if (ijEq+h < xCoordinate * yCoordinate ){
							system(ijEq, ijEq+h) = 0.0;
						}
						if (ijEq-h > 0 ){
							system(ijEq, ijEq-h) = 0.0;
						}		
					}	 
			
			
				}
			}
	}
	
	void remove_leech(Matrix &system, BDouble *b, Leech* leech) {
		 bool exists = false;
		 for (std::list<Leech>::iterator l = leeches.begin(); l != leeches.end(); ++l) {
			 if (l == *leech) {
				 exists = true;
			 }
		 }
		 
		 if (exists) {
			 leeches.remove(*leech);
			 build_system(system, b);			 
		 }
		 
	}*/

    // Invariante:
    // al menos 1 sanguijuela
    void simple_algorithm(const Matrix &system, BDouble *b) {
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
    void sherman_morrison(const Matrix &system, BDouble *b) {
		/*
        std::pair<Matrix, Matrix> factors = LU_factorization(system);
        BDouble minTemperature = 0.0;
        std::list<Leech>::iterator minPosition = leeches.begin();

        for (std::list<Leech>::iterator leech = ++(leeches.begin()); leech != leeches.end(); ++leech) {
            // Generamos una copia de la lista de sanguijuelas sin la actual
            std::list<Leech> temporal(leeches);
            temporal.erase(distance(leeches.begin(), leech));

            // Armamos la matriz del sistema

            // Actualizamos el mínimo
            if (minTemperature > ...) {
                minima = leech;
            }
        }

        // Nos fijamos si estamos muertos
        if (minTemperature < 235.0) {
            // No morimos!
        } else {
            // Morimos
        }*/

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
