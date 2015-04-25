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

		std::cout << method << std::endl;
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
                sherman_morrison_solution(system, b, temperatures);
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
/*
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
        temperatures = minTemperature;*/
    }

    // Invariante:
    // al menos 1 sanguijuela
    void sherman_morrison_solution(Matrix &system, BDouble *b, Matrix &temperatures) {
		std::list<Leech> singularLeeches;
		build_system(system, b, this->leeches);
		
		//Resolvemos el sistema sin resolver ninguna sanguijuela
		//Inicializamos el sistema	
		BDouble* b2 = new BDouble[system.columns()];
		for (int i = 0; i < system.columns(); i++) {
			b2[i] = b[i]; 
		}
		//Resolvemos el sistema
		std::pair<BDouble *, enum Solutions> solution = gaussian_elimination(system, b2);
		delete[] b2;
		
		//Solucion sacando una sanguijuela que minimiza la temperatura del punto critico
		BDouble* minX = solution.first;
		//for (int i = 0; i < system.lower_bandwidth() * system.upper_bandwidth(); i++) {
		//	std::cout << " " << minX[i];
		//} 
		//Valor del punto critico
		BDouble minCp = critic_point_temperature(system, minX);
		std::cout << std::endl << "Solucion original " << std::endl;
		std::cout << std::endl << "CPT: " << minCp << std::endl;
			
		//Resolvemos los sistemas dejando aparte las sanguijuelas de radio 1
		//Se resuelve el sistema cada vez, quitando una sanguijuela.
		int leech_index = 0 ;
		for (std::list<Leech>::iterator itLeech = leeches.begin(); itLeech != leeches.end(); ++itLeech) {
			Leech leech = *itLeech;
			if(	is_singular_leech(leech)){
				std::cout << "is singular!!" << std::endl;
				//La sanguijuela tiene radio 1
				//El sistema removiendo esta sanguijuela se resuelve aparte.
				singularLeeches.push_back(leech);
			} else {
				std::cout << "not singular!!" << std::endl;
				//Armamos una lista sin la sanguijuela
				std::list<Leech> leechesCopy = std::list<Leech>(leeches);
				std::list<Leech>::iterator itCopy = leechesCopy.begin();
				advance(itCopy, leech_index); 
				
				std::cout << "Solucion removiendo sanguijuela (no singular) " << std::endl;
				std::cout << itCopy->x << " " << itCopy->y << " " << itCopy->radio << " " <<  itCopy->temperature << std::endl;
				leechesCopy.erase(itCopy);
				
				// Inicializamos el sistema	
				BDouble* b2 = new BDouble[system.columns()];
				for (int i = 0; i < system.columns(); i++) {
					b2[i] = b[i]; 
				}
				clean_system(system);
				build_system(system, b2, leechesCopy);
			
				//Resolvemos el sistema
				std::pair<BDouble *, enum Solutions> solution = gaussian_elimination(system, b2);
				delete[] b2;
				BDouble* x2 = solution.first;
				BDouble cp = critic_point_temperature(system, x2);
				std::cout << "CPT: " << cp << std::endl;
				// Nos quedamos con la solucion si es mejor que la anterior 
				if(cp < minCp) {
					minX = x2;
					minCp = cp;
				} else {
					delete[] x2;
				}
			}
			leech_index++;
		}
		
		//Aprovechamos la identidad de Sherman-Morrison para resolver los sistemas
		//con sanguijuelas de radio 1
		std::pair<Matrix, Matrix> factors = LU_factorization(system);
        Matrix& L = factors.first;
        Matrix& U = factors.second;
		for (std::list<Leech>::iterator itLeech = leeches.begin(); itLeech != leeches.end(); ++itLeech) {
			Leech leech = *itLeech;
			
			std::cout << "Solucion SM (sanguijuela singular) " << std::endl;
			std::cout << leech.x << " " << leech.y << " " << leech.radio << " " <<  leech.temperature << std::endl;
			
			//TODO: Como calcular las coordenadas?
            int i = std::floor((leech.y + leech.radio)/h);
			int j = std::floor((leech.x + leech.radio)/h);
			
			std::pair<BDouble *, BDouble *> uv = generate_sherman_morrison_uv(system, i, j);
			BDouble* u = uv.first;
			BDouble* v = uv.second;
			
			std::pair<BDouble *, enum Solutions> solution = sherman_morrison(L, U, u, v, b);
			BDouble* x2 = solution.first;
			BDouble cp = critic_point_temperature(system, x2);
			std::cout << "CPT: " << cp << std::endl;
			delete[] u;
			delete[] v;

			// Nos quedamos con la solucion si es mejor que la anterior 
			if(cp < minCp) {
				delete[] minX;
				minX = x2;
				minCp = cp;
			} else {
				delete[] x2;
			}
			
			
		}
	}
	
	/**
	* Devuelve true si la sanguijuela solo afecta una ecuacion del sistema.
	*/
	bool is_singular_leech (Leech leech) {
		// Ponemos el rango que vamos a chequear
		BDouble topJ = std::min(leech.x + leech.radio, this->width - this->h)/h;
		BDouble bottomJ = std::max(leech.x - leech.radio, this->h)/h;
		BDouble topI = std::min(leech.y + leech.radio, this->height - this->h)/h;
		BDouble bottomI = std::max(leech.y - leech.radio, this->h)/h;

		int	coordinates_count = 0;
		for (int i = std::ceil(bottomI); BDouble(double(i)) <= topI; ++i) {
			BDouble iA = BDouble(double(i));

			for (int j = std::ceil(bottomJ); BDouble(double(j)) <= topJ; ++j) {
				BDouble iJ = BDouble(double(j));
				BDouble coef = std::pow(iA*this->h - leech.y, 2) + std::pow(iJ*this->h - leech.x, 2);

				if (coef <= std::pow(leech.radio, 2)) {
					coordinates_count ++;
				}
			}
		}
		return coordinates_count == 1;
	}
	
	
	std::pair<BDouble *, BDouble *> generate_sherman_morrison_uv (const Matrix& system, int leech_y, int leech_x) {
		int columns = system.upper_bandwidth();
        int rows = system.rows() / columns;
        int limit = columns * rows;
		//Construimos el vector columna con un vector canonico	
		//especificando la fila que corresponde a la ecuacion
		//donde hay una sanguijuela
		BDouble* u = new BDouble[system.rows()];
		for (int ijEq = 0; ijEq < limit; ijEq++) {
			int i = ijEq / columns;
            int j = ijEq % columns;
			if(i == leech_y && j == leech_x) {
				u[ijEq] = 1.0;
			} else {
				u[ijEq] = 0.0;
			}
		}
		
		//Armamos el vector fila con un vector especificando
		//las columnas donde colocaremos las componentes
		//que corresponden a las diferencias finitas 		
		BDouble* v = new BDouble[system.rows()];
		for (int ijEq = 0; ijEq < limit; ijEq++) {
			u[ijEq] = 0.0;
		}
		int i = leech_y;
        int j = leech_x;
		u[(i * columns) + j - 1] = -0.25;
        u[(i * columns) + j + 1] = -0.25;
        u[((i - 1) * columns) + j] = -0.25;
        u[((i + 1) * columns) + j] = -0.25;
		
		return	std::pair<BDouble *, BDouble *>(u, v);
		
	}
	
	BDouble critic_point_temperature (Matrix& system, BDouble* solution) {
		int cpX = std::ceil((this->width / 2.0) / h); 
		int cpY = std::ceil((this->height / 2.0) /h );
		int ijEq = (cpY * system.upper_bandwidth()) + cpX;
		std::cout << "cpX: " << cpX << std::endl;
		std::cout << "cpY: " << cpY << std::endl;
		std::cout << "ijEq: " << ijEq << std::endl;
		std::cout << "solution[ijEq]: " << solution[ijEq] << std::endl;
		return solution[ijEq];
	}
	
	void clean_system (Matrix& system) {
		int columns = system.upper_bandwidth();
        int rows = system.rows() / columns;
        int limit = columns * rows;
		 for (int ijEq = 0; ijEq < limit; ijEq++) {
            system(ijEq, ijEq) = 0.0;
			
			int bound = std::min(system.upper_bandwidth(), system.lower_bandwidth());
			for(int h = 1; h <= bound; h++){
				if(ijEq > h){
					system(ijEq, ijEq-h) = 0.0;
				}				
				if(ijEq + h < limit) {
					system(ijEq, ijEq-h) = 0.0;
				}
			}
        }	
	}
	
	void discretize_leech (Leech leech, int& topX, int& bottomX, int& topY, int& bottomY) {
		/*int columns = system.upper_bandwidth();
        int rows = system.rows() / columns;
        int limit = columns * rows;

		// Distribuimos las temperaturas de la sanguijuela
        topX = std::floor((leech.x + leech.radio)/h);
        bottomX = std::ceil((leech.x - leech.radio)/h);
        topY = std::floor((leech.y + leech.radio)/h);
        bottomY = std::ceil((leech.y - leech.radio)/h);	
			
        // Ponemos las coordenadas en rango
        topX = std::min(std::max(topX, 0), columns - 1);
        bottomX = std::min(std::max(bottomX, 0), columns - 1);

        topY = std::min(std::max(topY, 0), rows - 1);
        bottomY = std::min(std::max(bottomY, 0), rows - 1);*/
				
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