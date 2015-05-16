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
    double x;
    double y;
    double radio;
    double temperature;
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
            const double &width,
            const double &height,
            const double &h,
            std::list<Leech> &leeches)
            : width(width), height(height), h(h), leeches(leeches), method(method) {
    }

    Matrix run() {
        int rows = std::round(height / h) + 1;
        int columns = std::round(width / h) + 1;
        int dims = rows * columns;

        std::cerr << "[Problem] Rows: " << rows << std::endl;
        std::cerr << "[Problem] Columns: " << columns << std::endl;
        std::cerr << "[Problem] Dimension: " << dims << std::endl;
        std::cerr << "[Problem] Method: " << method << std::endl;
		
		Matrix system(dims, dims, columns, columns);
		double *b = new double[dims];
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
                sherman_morrison_solution(system, b, temperatures);
                break;
        }

        delete[] b;
        return temperatures;
    }
private:
    void load_temperature_matrix (double *x, Matrix &temperatures) {
        for (int i = 0; i < temperatures.rows(); ++i) {
            for (int j = 0; j < temperatures.columns(); ++j) {
                // Cargamos los datos en la matriz
                temperatures(i, j) = x[(i * temperatures.columns()) + j];
            }
        }
    }

    void band_gaussian_elimination(Matrix &system, double *b, Matrix &temperatures) {
        build_system(system, b, leeches);
		std::pair<double *, enum Solutions> solution = gaussian_elimination(system, b);	
		load_temperature_matrix(solution.first, temperatures);
		delete[] solution.first;
    }
	
	std::pair<double *, enum Solutions> lu_resolution (Matrix& L, Matrix& U, double *b) {
		//Resolvemos el sistema Ly = b
        std::pair<double *, enum Solutions> partialSolution = forward_substitution(L, b);

        //Resolvemos el sistema Ux = y
        std::pair<double *, enum Solutions> finalSolution = backward_substitution(U, partialSolution.first);
		delete[] partialSolution.first;
		return finalSolution;
	
	}

    void lu_factorization(Matrix &A, double *b, Matrix &temperatures) {
        build_system(A, b, this->leeches);

        // Sea A la matriz del sistema de ecuaciones,
        // factorizamos A = LU con L, U triangulares inferior/superior
        std::pair<Matrix, Matrix> factors = LU_factorization(A);
        Matrix& L = factors.first;
        Matrix& U = factors.second;
		std::pair<double *, enum Solutions> finalSolution = lu_resolution(L, U, b);
    
        //Cargamos la solucion en la matriz de temperaturas
        load_temperature_matrix(finalSolution.first, temperatures);

        // Liberamos la memoria que usamos.
        delete[] finalSolution.first;
    }

    std::map<std::pair<int, int>, double> generate_affected_positions(const std::list<Leech> &leeches) const {
        std::map<std::pair<int, int>, double> associations;

        for (auto &leech : leeches) {
            // Ponemos el rango que vamos a generar
            double topJ = std::min(leech.x + leech.radio, this->width - this->h)/h;
            double bottomJ = std::max(leech.x - leech.radio, this->h)/h;

            double topI = std::min(leech.y + leech.radio, this->height - this->h)/h;
            double bottomI = std::max(leech.y - leech.radio, this->h)/h;

            // Seteamos las temperaturas en la matriz.
            // Cabe destacar, la temperatura de cada sanguijuela es igual para todos los puntos que cubre.
            for (int i = std::ceil(bottomI); i <= std::floor(topI); ++i) {
                for (int j = std::ceil(bottomJ); j <= std::floor(topJ); ++j) {
                    double coef = std::pow(i*this->h - leech.y, 2) + std::pow(i*this->h - leech.x, 2);

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

        return associations;
    }

    /**
     * Construimos:
     * - system, la matriz de ecuaciones que representa la relación de las temperaturas.
     * - b, el vector de resultados que representa las condiciones del sistema.
     * */
    void build_system(Matrix &system, double *b, const std::list<Leech> &leeches) const {
        int columns = system.upper_bandwidth();
        int rows = system.rows() / columns;
        int limit = columns * rows;

        std::map<std::pair<int, int>, double> associations = generate_affected_positions(leeches);

        for (int ijEq = 0; ijEq < limit; ijEq++) {
            system(ijEq, ijEq) = 1.0;
            int i = ijEq / columns;
            int j = ijEq % columns;

            if (i == 0 || j == 0 || i == rows - 1 || j == columns - 1) {
                //Si esta en el borde el valor esta fijo en -100.0 y no hay que usar la ecuacion de laplace
                b[ijEq] = -100.0;
            } else {
                try {
                    // Si la posicion se encuentra en el radio de una sanguijuela la temperatura que afecta la posicion
                    // es la de la sanguijuela y no hay que usar la ecuacion de laplace
                    b[ijEq] = associations.at(std::pair<int, int>(i, j));
                } catch(...) {
                    // Finalmente si no es borde ni sanguijuela, hay que usar la ecuacion de laplace las posiciones de
                    // los bordes se ignoran porque figuran con -100.0 y fija el valor.
                    b[ijEq] = 	0.0;

                    // t[i-1][j] + t[i, j-1] - 4*t[i, j] + t[i+1, j] + t[i, j+1] = 0
                    // - t[i-1][j] - t[i, j-1] - t[i+1, j] - t[i, j+1] = 0 con t[i, j] = 0
                    system(ijEq, (i * columns) + j - 1) = -0.25;
                    system(ijEq, (i * columns) + j + 1) = -0.25;
                    system(ijEq, ((i - 1) * columns) + j) = -0.25;
                    system(ijEq, ((i + 1) * columns) + j) = -0.25;
                }
            }
        }
    }

    double critic_point_temperature(const Matrix &system, double *solution) {
        double centerJ = this->width / 2.0;
        double centerI = this->height / 2.0;

        double topJ = std::min(centerJ/this->h + 1.0, this->width/this->h - 1.0);
        double bottomJ = std::max(centerJ/this->h - 1, 1.0);

        double topI = std::min(centerI/this->h + 1.0, this->height/this->h - 1.0);
        double bottomI = std::max(centerI/this->h - 1.0, 1.0);

        double output = 0.0;
        double k = 0;

        for (int i = std::ceil(bottomI); i <= std::floor(topI); ++i) {
            for (int j = std::ceil(bottomJ); j <= std::floor(topJ); ++j) {
                output += solution[i*system.columns() + j];
                ++k;
            }
        }

        output /= k;

        return output;
    }

    void clean_system(Matrix &system) {
        int columns = system.upper_bandwidth();
        int rows = system.rows() / columns;
        int limit = columns * rows;
        int bound = std::min(system.upper_bandwidth(), system.lower_bandwidth());

        for (int ijEq = 0; ijEq < limit; ++ijEq) {
            system(ijEq, ijEq) = 0.0;

            for(int h = 1; h <= bound; ++h) {
                if (ijEq > this->h) {
                    system(ijEq, ijEq - this->h) = 0.0;
                }

                if (ijEq + this->h < limit) {
                    system(ijEq, ijEq + this->h) = 0.0;
                }
            }
        }
    }

	/**
	* Resuelve el problema por eliminacion gaussiana. En caso de que la temperatura del
	* punto critico sea mayor o igual a 235.0 grados de temperatura, resuelve el sistema
	* por cada sanguijuela, removiendo una de estas y se queda con la menor temperatura.
	**/
    void simple_algorithm(Matrix &system, double *b, Matrix &temperatures) {
		build_system(system, b, this->leeches);	
		
		//Solucion sin sacar sanguijuela
		std::pair<double *, enum Solutions> solution = gaussian_elimination(system, b);
		double* minX = solution.first;
		double minCp = critic_point_temperature(system, minX);
					
		if (minCp >= 235.0) {
			int leech_index = 0 ;
			for (std::list<Leech>::iterator itLeech = leeches.begin(); itLeech != leeches.end(); ++itLeech) {
				Leech leech = *itLeech;
			
				//Armamos una lista sin la sanguijuela
				std::list<Leech> leechesCopy = std::list<Leech>(leeches);
				std::list<Leech>::iterator itCopy = leechesCopy.begin();
				advance(itCopy, leech_index); 
				leechesCopy.erase(itCopy);
				
				std::cout << "====== REMOVED LEECH ======" << std::endl;
				std::cout << leech.x << " " << leech.y << " " << leech.radio << " " << leech.temperature << std::endl;			
				
				//Inicializamos el sistema sin la sanguijuela	
				double* b2 = new double[system.columns()];
				clean_system(system);
				//system = Matrix(system.columns(), system.columns(), system.upper_bandwidth(), system.upper_bandwidth());
				build_system(system, b2, leechesCopy);
							
				//Resolvemos el sistema
				solution = gaussian_elimination(system, b2);
				
				double cp = critic_point_temperature(system, solution.first);		
					
				//Liberamos memoria
				delete[] b2;
				
				// Nos quedamos con la solucion si es mejor que la anterior 
				if(cp < minCp) {
					minX = solution.first;
					minCp = cp;
				} else {
					delete[] solution.first;
				}			
				leech_index++;
			}
		}
		std::cout << "MIN CPT: " << minCp << std::endl;		
		load_temperature_matrix(minX, temperatures);	
    }

    /**
    * Devuelve true si la sanguijuela solo afecta una ecuacion del sistema.
    */
    bool is_singular_leech (const Leech &leech) {
        // Ponemos el rango que vamos a chequear
        double topJ = std::min(leech.x + leech.radio, this->width - this->h)/this->h;
        double bottomJ = std::max(leech.x - leech.radio, this->h)/this->h;

        double topI = std::min(leech.y + leech.radio, this->height - this->h)/this->h;
        double bottomI = std::max(leech.y - leech.radio, this->h)/this->h;

        int	count = 0;

        for (int i = std::ceil(bottomI); i <= std::floor(topI); ++i) {
            for (int j = std::ceil(bottomJ); j <= std::floor(topJ); ++j) {
                double coef = std::pow(i*this->h - leech.y, 2) + std::pow(j*this->h - leech.x, 2);

                if (coef <= std::pow(leech.radio, 2)) {
                    ++count;
                }
            }
        }

        return count == 1;
    }

	int singular_leeches_count() {
		int count = 0;

        for (auto &leech : leeches) {
            if (is_singular_leech(leech)) {
                ++count;
            }
        }

        return count;
	}
	
	std::pair<double *, enum Solutions> singular_leech_resolution(Matrix& system, Matrix& L, Matrix& U, double *b, std::list<Leech> &leeches, Leech removed_leech) {
		//Nos fijamos si otra sanguijuela afecta la posicion de esta
		int i = round(removed_leech.y / h);
		int j = round(removed_leech.x / h);
		
		//Tratamiento para sanguijuelas singulares (afectan una sola ecuacion)
		std::map<std::pair<int, int>, double> affected_positions = generate_affected_positions(leeches);	
		bool affected_position = affected_positions.count(std::pair<int, int>(i, j)) >= 1;
		
		std::pair<double *, enum Solutions> solution;

		if (affected_position) {
			//Otra sanguijuela afecta la posicion => No podemos aprovechar sherman-morrison.
			//Utilizamos unicamente la factorizacion LU.
			double newTemperature = affected_positions.at(std::pair<int, int>(i, j));
					
			// Inicializamos la solucion del sistema
			double* b2 = new double[system.columns()];
			std::copy(b, b + system.columns(), b2);
			//std::cout << "(i * system.upper_bandwidth()) + j == " << (i * system.upper_bandwidth()) + j << std::endl;
			//std::cout << "old temperature: " << b2[(i * system.upper_bandwidth()) + j]  << std::endl;
			//std::cout << "new temperature: " << newTemperature  << std::endl;
			b2[(i * system.columns()) + j] = newTemperature;
					
			// Resolvemos utilizando LU
			solution = lu_resolution(L, U, b2);
			delete[] b2;	
				
		} else {
			//Podemos aprovechar sherman-morrison!!
			std::pair<double *, double *> uv = generate_sherman_morrison_uv(system, i, j);
			double* u = uv.first;
			double* v = uv.second;
			double* b2 = generate_sherman_morrison_b(system, b, i, j);
					
			//Resolvemos utilizando sherman-morrison
			solution = sherman_morrison(L, U, u, v, b2);
					
			//Liberamos memoria
			delete[] u;
			delete[] v;
			delete[] b2;
					
		}

		return solution;
	}
	
	/**
	* En caso de que la cantidad de sanguijuelas singulares (afectan una sola ecuacion del sistema discretizado)
	* sea menor o igual a 1 resuelve el problema usando simple_algorithm. 
	* En caso contrario obtiene la factorizacion LU del sistema y separa el tratamiento de sanguijuelas normales
	* de las sanguijuelas singulares.
	* - Si la sanguijuela no es singular resuelve rehaciendo el sistema sin la sanguijuela como en simple_algorithm.
	* - Si la sanguijuela es singular a su vez separa en dos casos:
	* 	- Si la posicion se encuentra afectada por otra sanguijuela, simplemente modifica el valor del vector 
	*	correspondiente por el de mayor temperatura y resuelve utilizando la factorizacion LU.
	*	- Si la posicion no se encuentra afectada por otra sanguijuela, resuelve utilizando 
	**/
    void sherman_morrison_solution(Matrix &system, double *b, Matrix &temperatures) {
		std::list<Leech> singularLeeches;
		build_system(system, b, this->leeches);
		
		if (singular_leeches_count() < 2) {
			//Si la cantidad de sanguijuelas singulares es menor es 0 o 1
			//no tiene sentido obtener la factorizacion LU de la matriz.
			//Basta con utilizar la version simple del metodo
			simple_algorithm(system, b, temperatures);
			return;	
		}
		
		//Calculamos la factorizacion LU para aprovechar en las sanguijuelas singulares
		std::pair<Matrix, Matrix> factors = LU_factorization(system);
        Matrix& L = factors.first;
        Matrix& U = factors.second;			
	
		//Solucion sin sacar sanguijuela
		std::pair<double *, enum Solutions> solution = lu_resolution(L, U, b);
		double* minX = solution.first;
		double minCp = critic_point_temperature(system, minX);
			
		if (minCp >= 235.0) {
			int leech_index = 0 ;
			for (std::list<Leech>::iterator itLeech = leeches.begin(); itLeech != leeches.end(); ++itLeech) {
				Leech leech = *itLeech;
			
				//Armamos una lista sin la sanguijuela
				std::list<Leech> leechesCopy = std::list<Leech>(leeches);
				std::list<Leech>::iterator itCopy = leechesCopy.begin();
				advance(itCopy, leech_index); 
				leechesCopy.erase(itCopy);
				
				std::cout << "====== REMOVED LEECH ======" << std::endl;
				std::cout << leech.x << " " << leech.y << " " << leech.radio << " " << leech.temperature << std::endl;			
				
				if(is_singular_leech(leech) ){
					//Tratamos a las sanguijuelas singulares aparte
					solution = singular_leech_resolution(system, L, U, b, leechesCopy, leech);
					
				} else {				
					//Inicializamos el sistema sin la sanguijuela	
					double* b2 = new double[system.columns()];
					clean_system(system);
					build_system(system, b2, leechesCopy);
				
					//Resolvemos el sistema
					solution = gaussian_elimination(system, b2);
					
					//Liberamos memoria
					delete[] b2;
				}
				
				double cp = critic_point_temperature(system, solution.first);

				// Nos quedamos con la solucion si es mejor que la anterior 
				if(cp < minCp) {
					minX = solution.first;
					minCp = cp;
				} else {
					delete[] solution.first;
				}

				leech_index++;
			}
		}
		std::cout << "MIN CPT: " << minCp << std::endl;	
		load_temperature_matrix(minX, temperatures);		
	}
	
	double * generate_sherman_morrison_b (const Matrix& system, double* b, int leech_y, int leech_x) {
		int columns = system.upper_bandwidth();
		double* b2 = new double[system.columns()];
		std::copy(b, b + system.columns(), b2);
		//std::cout << "b2[" << (leech_y * columns) + leech_x << "] = " << b2[(leech_y * columns) + leech_x] << std::endl;
		b2[(leech_y * columns) + leech_x] = 0.0;	
		return b2;
	}
	
	
	std::pair<double *, double *> generate_sherman_morrison_uv (const Matrix& system, int leech_y, int leech_x) {
		int columns = system.upper_bandwidth();
        int rows = system.rows() / columns;
        int limit = columns * rows;
		
		//std::cout << "columns: " << columns << std::endl;
		//std::cout << "rows: " << rows << std::endl;
		//std::cout << "limit: " << limit << std::endl;
		
		//Construimos el vector columna con un vector canonico	
		//especificando la fila que corresponde a la ecuacion
		//donde hay una sanguijuela
		double* u = new double[system.rows()];
		for (int ijEq = 0; ijEq < limit; ijEq++) {
			int i = ijEq / columns;
            int j = ijEq % columns;
			if(i == leech_y && j == leech_x) {
				u[ijEq] = 1.0;
				//std::cout << "u[" << ijEq << "] = " << 1.0 << std::endl; 
			} else {
				u[ijEq] = 0.0;
			}
		}
		
		//Armamos el vector fila con un vector especificando
		//las columnas donde colocaremos las componentes
		//que corresponden a las diferencias finitas 		
		double* v = new double[system.rows()];
		for (int ijEq = 0; ijEq < limit; ijEq++) {
			v[ijEq] = 0.0;
		}
		int i = leech_y;
        int j = leech_x;
		v[(i * columns) + j - 1] = -0.25;
        v[(i * columns) + j + 1] = -0.25;
        v[((i - 1) * columns) + j] = -0.25;
        v[((i + 1) * columns) + j] = -0.25;
		//std::cout << "v[" << (i * columns) + j - 1 << "] = " << v[(i * columns) + j - 1] << std::endl; 
		//std::cout << "v[" << (i * columns) + j + 1 << "] = " << v[(i * columns) + j + 1] << std::endl; 
		//std::cout << "v[" << ((i - 1) * columns) + j<< "] = " << v[((i - 1) * columns) + j] << std::endl; 
		//std::cout << "v[" << ((i + 1) * columns) + j << "] = " << v[((i + 1) * columns) + j] << std::endl; 
		
		return	std::pair<double *, double *>(u, v);
		
	}

    double width;
    double height;
    double h;
    std::list<Leech> leeches;
    enum Method method;
};


#endif //_TP1_PROBLEM_H_