#ifndef _TP1_MATRIX_H_
#define _TP1_MATRIX_H_ 1

#include <stdexcept>
#include <cassert>
#include <utility>
#include <iostream>
#include "BDouble.h"

enum Solutions {
    INFINITE,
    SINGLE,
    NONE
};

// Este magic number nos dice cuándo convertir automáticamente una matriz banda en una matriz normal.
#define MAGIC_NUMBER 5698904538L

/*
* Matriz Banda.
*/
class Matrix {
    friend std::ostream &operator<<(std::ostream &, const Matrix &);
public:
    Matrix(const Matrix &m)
            : N(m.rows()), M(m.columns()), uband(m.upper_bandwidth()), lband(m.lower_bandwidth()) {
        std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;
        this->matrix = new BDouble*[this->rows()];

        for (std::size_t i = 0; i < this->rows(); ++i) {
            this->matrix[i] = new BDouble[bound];

            for (std::size_t j = 0; j < bound; ++j) {
                this->matrix[i][j] = m.matrix[i][j];
            }
        }
    }

    Matrix(std::size_t N, std::size_t M, std::size_t lband = MAGIC_NUMBER, std::size_t uband = MAGIC_NUMBER)
            : N(N), M(M), uband(uband), lband(lband) {
        if (this->rows() == 0 || this->columns() == 0) {
            throw new std::out_of_range("Invalid matrix dimension");
        }

        if (lband > N) {
            this->lband = N-1;
        }

        if (uband > M) {
            this->uband = M-1;
        }

        std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;
        this->matrix = new BDouble*[this->rows()];

        for (std::size_t i = 0; i < this->rows(); ++i) {
            this->matrix[i] = new BDouble[bound];

            for (std::size_t j = 0; j < bound; ++j) {
                this->matrix[i][j] = 0.0;
            }
        }
    }

    inline std::size_t rows() const {
        return this->N;
    }

    inline std::size_t columns() const {
        return this->M;
    }

    inline std::size_t upper_bandwidth() const {
        return this->uband;
    }

    inline std::size_t lower_bandwidth() const {
        return this->lband;
    }

    BDouble &operator()(const std::size_t &i, const std::size_t &j) {
        if (i >= this->rows() || j >= this->columns()) {
            throw new std::out_of_range("Index access out of range");
        }

        if (i <= j + this->lower_bandwidth() && j <= i + this->upper_bandwidth()) {
            return matrix[i][j - i + this->lower_bandwidth()];
        } else {
            throw new std::out_of_range("Out of modifiable range");
        }
    }

    const BDouble &operator()(const std::size_t &i, const std::size_t &j) const {
        if (i >= this->rows() || j >= this->columns()) {
            throw new std::out_of_range("Index access out of range");
        }

        if (i > j + this->lower_bandwidth()) {
            return zero;
        } else if (j > i + this->upper_bandwidth()) {
            return zero;
        } else {
            return matrix[i][j - i + this->lower_bandwidth()];
        }
    }

    Matrix &operator=(const Matrix &m) {
        if (*this != m) {
            // Limpiar memoria.
            for (std::size_t i = 0; i < this->rows(); ++i) {
                delete[] this->matrix[i];
            }

            delete[] this->matrix;

            // Poner información de la representación interna
            this->N = m.rows();
            this->M = m.columns();
            this->lband = m.lower_bandwidth();
            this->uband = m.upper_bandwidth();

            // Crear matriz nueva
            std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;
            this->matrix = new BDouble*[m.rows()];

            for (std::size_t i = 0; i < this->rows(); ++i) {
                this->matrix[i] = new BDouble[bound];

                for (std::size_t j = 0; j < bound; ++j) {
                    // Copiar los valores de la matriz
                    this->matrix[i][j] = m.matrix[i][j];
                }
            }
        }

        return *this;
    }

    bool operator==(const Matrix &m) const {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            return false;
        } else {
            for (std::size_t i = 0; i < this->rows(); i++) {
                for (std::size_t j = 0; j < this->columns(); j++) {
                    if ((*this)(i, j) != m(i, j)) {
                        return false;
                    }
                }
            }

            return true;
        }
    }

    bool operator!=(const Matrix &m) const {
        return !(*this == m);
    }

    Matrix &operator+=(const Matrix &m) {
        if (this->rows() == m.rows() && this->columns() == m.columns()) {
            // Si podemos sumar

            if (this->lower_bandwidth() == m.lower_bandwidth() && this->upper_bandwidth() == m.upper_bandwidth()) {
                // Si tenemos dos matrices banda con los mismos anchos de banda, simplemente sumamos la matriz miembro a miembro.
                std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

                for (std::size_t i = 0; i < this->rows(); ++i) {
                    for (std::size_t j = 0; j < bound; ++j) {
                        this->matrix[i][j] += m.matrix[i][j];
                    }
                }
            } else {
                // Si no, nos fijamos cuales son los nuevos anchos de banda
                std::size_t new_lband = std::max(this->lower_bandwidth(), m.lower_bandwidth());
                std::size_t new_uband = std::max(this->upper_bandwidth(), m.lower_bandwidth());
                std::size_t new_bound = new_lband + new_uband + 1;

                // Creamos una nueva matriz que guarda directamente la suma
                BDouble **output = new BDouble*[this->rows()];

                for (std::size_t i = 0; i < this->rows(); ++i) {
                    output[i] = new BDouble[new_bound];

                    for (std::size_t j = 0; j < new_bound; ++j) {
                        output[i][j] = this->matrix[i][j + i - this->lower_bandwidth()] + m.matrix[i][j + i - m.lower_bandwidth()];
                    }
                }

                // Nos aseguramos que los cambios sean efectivos
                this->lband = new_lband;
                this->uband = new_uband;

                // Borramos la matriz vieja
                for (std::size_t i = 0; i < this->rows(); ++i) {
                    delete[] this->matrix[i];
                }

                delete[] this->matrix;

                // Terminamos de fijar los cambios
                this->matrix = output;
            }
        } else {
            // No podemos sumar
            throw new std::out_of_range("Different dimensions for matrix sum");
        }

        return *this;
    }

    Matrix &operator*=(const BDouble &c) {
        std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

        for (std::size_t i = 0; i < this->rows(); ++i) {
            for (std::size_t j = 0; j < bound; ++j) {
                this->matrix[i][j] *= c;
            }
        }

        return *this;
    }

    std::pair<BDouble *, enum Solutions> gaussian_elimination(BDouble *b) {
        Matrix workspace(*this);

        std::size_t diagonal = std::min(workspace.columns(), workspace.rows());

        for (std::size_t d = 0; d < diagonal; ++d) {
            std::size_t i = d + 1;

            if (workspace(d, d) == 0.0) {
                // Hay un cero en la base
                bool swap = false;

                for (i = d + 1; i < std::min(d + workspace.lower_bandwidth(), workspace.rows()); ++i) {
                    if (workspace(i, d) != 0.0) {
                        // Encontramos una fila más abajo que es distinta de 0
                        swap = true;
                        break;
                    }
                }

                if (swap) {
                    // Swappeamos las filas, sólo entre los elementos posibles.
                    for (std::size_t j = d - workspace.lower_bandwidth(); j < i + workspace.upper_bandwidth(); ++j) {
                        BDouble tmp = workspace(d, j);
                        workspace(d, j) = workspace(i, j);
                        workspace(i, j) = tmp;
                    }

                    // Realizamos el mismo cambio en la solución del sistema
                    BDouble tmp = b[d];
                    b[d] = b[i];
                    b[i] = tmp;
                } else {
                    ++d;
                }
            } else {
                // Tenemos algo distinto de cero en la base
                for (i = d + 1; i < d + workspace.lower_bandwidth(); ++i) {
                    if (workspace(i, d) != 0.0) {
                        // Tenemos algo distinto de cero en alguna fila más abajo
                        BDouble coefficient = workspace(i, d)/workspace(d, d);

                        // Setear esto en 0 debería reducir el error posible (por ejemplo, restando números muy chicos)
                        workspace(i, d) = 0.0;

                        // Realizamos el mismo cambio en la solución del sistema
                        b[i] -= coefficient * b[d];

                        for (std::size_t j = d + 1; j < std::min(i + workspace.upper_bandwidth(), workspace.columns()); ++j) {
                            // Realizamos la resta a toda la fila.
                            workspace(i, j) -=  coefficient * workspace(d, j);
                        }
                    }
                }

                ++d;
            }
        }

        // Workspace esta triangulado, b siguió igual
        return backward_substitution(workspace, b);
    }


	std::pair<Matrix *, Matrix *> LU_factorization(BDouble *b) {
	
		//Original matrix
		Matrix* pA = this;
		std::size_t N = std::min(pA->rows(), pA->columns());

		//Lower and Upper triangular matrix
		Matrix* pL = new Matrix(pA->rows(), pA->columns(), pA->upper_bandwidth(), 0);
		Matrix* pU = new Matrix(pA->rows(), pA->columns(), 0, pA->lower_bandwidth());

		//Syntactic sugar on pointers
		Matrix& A = *pA;
		Matrix& L = *pL;
		Matrix& U = *pU;


		//Init L and U as identity matrixs
		for (std::size_t i = 0; i < std::min(A.rows(), A.columns()); i++) {
			L(i,i) = 1.0;
			U(i,i) = 1.0;
        }
		
		
		
		//Arbitraty choose that satisfy L(0,0) * U(0,0) = A(0,0)
		U(0,0) = A(0,0);	
		L(0,0) = 1.0;

		if ( U (0,0) == 0.0) {
			  // No podemos factorizar
            throw new std::out_of_range("Factorization impossible");
		}

		//Set first row of U and firt column of L
		for (std::size_t i = 1; i < N; i++) {
			U(0,i) = A(0,i) / L(0,0);
			L(i,0) = A(i,0) / U(0,0);
        }

		//Set rows/columns from 1 to n-1
		for (std::size_t i = 1; i < N - 1; i++) {
			
			U(i,i) = A(i,i);
			//Aprovechamos banda
			std::size_t bound = std::min(A.lower_bandwidth(), A.upper_bandwidth());
			for (std::size_t h = 1; h <= bound; h++) {
				if ( i >=h ) {
					U(i,i) -= L(i, i - h) * U(i - h,i);
				}
			}
			U(i,i) /= L(i,i);

			if ( U (i,i) == 0.0) {
				// No podemos factorizar
				throw new std::out_of_range("Factorization impossible");
			}
			
			
			for (std::size_t j = i+1; j < N; j++) {
				U(i,j) = A(i,j);
				L(j,i) = A(j,i);
				
				//Aprovechamos banda
				std::size_t bound = std::min(A.lower_bandwidth(), A.upper_bandwidth());
				for (std::size_t h = 1; h <= bound; h++) {
					if ( i >= h) {
						U(i,j) -= L(i,i-h) * U(i-h,j); // iº ROW OF U
						L(j,i) -= L(j,i-h) * U(i-h,i); // jº COLUMN OF L
					}
				}
				U(i,j) /= L(i,i);
				L(j,i) /= U(i,i);

			}

		}
		
		//Set last position
		U(N-1,N-1) = A(N-1,N-1);
		
		std::size_t bound = std::min(A.lower_bandwidth(), A.upper_bandwidth());
		for (std::size_t h = 1; h <= bound; h++) {
			U(N-1,N-1) -= L(N-1,N-1-h) * U(N-1-h,N-1);
		}
		U(N-1,N-1) /= L(N-1,N-1);

		return std::pair<Matrix*, Matrix*>(pL, pU);
	}
	
	/**
	* - L matriz triangular inferior de la descomposicion.
	* - U matriz triangular superior de la descomposicion.
	* - i fila del elemento a modificar.
	* - j columna del elemento a modificar.
	* - a nuevo valor de la posicion (i, j)
	* A matrix original del sistema.
	**/
	std::pair<BDouble *, enum Solutions> sherman_morrison(Matrix* pL, 
												Matrix* pU, 
												std::size_t i, 
												std::size_t j,
												BDouble a,
												BDouble *b) {
	
		//Original matrix
		Matrix* pA = this;
		std::size_t N = std::min(pA->rows(), pA->columns());

		//Syntactic sugar on pointers
		Matrix& A = *pA;
		Matrix& L = *pL;
		Matrix& U = *pU;
		
		//Sherman-Morrison formula vectors
		//Altered system: A2 = (A + uv')
		BDouble* u = new BDouble[N];
		BDouble* v = new BDouble[N];
		
		for (std::size_t k = 0; k < N; k++) {
			u[k] = 0.0;
			v[k] = 0.0;
		}
		//Column vector
		u[i] = 1.0;
		
		//Row vector
		v[j] = a - A(i,j);
		
		//From Sherman-Morrison
		// A^-1 b = y <=> Ay = b
		// A^-1 u = z <=> Az = u
		
		//First we solve:
		// L y2 = b and L z2 = u
		BDouble* y2;
		BDouble* z2;
		
		std::pair<BDouble *, enum Solutions> solution; 
		solution = L.forward_substitution(L, b);
		y2 = solution.first;
		solution = L.forward_substitution(L, u);
		z2 = solution.first;
		
		//Then we solve:
		// U y = y2 and U z = z2
		BDouble* y;
		BDouble* z;
		solution = L.backward_substitution(U, y2);
		y = solution.first;
		solution = L.backward_substitution(U, z2);
		z = solution.first;
		delete[] y2;
		delete[] z2;
		
		//Finally x = y - z * [(v' y)/(1 + v' z)]
		BDouble* x = new BDouble[N];
		
		//First we calculate k = (v' y)/(1 + v' z) (scalar value)
		BDouble vy = 0.0;
		BDouble vz = 1.0;
		for (std::size_t h = 0; h < N; h++) {
			vy += v[h] * y[h];
			vz += v[h] * z[h];
        }
		
		//Finally we calculate x = y + z * k
		for (std::size_t h = 0; h < N; h++) {
			x[h] = y[h] - (z[h] * (vy / vz));
        }
		delete[] y;
		delete[] z;
		
		return  std::pair<BDouble *, enum Solutions>(x, SINGLE);
	}
	
	// m tiene que estar triangulada
    // el usuario libera la memoria
    std::pair<BDouble *, enum Solutions> forward_substitution(BDouble *b) {
		return forward_substitution(*this, b); 
    }
	
	std::pair<BDouble *, enum Solutions> backward_substitution(BDouble *b) {
		return backward_substitution(*this, b); 
	}


    ~Matrix() {
        for (std::size_t i = 0; i < this->rows(); ++i) {
            delete[] this->matrix[i];
        }

        delete[] this->matrix;
    }
private:
    // m tiene que estar triangulada
    // el usuario libera la memoria
    std::pair<BDouble *, enum Solutions> backward_substitution(Matrix &m, BDouble *b) {
        // TODO: caso en el que no hay solución. No queda claro cómo detectarlo.
        BDouble *x = new BDouble[m.columns()];
        enum Solutions solution = SINGLE;

		std::size_t N = std::min(m.rows(), m.columns());
        std::size_t i = N - 1;

        while (true) {
            if (m(i, i) == 0.0) {
                x[i] = 1.0;
                solution = INFINITE;
            } else {
                std::size_t bound = std::min(m.columns(), m.upper_bandwidth());
                x[i] = b[i];
				
				//std::cout << "x" << i << "= [b" << i;
                for (std::size_t h = 1; h <= bound; h++) {
                    if(i + h < N) {
						//std::cout << "- M(" << i << ", " << i+h << ") * x" << i+h;
                        x[i] -= m(i, i+h) * x[i+h];
                    }
                }

                x[i] /= m(i, i);
                //std::cout << "] / M(" << i << ", " << i << ") " << std::endl; 
                
            }

            // When i = 0, decreasing i will land it to MAX_SIZE, which is higher than 0, producing an error.
            if (i == 0) {
                break;
            }

            --i;
        }

        return std::pair<BDouble *, enum Solutions>(x, solution);
    }
    
    // m tiene que estar triangulada
    // el usuario libera la memoria
    std::pair<BDouble *, enum Solutions> forward_substitution(Matrix &m, BDouble *b) {
        // TODO: caso en el que no hay solución. No queda claro cómo detectarlo.
        BDouble *x = new BDouble[m.columns()];
        enum Solutions solution = SINGLE;

		std::size_t N = std::min(m.rows(), m.columns());

        for (std::size_t i = 0; i < N; ++i) {
            if (m(i, i) == 0.0) {
                x[i] = 1.0;
                solution = INFINITE;
            } else {
                std::size_t bound = std::min(i, m.lower_bandwidth());
                //std::cout << "bound: " << bound << std::endl;
                x[i] = b[i];
                //std::cout << "x" << i << "= [b" << i;
                for (std::size_t h = 1; h <= bound; h++) {
                    //if (m(i, i+h) != 0.0) {
                    if(i >= h) {
						//std::cout << "- M(" << i << ", " << i-h << ") * x" << i-h;
						x[i] -= m(i, i-h) * x[i-h];
					}
                        //x[i] -= m(i, i+h) * x[i+h];
                    //}
                }

                x[i] /= m(i, i);
                //std::cout << "] / M(" << i << ", " << i << ") " << std::endl; 
            }

            // When i = 0, decreasing i will land it to MAX_SIZE, which is higher than 0, producing an error.
            //if (i == 0) {
            //    break;
            //}
        }

        return std::pair<BDouble *, enum Solutions>(x, solution);
    }


    // Matrix
    std::size_t N;
    std::size_t M;
    std::size_t uband;
    std::size_t lband;
    BDouble **matrix;
};

std::ostream &operator<<(std::ostream &os, const Matrix &m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.columns(); ++j) {
            if (i <= j + m.lower_bandwidth() && j <= i + m.upper_bandwidth()) {
                os << "\033[1;31m" << m(i, j) << "\033[0m" << '\t';
            } else {
                os << m(i, j) << '\t';
            }

        }

        os << std::endl;
    }

    os << std::endl;

    return os;
}

Matrix operator+(const Matrix &m, const Matrix &n) {
    Matrix output(m);
    output += n;
    return output;
}

Matrix operator*(const Matrix &m, const BDouble &c) {
    Matrix output(m);
    output *= c;
    return output;
}

// TODO: hay que revisar absolutamente TODOS los FOR que vayan por la diagonal:
// TODO: restar puede llevar a irse a la mierda con el indice por ir a valores "negativos"
// TODO: sumar puede terminar por irse a valores positivos chicos.
// TODO: esto no funciona.
Matrix operator*(const Matrix &m, const Matrix &n) {
    if (m.columns() == n.rows()) {
        std::size_t max_lower_upper = std::max(m.lower_bandwidth(), n.upper_bandwidth());
        std::size_t max_upper_lower = std::max(m.upper_bandwidth(), n.lower_bandwidth());

        Matrix output(m.rows(), n.columns(), max_lower_upper, max_upper_lower);

        std::size_t diagonal = std::min(output.columns(), output.rows());

        for (std::size_t d = 0; d < diagonal; ++d) {
            for (std::size_t j = d - output.lower_bandwidth(); j < std::min(d + output.upper_bandwidth(), output.columns()); ++j) {

                for (std::size_t k = 0; k < output.columns();  ++k) {
                    output(d, j) += m(d, k) * n(k, j);
                }
            }
        }

        return output;
    } else {
        throw new std::out_of_range("Matrix product between incompatible matrices.");
    }
}

#endif //_TP1_MATRIX_H_
