//
// Created by julian on 3/21/15.
//

#ifndef _TP1_MATRIX_H_
#define _TP1_MATRIX_H_

#include <cstdint>
#include <cassert>
#include <iostream>
#include <stdexcept>

/*
* La matriz se crea como cualquier otro objeto. N es la cantidad de filas, M la cantidad de columnas.
* Para acceder a un elemento de la matriz se usa, suponiendo que la matriz es A, (i, j) el índice: A(i, j).
* Cabe destacar que los algorítmos por defecto funcionan seguro, pero van a ser super ineficientes, ya que no van
* a explotar propiedades de la matriz o su representación.
*/
template <class F>
class Matrix {
public:
    Matrix<F>(std::size_t N, std::size_t M) : N(N), M(M) {
        assert(N > 0);
        assert(M > 0);
    }

    std::size_t rows() const {
        return this->N;
    };

    std::size_t columns() const {
        return this->M;
    };

    // Asignación
    Matrix<F> &operator=(const Matrix<F> &m) {
        if (this == &m) {
            return *this;
        } else if (this->rows() != m.rows() || this->columns() != m.columns()) {
            throw new std::out_of_range("Different row or column number");
        } else {
            for (std::size_t i = 0; i < this->rows(); ++i) {
                for (std::size_t j = 0; j < this->columns(); ++j) {
                    (*this)(i, j) = m(i, j);
                }
            }

            return *this;
        }
    }

    // Asignación a miembro
    virtual F &operator()(std::size_t, std::size_t) = 0;

    // Lectura de miembro
    virtual const F &operator()(std::size_t, std::size_t) const = 0;

    // Suma de matrices
    Matrix<F> &operator+=(const Matrix<F> &m) {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            throw new std::out_of_range("Different row or column number");
        } else {
            for (std::size_t i = 0; i < this->rows(); ++i) {
                for (std::size_t j = 0; j < this->columns(); ++j) {
                    (*this)(i, j) += m(i, j);
                }
            }

            return *this;
        }
    }

    // Igualdad
    bool operator==(const Matrix<F> &m) const {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            return false;
        } else {
            for (std::size_t i = 0; i < this->rows(); ++i) {
                for (std::size_t j = 0; j < this->columns(); ++j) {
                    if ((*this)(i, j) != m(i, j)) {
                        return false;
                    }
                }
            }

            return true;
        }
    };

    bool operator!=(const Matrix<F> &m) const {
        return !(*this == m);
    }

    // Producto por una constante
    Matrix<F> &operator*=(const F &c) {
        for (std::size_t i = 0; i < this->rows(); ++i) {
            for (std::size_t j = 0; j < this->columns(); ++j) {
                (*this)(i, j) *= c;
            }
        }

        return *this;
    };

    // Producto por otra matriz
    Matrix<F> &operator*=(const Matrix<F> &m) {
        // TODO: Ver como hacer esto bien. Tiene que ser destructivo pero hacer el producto igual... Si esto no sale, cambiar el de abajo.
    };

    virtual ~Matrix() { }
private:
    size_t N;
    size_t M;
};

// Impresión en pantalla
template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.columns(); ++j) {
            os << m(i, j) << ' ';
        }

        os << std::endl;
    }

    os << std::endl;

    return os;
}

/*
* Matriz Banda.
*/
template <class F, F zero>
class BandMatrix : public Matrix<F> {
public:
    // Constructor
    BandMatrix(std::size_t lband, std::size_t uband, std::size_t N, std::size_t M)
            : Matrix<F>(N, M), lband(lband), uband(uband) {
        assert(lband < N);
        assert(uband < M);

        std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

        this->matrix = new F*[this->rows()];

        for (std::size_t i = 0; i < this->rows(); ++i) {
            this->matrix[i] = new F[bound];
            for (std::size_t j = 0; j < bound; ++j) {
                this->matrix[i][j] = zero;
            }
        }
    }

    // Constructor por copia
    BandMatrix<F, zero>(const BandMatrix<F, zero> &m)
            : Matrix<F>(m.rows(), m.columns()), lband(m.lower_bandwidth()), uband(m.upper_bandwidth()) {
        std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

        this->matrix = new F*[this->rows()];

        for (std::size_t i = 0; i < this->rows(); ++i) {
            this->matrix[i] = new F[bound];

            for (std::size_t j = 0; j < bound; ++j) {
                this->matrix[i][j] = m.matrix[i][j];
            }
        }
    }

    std::size_t upper_bandwidth() const {
        return uband;
    }

    std::size_t lower_bandwidth() const {
        return lband;
    }

    // Lector de indice
    F &operator()(std::size_t i, std::size_t j) {
        assert(j > 0 && j < this->columns());
        assert(i > 0 && i < this->rows());

        if (i > j + this->lower_bandwidth()) {
            return zero;
        } else if (j > i + this->upper_bandwidth()) {
            return zero;
        } else {
            return matrix[i][j - i + this->lower_bandwidth()]; // +1
        }
    }

    // Lector de indice constante
    const F &operator()(std::size_t i, std::size_t j) const {
        assert(j > 0 && j < this->columns());
        assert(i > 0 && i < this->rows());

        if (i > j + this->lower_bandwidth()) {
            return zero;
        } else if (j > i + this->upper_bandwidth()) {
            return zero;
        } else {
            return matrix[i][j - i + this->lower_bandwidth()]; // +1
        }
    }

    F *triangulate(F *b) {
        F *x = new F[this->rows()];
        std::size_t min = (this->rows() <= this->columns())? this->rows():this->columns();
        std::size_t d = 0;

        while (d < min) {
            std::size_t i = d + 1;

            if ((*this)(d, d) == 0) {
                // Hay un cero en la base
                bool swap = false;

                for (i = d + 1; i < d + this->lower_bandwidth(); ++i) {
                    if ((*this)(i, d) != 0) {
                        // Encontramos una fila más abajo que es distinta de 0
                        swap = true;
                        break;
                    }
                }

                if (swap) {
                    // Swappeamos las filas, sólo entre los elementos posibles.
                    // TODO: cuidado, esto no podría romper la estructura de banda?
                    for (std::size_t j = d - lower_bandwidth(); j < i + this->upper_bandwidth(); ++j) {
                        F tmp = (*this)(d, j);
                        (*this)(d, j) = (*this)(i, j);
                        (*this)(i, j) = tmp;
                    }

                    // Realizamos el mismo cambio en la solución del sistema
                    F tmp = b[d];
                    b[d] = b[i];
                    b[i] = b[d];
                } else {
                    ++d;
                }
            } else {
                // Tenemos algo distinto de cero en la base
                for (i = d + 1; i < d + this->lower_bandwidth(); ++i) {
                    if ((*this)(i, d) != 0) {
                        // Tenemos algo distinto de cero en alguna fila más abajo
                        F coefficient = (*this)(i, d)/(*this)(d, d);

                        // Setear esto en 0 debería reducir el error posible (por ejemplo, restando números muy chicos)
                        (*this)(i, d) = 0;

                        // Realizamos el mismo cambio en la solución del sistema
                        b[i] -= coefficient * b[d];

                        for (std::size_t j = d + 1; j < i + this->upper_bandwidth(); ++j) {
                            // Realizamos la resta a toda la fila.
                            (*this)(i, j) -=  coefficient * (*this)(d, j);
                        }
                    }
                }

                ++d;
            }
        }
    }

    // Destructor
    virtual ~BandMatrix() {
        for (std::size_t i = 0; i < this->rows(); ++i) {
            delete[] this->matrix[i];
        }

        delete[] this->matrix;
    }

    // Descomposición LU para matrices bandas
    BandMatrix<F> LUDecomposition( BandMatrix<F> &m, int M, double b[M] ){

        // Primero no habria que aplicar eliminacion gaussiana para despues encontrar las soluciones?
        int i = M;
        int j = M;
        int y = b[i] / m[i][i];
        int x[M];
        int x[i] = b[i] / m[i][i];

        for(int i = M - 1; 0 < i; i--){
            x[i] = (b[i] - (y * m[i][j])) / m[i - 1][j];
            j--;
            y = y + (m[i - 1][j] * x[i]);
        }
        return m;
    }


private:
    // Matrix
    std::size_t uband;
    std::size_t lband;
    F **matrix;
};


template <class F>
const Matrix<F> operator+(const Matrix<F> &m, const Matrix<F> &n) {
    Matrix<F> output(m);
    output += n;
    return output;
}

template <class F>
const Matrix<F> operator*(const Matrix<F> &m, const F &c) {
    throw new std::runtime_error("Must implement operator for Matrix instance");
}

template <class F>
const Matrix<F> operator*(const Matrix<F> &m, const Matrix<F> &n) {
    throw new std::runtime_error("Must implement operator for Matrix instance");
}

#endif //_TP1_MATRIX_H_
