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
* Matriz ineficiente. Todos los algorítmos son malísimos.
*/
template <class F, F zero>
class InefficientMatrix : public Matrix<F> {
public:
    // Constructor
    InefficientMatrix(std::size_t N, std::size_t M) : Matrix<F>(N, M) {
        this->matrix = new F*[N];

        for (std::size_t i = 0; i < N; ++i) {
            this->matrix[i] = new F[M];

            for (std::size_t j = 0; j < M; ++j) {
                this->matrix[i][j] = zero;
            }
        }
    }

    // Constructor por copia
    InefficientMatrix<F, zero>(const InefficientMatrix<F, zero> &m) : Matrix<F>(m.rows(), m.columns()) {
        this->matrix = new F*[this->rows()];

        for (std::size_t i = 0; i < this->rows(); ++i) {
            this->matrix[i] = new F[this->columns()];

            for (std::size_t j = 0; j < this->columns(); ++j) {
                this->matrix[i][j] = m(i, j);
            }
        }
    }

    // Lector de indice
    F &operator()(std::size_t i, std::size_t j) {
        assert(j >= 0 && j < this->columns());
        assert(i >= 0 && i < this->rows());

        return this->matrix[i][j];
    }

    // Lector de indice constante
    const F &operator()(std::size_t i, std::size_t j) const {
        assert(j >= 0 && j < this->columns());
        assert(i >= 0 && i < this->rows());

        return this->matrix[i][j];
    }

    // TODO: lograr borrar estos dos metodos.
    virtual const InefficientMatrix<F, zero> operator+(const InefficientMatrix<F, zero> &m) const {
        InefficientMatrix<F, zero> output(*this);
        output += m;
        return output;
    };

    const InefficientMatrix<F, zero> operator*(const F &c) {
        InefficientMatrix<F, zero> output(*this);
        output *= c;
        return output;
    };

    // Destructor
    virtual ~InefficientMatrix() {
        for (std::size_t i = 0; i < this->rows(); ++i) {
            delete[] this->matrix[i];
        }

        delete[] this->matrix;
    }
private:
    // Matrix
    F **matrix;
};

/*
* Matriz Banda.
*/
template <class F, F zero>
class BandMatrix : public Matrix<F> {
public:
    // Constructor
    /*
Formally, consider an n×n matrix A=(ai,j ). If all matrix elements are zero outside a diagonally bordered band whose range
    is determined by constants k1 and k2:

a_{i,j}=0 \quad\mbox{if}\quad j<i-k_1 \quad\mbox{ or }\quad j>i+k_2; \quad k_1, k_2 \ge 0.\,
then the quantities k1 and k2 are called the lower and upper bandwidth, respectively.[1] The bandwidth of the matrix is
    the maximum of k1 and k2; in other words, it is the number k such that  a_{i,j}=0  if  |i-j| > k .[2]

A matrix is called a band matrix or banded matrix if its bandwidth is reasonably small.

A band matrix with k1 = k2 = 0 is a diagonal matrix; a band matrix with k1 = k2 = 1 is a tridiagonal matrix; when
k1 = k2 = 2 one has a pentadiagonal matrix and so on. If one puts k1 = 0, k2 = n−1, one obtains the definition of an
upper triangular matrix; similarly, for k1 = n−1, k2 = 0 one obtains a lower triangular matrix.
    */
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
    // TODO: esto SEGURO que está mal. Si te fuiste de rango de lo que "está definido", quiero que recibas una referencia que no haga nada.
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
private:
    // Matrix
    std::size_t uband;
    std::size_t lband;
    F **matrix;
};

#endif //_TP1_MATRIX_H_
