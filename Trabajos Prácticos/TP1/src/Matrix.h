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
* Matriz Banda.
*/
template <class F, F zero>
class BandMatrix {
public:
    // Constructor
    BandMatrix(std::size_t lband, std::size_t uband, std::size_t N, std::size_t M)
            : N(N), M(M), lband(lband), uband(uband) {
        assert(N > 0 && M > 0);
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
            : N(m.rows()), M(m.columns()), lband(m.lower_bandwidth()), uband(m.upper_bandwidth()) {
        std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;
        this->matrix = new F*[this->rows()];

        for (std::size_t i = 0; i < this->rows(); ++i) {
            this->matrix[i] = new F[bound];

            for (std::size_t j = 0; j < bound; ++j) {
                this->matrix[i][j] = m.matrix[i][j];
            }
        }
    }

    std::size_t rows() const {
        return this->N;
    };

    std::size_t columns() const {
        return this->M;
    }

    std::size_t upper_bandwidth() const {
        return this->uband;
    }

    std::size_t lower_bandwidth() const {
        return this->lband;
    }

    // Asignación
    BandMatrix<F, zero> &operator=(const BandMatrix<F, zero> &m) {
        if (this->rows() == m.rows() && this->columns() == m.columns()) {
            std::size_t k = std::min(this->rows(), this->columns());

            if (this->upper_bandwidth() == m.upper_bandwidth() && this->lower_bandwidth() == m.lower_bandwidth()) {
                for (std::size_t d = 0; d < k; ++d) {
                    for (std::size_t j = d - this->lower_bandwidth(); j < d + this->upper_bandwidth(); ++j) {
                        this->matrix[d][j] = m.matrix[d][j];
                    }
                }
            } else {
                std::size_t
                std::size_t bound = std::max(this->lower_bandwidth(), m.lower_bandwidth()) + std::max(this->upper_bandwidth(), m.upper_bandwidth()) + 1;
                F **output = new F*[this->rows()];

                // TODO: terminar asignación
            }
        } else {
            throw new std::out_of_range("Different row or column number");
        }
    }

    // Igualdad
    bool operator==(const BandMatrix<F, zero> &m) const {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            return false;
        } else {
            std::size_t diagonal = std::min(this->rows(), this->columns());
            std::size_t lower = std::max(this->lower_bandwidth(), m.lower_bandwidth());
            std::size_t upper = std::max(this->upper_bandwidth(), m.upper_bandwidth());

            for (std::size_t d = 0; d < diagonal; ++d) {
                for (std::size_t j = d - lower; j < d + upper; ++j) {
                    if ((*this)(i, j) != m(i, j)) {
                        return false;
                    }
                }
            }

            return true;
        }
    };

    // Desigualdad
    bool operator!=(const BandMatrix<F, zero> &m) const {
        return !(*this == m);
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
            return matrix[i][j - i + this->lower_bandwidth()]; // TODO: +1. Testear extensivamente.
        }
    }

    // Lector de indice constante
    const F &operator()(std::size_t i, std::size_t j) const {
        return this->operator()(i, j);
    }

    // Suma de matrices
    BandMatrix<F, zero> &operator+=(const BandMatrix<F, zero> &m) {
        if (this->rows() == m.rows() && this->columns() == m.columns()) {
            // Si podemos sumar

            if (this->lower_bandwidth() == m.lower_bandwidth() && this->upper_bandwidth() == m.upper_bandwidth()) {
                // Si tenemos dos matrices banda con los mismos anchos de banda, simplemente sumamos la matriz miembro a miembro.
                std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

                // TODO: cuidado, hay posiciones que van a ser 0 y podríamos estar accediendolas.
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
                F **output = new F*[this->rows()];

                for (std::size_t i = 0; i < this->rows(); ++i) {
                    output[i] = new F[new_bound];

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

            return *this;
        } else {
            // No podemos sumar
            throw new std::out_of_range("Different dimensions for matrix sum");
        }
    }

    // Producto por una constante
    Matrix<F> &operator*=(const F &c) {
        std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

        for (std::size_t i = 0; i < this->rows(); ++i) {
            for (std::size_t j = 0; j < bound; ++j) {
                this->matrix[i][j] *= c;
            }
        }

        return *this;
    };

    // Destructor
    virtual ~BandMatrix() {
        for (std::size_t i = 0; i < this->rows(); ++i) {
            delete[] this->matrix[i];
        }

        delete[] this->matrix;
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

    // Descomposición LU para matrices bandas
    std::tuple<BandMatrix<F, zero>, BandMatrix<F, zero>, BandMatrix<F, zero>>
        PLUDecomposition(BandMatrix<F> &m, int M, double b[M]) {

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
    std::size_t N;
    std::size_t M;
    std::size_t uband;
    std::size_t lband;
    F **matrix;
};

// Impresión en pantalla
template <class F, F zero>
std::ostream &operator<<(std::ostream &os, const BandMatrix<F, zero> &m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.columns(); ++j) {
            os << m(i, j) << ' ';
        }

        os << std::endl;
    }

    os << std::endl;

    return os;
}

// TODO: todo esto
template <class F, F zero>
const BandMatrix<F, zero> operator+(const BandMatrix<F, zero> &m, const BandMatrix<F, zero> &n) {
    Matrix<F> output(m);
    output += n;
    return output;
}

template <class F, F zero>
const BandMatrix<F, zero> operator*(const BandMatrix<F, zero> &m, const F &c) {
    throw new std::runtime_error("Must implement operator for Matrix instance");
}

template <class F, F zero>
const BandMatrix<F, zero> operator*(const BandMatrix<F, zero> &m, const BandMatrix<F, zero> &n) {
    throw new std::runtime_error("Must implement operator for Matrix instance");
}

#endif //_TP1_MATRIX_H_
