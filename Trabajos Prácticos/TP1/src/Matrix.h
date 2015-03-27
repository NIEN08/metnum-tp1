//
// Created by julian on 3/21/15.
//

#ifndef _TP1_MATRIX_H_
#define _TP1_MATRIX_H_ 1

#include <iostream>
#include "BDouble.h"
#include <limits>
#include <cassert>
#include <utility>

class BandMatrix;
std::ostream &operator<<(std::ostream &os, const BandMatrix &m);

enum Solutions {
    INFINITE,
    SINGLE,
    NONE
};

/*
* Matriz Banda.
*/
class BandMatrix {
public:
    BandMatrix(std::size_t N, std::size_t M, std::size_t lband = SIZE_MAX, std::size_t uband = SIZE_MAX)
            : N(N), M(M), uband(uband), lband(lband) {
        assert(N > 0 && M > 0);

        if (lband == SIZE_MAX) {
            this->lband = N;
        }

        if (uband == SIZE_MAX) {
            this->uband = M;
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

    BandMatrix(const BandMatrix &m)
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

    std::size_t rows() const {
        return this->N;
    }

    std::size_t columns() const {
        return this->M;
    }

    std::size_t upper_bandwidth() const {
        return this->uband;
    }

    std::size_t lower_bandwidth() const {
        return this->lband;
    }

    BDouble &operator()(std::size_t i, std::size_t j) {
#ifdef DEBUG
        if (j < 0 || i < 0 || j >= this->columns() || i >= this->rows()) {
            throw new std::out_of_range("Index access out of range");
        }
        #endif

        if (i <= j + this->lower_bandwidth() && j <= i + this->upper_bandwidth()) {
            return matrix[i][j - i + this->lower_bandwidth()];
        } else {
            throw new std::out_of_range("Out of modifiable range");
        }
    }

    const BDouble &operator()(std::size_t i, std::size_t j) const {
#ifdef DEBUG
        if (j < 0 || i < 0 || j >= this->columns() || i >= this->rows()) {
            throw new std::out_of_range("Index access out of range");
        }
        #endif

        if (i > j + this->lower_bandwidth()) {
            return zero;
        } else if (j > i + this->upper_bandwidth()) {
            return zero;
        } else {
            return matrix[i][j - i + this->lower_bandwidth()];
        }
    }

    BandMatrix &operator=(const BandMatrix &m) {
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

    bool operator==(const BandMatrix &m) const {
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

    bool operator!=(const BandMatrix &m) const {
        return !(*this == m);
    }

    BandMatrix &operator+=(const BandMatrix &m) {
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

    BandMatrix &operator*=(const BDouble &c) {
        std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

        for (std::size_t i = 0; i < this->rows(); ++i) {
            for (std::size_t j = 0; j < bound; ++j) {
                this->matrix[i][j] *= c;
            }
        }

        return *this;
    }

    std::pair<BDouble *, enum Solutions> gaussian_elimination(BDouble *b) {
        BandMatrix workspace(*this);

        std::size_t diagonal = std::min(workspace.columns(), workspace.rows());

        for (std::size_t d = 0; d < diagonal; ++d) {
            std::size_t i = d + 1;

            if (workspace(d, d) == 0.0) {
                // Hay un cero en la base
                bool swap = false;

                for (i = d + 1; i < d + workspace.lower_bandwidth(); ++i) {
                    if (workspace(i, d) != 0.0) {
                        // Encontramos una fila más abajo que es distinta de 0
                        swap = true;
                        break;
                    }
                }

                if (swap) {
                    // Swappeamos las filas, sólo entre los elementos posibles.
                    // TODO: cuidado, esto no podría romper la estructura de banda?
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

                        for (std::size_t j = d + 1; j < i + this->upper_bandwidth(); ++j) {
                            // Realizamos la resta a toda la fila.
                            workspace(i, j) -=  coefficient * workspace(d, j);
                        }
                    }
                }

                ++d;
            }
        }

        std::cout << workspace;

        // Workspace esta triangulado, b siguió igual
        return backward_substitution(workspace, b);
    }


    ~BandMatrix() {
        for (std::size_t i = 0; i < this->rows(); ++i) {
            delete[] this->matrix[i];
        }

        delete[] this->matrix;
    }
private:
    // m tiene que estar triangulada
    std::pair<BDouble *, enum Solutions> backward_substitution(BandMatrix &m, BDouble *b) {
        BDouble *x = new BDouble[m.rows()];
        enum Solutions solution = SINGLE;

        for (std::size_t i = m.rows() - 1; i > 0; --i) {
            if (m(i, i) == 0.0) {
                solution = INFINITE;
            } else {
                std::cout << m;
                x[i] = b[i];

                for (std::size_t j = m.columns() - 1; j >= i; --j) {
                    x[i] -= m(i, j + 1) * x[i+1];
                }

                x[i] /= m(i, i);
            }
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

std::ostream &operator<<(std::ostream &os, const BandMatrix &m) {
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.columns(); ++j) {
            os << m(i, j) << ' ';
        }

        os << std::endl;
    }

    os << std::endl;

    return os;
}

BandMatrix operator+(const BandMatrix &m, const BandMatrix &n) {
    BandMatrix output(m);
    output += n;
    return output;
}

BandMatrix operator*(const BandMatrix &m, const BDouble &c) {
    BandMatrix output(m);
    output *= c;
    return output;
}

// TODO: todo esto
BandMatrix operator*(const BandMatrix &m, const BandMatrix &n) {
    throw new std::runtime_error("Must implement operator for Matrix instance");
}

#endif //_TP1_MATRIX_H_
