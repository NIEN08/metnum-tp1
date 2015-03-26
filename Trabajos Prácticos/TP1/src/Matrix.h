//
// Created by julian on 3/21/15.
//

#ifndef _TP1_MATRIX_H_
#define _TP1_MATRIX_H_

#include <cstdint>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include "BDouble.h"

/*
* Matriz Banda.
*/
class BandMatrix {
public:
    // Constructor
    BandMatrix(std::size_t N, std::size_t M, std::size_t lband = -1, std::size_t uband = -1)
            : N(N), M(M), uband(uband), lband(lband) {
        assert(N > 0 && M > 0);

        if (lband == -1) {
            this->lband = N;
        }

        if (uband == -1) {
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

    // Constructor por copia
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

    // Destructor
    virtual ~BandMatrix() {
        for (std::size_t i = 0; i < this->rows(); ++i) {
            delete[] this->matrix[i];
        }

        delete[] this->matrix;
    }

    // Asignación
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


 // Igualdad
    bool operator==(const BandMatrix &m) const {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            return false;
        } else {
                std::size_t diagonal = std::min(this->rows(), this->columns());
                std::size_t lower = std::max(this->lower_bandwidth(), m.lower_bandwidth());
                std::size_t upper = std::max(this->upper_bandwidth(), m.upper_bandwidth());
            for (std::size_t d = 0; d < diagonal; ++d) {
                for (std::size_t j = d - lower; j < d + upper; ++j) {
                    if ((*this)(d, j) != m(d, j)) {
                        return false;
                    }
                }
            }
            return true;
            }
        }


    // Igualdad  
    // Esta version anda bien pero no se si es eficiente. 
    // Recorre todas las posiciones de la matriz y se fija si todos los elementos son iguales
    // si NO guardamos los 0 entonces tendria que servir este algoritmo.
    /*bool operator==(const BandMatrix &m) const {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            return false;
        } else {
            std::size_t diagonal = std::min(this->rows(), this->columns());
            std::size_t lower = std::max(this->lower_bandwidth(), m.lower_bandwidth());
            std::size_t upper = std::max(this->upper_bandwidth(), m.upper_bandwidth());


            int ancho = this->columns();
            int alto = this->rows();

            for (std::size_t i = 0; i < alto; i++) {
                for (std::size_t j = 0; j < ancho; j++) {
                    if ((*this)(i, j) != m(i, j)) {
                        return false;
                    }
                }
            }

            return true;
        }
    }*/

    // Desigualdad
    bool operator!=(const BandMatrix &m) const {
        return !(*this == m);
    }

    // Lector de indice
    BDouble &operator()(std::size_t i, std::size_t j) {
        assert(j >= 0 && j < this->columns());
        assert(i >= 0 && i < this->rows());

        if (i <= j + this->lower_bandwidth() && j <= i + this->upper_bandwidth()) {
            return matrix[i][j - i + this->lower_bandwidth()];
        } else {
            throw new std::out_of_range("Out of modifiable range");
        }
    }

    // Lector de indice constante
    const BDouble &operator()(std::size_t i, std::size_t j) const {
        assert(j >= 0 && j < this->columns());
        assert(i >= 0 && i < this->rows());

        if (i > j + this->lower_bandwidth()) {
            return zero;
        } else if (j > i + this->upper_bandwidth()) {
            return zero;
        } else {
            return matrix[i][j - i + this->lower_bandwidth()];
        }
    }

    // Suma de matrices
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

    // Producto por una constante
    BandMatrix &operator*=(const BDouble &c) {
        std::size_t bound = this->lower_bandwidth() + this->upper_bandwidth() + 1;

        for (std::size_t i = 0; i < this->rows(); ++i) {
            for (std::size_t j = 0; j < bound; ++j) {
                this->matrix[i][j] *= c;
            }
        }

        return *this;
    };

    // Producto de matrices 
    // ANDA MAL: CORREGIR
    BandMatrix &operator*(const BandMatrix &m) {
        BDouble tmp;
        if(this->columns() == m.rows()){
            for(int i = 0; i < this->rows(); i++){ // i se mueve por las filas i de la matriz A
                for(int j = 0; j < m.columns(); j++){ // j se mueve en las columnas j de la matriz B
                    for(int l = 0; l < this->columns(); l++){ // Ej: A[1][1]*B[1][1] + A[1][2]*B[2][1] + .. + A[1][l]*B[l][1]
                       tmp += this->matrix[i][l] * m.matrix[l][j];
                       //std::cout << this->matrix[i][l] << m.matrix[l][j] << tmp << std::endl;
                    }
                     this->matrix[i][j] = tmp; // una vez q ya se cuando da la sumatoria de la fila i por la columna j, pongo el valor en M[i][j] 
                     tmp = 0;       // acá guardo el valor del producto de matrices en el índice i,j
                }
            }

        } else{
            // No podemos hacer el producto de matrices
            throw new std::out_of_range("Different dimensions for matrix product");
        }
        return *this;
    }


    BDouble *triangulate(BDouble *b) {
        //BDouble *x = new BDouble[this->rows()]; // TODO: esto
        std::size_t min = (this->rows() <= this->columns())? this->rows():this->columns();
        std::size_t d = 0;

        while (d < min) {
            std::size_t i = d + 1;

            if ((*this)(d, d) == 0.0) {
                // Hay un cero en la base
                bool swap = false;

                for (i = d + 1; i < d + this->lower_bandwidth(); ++i) {
                    if ((*this)(i, d) != 0.0) {
                        // Encontramos una fila más abajo que es distinta de 0
                        swap = true;
                        break;
                    }
                }

                if (swap) {
                    // Swappeamos las filas, sólo entre los elementos posibles.
                    // TODO: cuidado, esto no podría romper la estructura de banda?
                    for (std::size_t j = d - lower_bandwidth(); j < i + this->upper_bandwidth(); ++j) {
                        BDouble tmp = (*this)(d, j);
                        (*this)(d, j) = (*this)(i, j);
                        (*this)(i, j) = tmp;
                    }

                    // Realizamos el mismo cambio en la solución del sistema
                    BDouble tmp = b[d];
                    b[d] = b[i];
                    b[i] = b[d];
                } else {
                    ++d;
                }
            } else {
                // Tenemos algo distinto de cero en la base
                for (i = d + 1; i < d + this->lower_bandwidth(); ++i) {
                    if ((*this)(i, d) != 0.0) {
                        // Tenemos algo distinto de cero en alguna fila más abajo
                        BDouble coefficient = (*this)(i, d)/(*this)(d, d);

                        // Setear esto en 0 debería reducir el error posible (por ejemplo, restando números muy chicos)
                        (*this)(i, d) = 0.0;

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


private:
    // Matrix
    std::size_t N;
    std::size_t M;
    std::size_t uband;
    std::size_t lband;
    BDouble **matrix;
};

// Impresión en pantalla
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

const BandMatrix operator+(const BandMatrix &m, const BandMatrix &n) {
    BandMatrix output(m);
    output += n;
    return output;
}

const BandMatrix operator*(const BandMatrix &m, const BDouble &c) {
    BandMatrix output(m);
    output *= c;
    return output;
}

// TODO: todo esto
const BandMatrix operator*(const BandMatrix &m, const BandMatrix &n) {
    throw new std::runtime_error("Must implement operator for Matrix instance");
}

#endif //_TP1_MATRIX_H_
