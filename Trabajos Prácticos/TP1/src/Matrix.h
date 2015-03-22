//
// Created by julian on 3/21/15.
//

#ifndef _TP1_MATRIX_H_
#define _TP1_MATRIX_H_

#include <cstdint>
#include <iostream>
#include <stdexcept>

/*
* La matriz se crea como cualquier otro objeto. N es la cantidad de filas, M la cantidad de columnas.
* Para acceder a un elemento de la matriz se usa, suponiendo que la matriz es A, (i, j) el índice: A(i, j).
* Cabe destacar que los algorítmos por defecto funcionan seguro, pero van a ser super ineficientes, ya que no van
* a explotar propiedades de la matriz o su representación.
*/
template <typename F>
class Matrix {
public:
    Matrix<F>(std::size_t N, std::size_t M) : N(N), M(M) { }

    // Constructor por copia
    virtual Matrix<F>(const Matrix<F> &) = 0;

    std::size_t rows() const {
        return N;
    };

    std::size_t columns() const {
        return M;
    };

    // Asignación
    Matrix<F> &operator=(const Matrix<F> &m) {
        if (this == &m) {
            return *this;
        } else if (this->rows() != m.rows() || this->columns() != m.columns()) {
            throw new std::out_of_range("Different row or column number");
        } else {
            for (auto i = 0; i < this->rows(); ++i) {
                for (auto j = 0; j < this->columns(); ++j) {
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

    // Impresión en pantalla
    std::ostream &operator<<(std::ostream &os, Matrix<F> const &m) {
        for (auto i = 0; i < this->rows(); ++i) {
            for (auto j = 0; j < this->columns(); ++j) {
                os << (*this)(i, j) << ' ';
            }

            os << '\n';
        }

        os << '\n';

        return os;
    }

    // Suma de matrices
    Matrix<F> &operator+=(const Matrix<F> &m) {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            throw new std::out_of_range("Different row or column number");
        } else {
            for (auto i = 0; i < this->rows(); ++i) {
                for (auto j = 0; j < this->columns(); ++j) {
                    (*this)(i, j) += m(i, j);
                }
            }

            return *this;
        }
    }

    const Matrix<F> operator+(const Matrix<F> &m) const {
        Matrix<F> output(*this);
        output += m;
        return output;
    };

    // Igualdad
    bool operator==(const Matrix<F> &m) const {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            return false;
        } else {
            for (auto i = 0; i < this->rows(); ++i) {
                for (auto j = 0; j < this->columns(); ++j) {
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
        for (auto i = 0; i < this->rows(); ++i) {
            for (auto j = 0; j < this->columns(); ++j) {
                (*this)(i, j) *= c;
            }
        }

        return *this;
    };

    const Matrix<F> operator*(const F &c) {
        Matrix<F> output(*this);
        output *= c;
        return output;
    };

    // Producto por otra matriz
    Matrix<F> &operator*=(const Matrix<F> &m) {
        // TODO: Ver como hacer esto bien. Tiene que ser destructivo pero hacer el producto igual... Si esto no sale, cambiar el de abajo.
    };

    const Matrix<F> operator*(const Matrix<F> &m) {
        Matrix<F> output(*this);
        output *= m;
        return output;
    };

    virtual ~Matrix() { }
private:
    size_t N;
    size_t M;
};

/*
* Matriz ineficiente. Todos los algorítmos son malísimos.
*/
template <typename F>
class InefficientMatrix : public Matrix<F> {
public:
    // Constructor
    InefficientMatrix(std::size_t N, std::size_t M) : Matrix(N, M) {
        matrix = new F[N][M];

        for (auto i = 0; i < N; ++i) {
            for (auto j = 0; j < M; ++j) {
                (*this)(i, j) = 0;
            }
        }
    }

    // Constructor por copia
    InefficientMatrix<F>(const InefficientMatrix<F> &m) : Matrix(m.rows(), m.columns()) {
        matrix = new F[this->rows()][this->columns()];

        for (auto i = 0; i < this->rows(); ++i) {
            for (auto j = 0; j < this->columns(); ++j) {
                (*this)(i, j) = m(i, j);
            }
        }
    }

    // Lector de indice
    F &operator()(std::size_t i, std::size_t j) {
        return matrix[i][j];
    }

    // Lector de indice constante
    const F &operator()(std::size_t i, std::size_t j) const {
        return matrix[i][j];
    }

    // Destructor
    virtual ~InefficientMatrix() {
        delete[] this->matrix;
    }
private:
    // Matrix
    F **matrix;
};

#endif //_TP1_MATRIX_H_
