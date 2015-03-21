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
*/
template <typename T>
class Matrix {
public:
    Matrix<T>(std::size_t N, std::size_t M) : N(N), M(M) { }

    // Constructor por copia
    virtual Matrix<T>(const Matrix<T> &) = 0;

    std::size_t rows() const {
        return N;
    };

    std::size_t columns() const {
        return M;
    };

    // Asignación a miembro
    virtual T &operator()(std::size_t, std::size_t) = 0;
    virtual const T &operator()(std::size_t, std::size_t) const = 0;

    // Printing
    virtual std::ostream &operator<<(std::ostream &, Matrix<T> const &) = 0;

    // Asignación
    virtual Matrix<T> &operator=(const Matrix<T> &) = 0;

    // Igualdad
    virtual bool operator==(const Matrix<T> &) const = 0;
    virtual bool operator!=(const Matrix<T> &) const = 0;

    // Suma de matrices
    virtual Matrix<T> &operator+=(const Matrix<T> &) = 0;
    virtual const Matrix<T> operator+(const Matrix<T> &) const = 0;

    // Producto por una constante
    virtual Matrix<T> &operator*=(const T &) = 0;
    virtual const Matrix<T> operator*(const T &) const = 0;

    // Producto por una matriz
    virtual Matrix<T> &operator*=(const Matrix<T> &) = 0;
    virtual const Matrix<T> operator*(const Matrix<T> &) const = 0;
    virtual ~Matrix() { }
private:
    size_t N;
    size_t M;
};

/*
* Matriz ineficiente. Todos los algorítmos son malísimos.
*/
template <typename T>
class InefficientMatrix : public Matrix<T> {
public:
    // Constructor
    InefficientMatrix(std::size_t N, std::size_t M) : Matrix(N, M) {
        matrix = new T[N][M];

        for (auto i; i < N; ++i) {
            for (auto j; j < M; ++j) {
                (*this)(i, j) = 0;
            }
        }
    }

    // Constructor por copia
    InefficientMatrix<T>(const InefficientMatrix<T> &m) : Matrix(m.rows(), m.columns()) {
        matrix = new T[this->rows()][this->columns()];

        for (auto i; i < this->rows(); ++i) {
            for (auto j; j < this->columns(); ++j) {
                (*this)(i, j) = m(i, j);
            }
        }
    }

    // Lector de indice
    T &operator()(std::size_t i, std::size_t j) {
        return matrix[i][j];
    }

    // Lector de indice constante
    const T &operator()(std::size_t i, std::size_t j) const {
        return matrix[i][j];
    }

    // Printing
    std::ostream &operator<<(std::ostream &os, InefficientMatrix<T> const &m) {
        for (auto i; i < this->rows(); ++i) {
            for (auto j; j < this->columns(); ++j) {
                os << (*this)(i, j) << ' ';
            }

            os << '\n';
        }

        os << '\n';

        return os;
    }

    // Asignación
    InefficientMatrix<T> &operator=(const InefficientMatrix<T> &m) {
        if (this == &m) {
            return *this;
        } else if (this->rows() != m.rows() || this->columns() != m.columns()) {
            throw new std::out_of_range("Different row or column number");
        } else {
            for (auto i; i < this->rows(); ++i) {
                for (auto j; j < this->columns(); ++j) {
                    (*this)(i, j) = m(i, j);
                }
            }
        }

        return *this;
    }

    // Suma
    InefficientMatrix<T> &operator+=(const InefficientMatrix<T> &m) {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            throw new std::out_of_range("Different row or column number");
        } else {
            for (auto i; i < this->rows(); ++i) {
                for (auto j; j < this->columns(); ++j) {
                    (*this)(i, j) += m(i, j);
                }
            }

            return *this;
        }
    }

    const InefficientMatrix<T> operator+(const InefficientMatrix<T> &m) const {
        InefficientMatrix<T> output = *this;
        output += m;
        return output;
    };

    // Igualdad
    bool operator==(const InefficientMatrix<T> &m) const {
        if (this->rows() != m.rows() || this->columns() != m.columns()) {
            return false;
        } else {
            for (auto i; i < this->rows(); ++i) {
                for (auto j; j < this->columns(); ++j) {
                    if ((*this)(i, j) != m(i, j)) {
                        return false;
                    }
                }
            }

            return true;
        }
    };

    bool operator!=(const InefficientMatrix<T> &m) const {
        return !(*this == m);
    }

    // Producto por una constante
    InefficientMatrix<T> &operator*=(const T &c) {
        for (auto i; i < this->rows(); i++) {
            for (auto j; j < this->columns(); j++) {
                (*this)(i, j) *= c;
            }
        }

        return *this;
    };

    const InefficientMatrix<T> operator*(const T &c) {
        InefficientMatrix<T> output = *this;
        output *= c;
        return output;
    };

    // Producto por otra matriz
    InefficientMatrix<T> &operator*=(const InefficientMatrix<T> &m) {
        // TODO: Ver como hacer esto bien. Tiene que ser destructivo pero hacer el producto igual... Si esto no sale, cambiar el de abajo.
    };

    const InefficientMatrix<T> operator*(const InefficientMatrix<T> &m) {
        InefficientMatrix<T> output = *this;
        output *= m;
        return output;
    };

    virtual ~InefficientMatrix() {
        delete[] this->matrix;
    }
private:
    // Matrix
    T **matrix;
};

#endif //_TP1_MATRIX_H_
