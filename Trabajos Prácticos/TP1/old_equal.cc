BandMatrix &BandMatrix::operator*(const BandMatrix &m) {
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

bool BandMatrix::operator==(const BandMatrix &m) const {
    if (this->rows() != m.rows() || this->columns() != m.columns()) {
        return false;
    } else {
        // Queremos: empezar por los numeros que estan fuera de la interseccion de las bandas
        // Asumimos que las chances son que ahi alla algo distinto de 0, sino simplemente podriamos
        // cambiar la cantidad de bandas de la matriz a la cantidad que sea menor!
        long int diagonal = static_cast<long int>(std::min(this->rows(), this->columns()));

        long int max_lower = static_cast<long int>(std::max(this->lower_bandwidth(), m.lower_bandwidth()));
        long int min_lower = static_cast<long int>(std::min(this->lower_bandwidth(), m.lower_bandwidth()));

        // Chequeamos la parte de la izquierda de la diagonal - max_lower
        if (max_lower > min_lower) {
            for (long int d = 0L; d < diagonal; ++d) {
                // No queremos irnos fuera de rango
                long int bound = std::max(d - min_lower, 0L);

                for (long int j = std::max(d - max_lower, 0L); j < bound; ++j) {
                    if ((*this)(d, j) != m(d, j)) {
                        return false;
                    }
                }
            }
        }

        long int max_upper = static_cast<long int>(std::max(this->upper_bandwidth(), m.upper_bandwidth()));
        long int min_upper = static_cast<long int>(std::min(this->upper_bandwidth(), m.upper_bandwidth()));

        // Chequeamos la parte de la derecha de la diagonal + min_lower
        if (max_upper > min_upper) {
            for (long int d = 0L; d < diagonal; ++d) {
                // No queremos irnos fuera de rango
                long int bound = std::min(d + max_upper, static_cast<long int>(this->columns()));

                for (long int j = std::min(d + min_upper, static_cast<long int>(this->columns())); j < bound; ++j) {
                    if ((*this)(d, j) != m(d, j)) {
                        return false;
                    }
                }
            }
        }

        // Si llegamos hasta aca, es porque los numeros que estan fuera de la interseccion son todos iguales,
        // o porque tienen los mismos valores de upper y lower bandwidth
        // Chequeamos la parte de la intersección de la banda
        for (long int d = 0L; d < diagonal; ++d) {
            long int bound = std::min(d + static_cast<long int>(this->upper_bandwidth()), static_cast<long int>(this->columns()));

            for (long int j = std::max(d - static_cast<long int>(this->lower_bandwidth()), 0L); j < bound; ++j) {
                if ((*this)(d, j) != m(d, j)) {
                    return false;
                }
            }
        }

        return true;
    }
}