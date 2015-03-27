// Igualdad
// Esta version anda bien pero no se si es eficiente.
// Recorre todas las posiciones de la matriz y se fija si todos los elementos son iguales
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