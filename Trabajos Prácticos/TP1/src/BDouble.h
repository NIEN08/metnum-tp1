//
// Created by Julian Bayardo on 3/25/15.
//

#ifndef _TP1_BDOUBLE_H_
#define _TP1_BDOUBLE_H_

#include <cmath>
#include <limits>
#include <utility>

class BDouble {
public:
    BDouble() : x(0.0) { }
    BDouble(double x) : x(x) { }
    BDouble(float x): x(x) { }
    BDouble(const BDouble &d) : x(d.x) { }

    double get() const { return this->x; }

    BDouble &operator=(const BDouble &rhs) {
        this->x = rhs.x;
        return *this;
    }

    BDouble &operator+=(const BDouble &rhs) {
        this->x += rhs.x;
        return *this;
    }

    BDouble &operator-=(const BDouble &rhs) {
        this->x -= rhs.x;
        return *this;
    }

    BDouble &operator*=(const BDouble &rhs) {
        this->x *= rhs.x;
        return *this;
    }

    BDouble &operator/=(const BDouble &rhs) {
        this->x /= rhs.x;
        return *this;
    }

    BDouble &operator=(const double &rhs) {
        this->x = rhs;
        return *this;
    }

    BDouble &operator+=(const double &rhs) {
        this->x += rhs;
        return *this;
    }

    BDouble &operator-=(const double &rhs) {
        this->x -= rhs;
        return *this;
    }

    BDouble &operator*=(const double &rhs) {
        this->x *= rhs;
        return *this;
    }

    BDouble &operator/=(const double &rhs) {
        this->x /= rhs;
        return *this;
    }

    BDouble &operator++() {
        this->x++;
        return *this;
    }

    BDouble &operator--() {
        this->x--;
        return *this;
    }

    bool operator==(const BDouble &rhs) const {
        return std::fabs(this->x - rhs.x) < std::numeric_limits<double>::epsilon();
    }

    bool operator<(const BDouble &rhs) const {
        return rhs.x - this->x < std::numeric_limits<double>::epsilon();
    }

    bool operator!=(const BDouble &m) const {
        return !(*this == m);
    }

    bool operator==(const double &rhs) const {
        return std::fabs(this->x - rhs) < std::numeric_limits<double>::epsilon();
    }

    bool operator<(const double &rhs) const {
        return rhs - this->x < std::numeric_limits<double>::epsilon();
    }

    bool operator!=(const double &m) const {
        return !(*this == m);
    }
private:
    double x;
};

BDouble operator+(const BDouble &lhs, const BDouble &rhs) {
    BDouble output(lhs);
    output += rhs;
    return output;
}

BDouble operator-(const BDouble &lhs, const BDouble &rhs) {
    BDouble output(lhs);
    output -= rhs;
    return output;
}

BDouble operator*(const BDouble &lhs, const BDouble &rhs) {
    BDouble output(lhs);
    output *= rhs;
    return output;
}

BDouble operator/(const BDouble &lhs, const BDouble &rhs) {
    BDouble output(lhs);
    output /= rhs;
    return output;
}

BDouble operator+(const BDouble &lhs, const double &rhs) {
    return lhs + BDouble(rhs);
}

BDouble operator-(const BDouble &lhs, const double &rhs) {
    return lhs - BDouble(rhs);
}

BDouble operator*(const BDouble &lhs, const double &rhs) {
    return lhs * BDouble(rhs);
}

BDouble operator/(const BDouble &lhs, const double &rhs) {
    return lhs / BDouble(rhs);
}

std::ostream &operator<<(std::ostream &os, const BDouble &m) {
    os << m.get();
    return os;
}

#endif //_TP1_BOXEDDOUBLE_H_
