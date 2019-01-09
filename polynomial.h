#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <QVector>

#include "range.h"

namespace mp {
class Polynomial
{
    Range _terms;

public:
    Polynomial();
    Polynomial(double d);
    Polynomial(double d, int degree);
    explicit Polynomial(const Range &coefs);

    static Polynomial bessel(int n);
    static Polynomial hermite(int n);
    static Polynomial laguerre(int n);
    static Polynomial legendre(int n);
    static Polynomial bernstein(int i, int n); //i <= n
    static Polynomial chebyshev(int n, const int &kind = 1);

    int degree() const;
    double coefficient(int term) const;

    static Polynomial add(const Polynomial &l, const Polynomial &r);
    static Polynomial sub(const Polynomial &l, const Polynomial &r);
    static Polynomial mul(const Polynomial &l, const Polynomial &r);

    QString toString(char f = 'g', int precision = 6) const;

    double operator ()(double x) const;
    Range operator ()(const Range &x) const;

    Polynomial derivate() const;
    Polynomial integrate(double constant = 0) const;

    Polynomial operator +(const Polynomial &p) const { return Polynomial::add(*this, p); }
    Polynomial operator -(const Polynomial &p) const { return Polynomial::sub(*this, p); }
    Polynomial operator *(const Polynomial &p) const { return Polynomial::mul(*this, p); }

};
}

#endif // POLYNOMIAL_H
