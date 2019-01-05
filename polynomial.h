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
};
}

#endif // POLYNOMIAL_H
