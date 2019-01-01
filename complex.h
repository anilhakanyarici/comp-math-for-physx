#ifndef COMPLEX_H
#define COMPLEX_H

#include "range.h"

namespace mp {
class Complex
{
    friend class ComplexRange;
    double _r, _i;
public:
    Complex(double r);
    explicit Complex(double r = 0, double i = 0);

    double real() const;
    double &real();
    double imaginer() const;
    double &imaginer();

    double magnitude() const;
    Complex conjugate() const;

    Complex operator +(const double &r);
    Complex operator -(const double &r);
    Complex operator *(const double &r);
    Complex operator /(const double &r);

    Complex operator +(const Complex &c);
    Complex operator -(const Complex &c);
    Complex operator *(const Complex &c);
    Complex operator /(const Complex &c);

    static Complex add(const double &r, const Complex &c);
    static Complex sub(const double &r, const Complex &c);
    static Complex mul(const double &r, const Complex &c);
    static Complex div(const double &r, const Complex &c);

    friend Complex operator +(const double &r, const Complex &c) { return Complex::add(r, c); }
    friend Complex operator -(const double &r, const Complex &c) { return Complex::sub(r, c); }
    friend Complex operator *(const double &r, const Complex &c) { return Complex::mul(r, c); }
    friend Complex operator /(const double &r, const Complex &c) { return Complex::div(r, c); }
};

class ComplexRange
{
    Range _r, _i;
public:
    ComplexRange(const Range &r);
    explicit ComplexRange(int size = 0);
    ComplexRange(const Range &r, const Range &i);

    int size() const { return this->_r.size(); }
    void append(const Complex &c);
    void append(const double &r, const double &i);

    Range real() const;
    Range &real();
    Range imaginer() const;
    Range &imaginer();

    Range magnitude() const;
    ComplexRange conjugate() const;

    Complex operator[](int i) { return Complex(this->_r[i], this->_i[i]); }

    ComplexRange operator +(const double &r);
    ComplexRange operator -(const double &r);
    ComplexRange operator *(const double &r);
    ComplexRange operator /(const double &r);

    ComplexRange operator +(const Complex &c);
    ComplexRange operator -(const Complex &c);
    ComplexRange operator *(const Complex &c);
    ComplexRange operator /(const Complex &c);

    ComplexRange operator +(const ComplexRange &c);
    ComplexRange operator -(const ComplexRange &c);
    ComplexRange operator *(const ComplexRange &c);
    ComplexRange operator /(const ComplexRange &c);

    static ComplexRange add(const double &r, const ComplexRange &c);
    static ComplexRange sub(const double &r, const ComplexRange &c);
    static ComplexRange mul(const double &r, const ComplexRange &c);
    static ComplexRange div(const double &r, const ComplexRange &c);

    static ComplexRange add(const Complex &c, const ComplexRange &cr);
    static ComplexRange sub(const Complex &c, const ComplexRange &cr);
    static ComplexRange mul(const Complex &c, const ComplexRange &cr);
    static ComplexRange div(const Complex &c, const ComplexRange &cr);

    friend ComplexRange operator +(const double &r, const ComplexRange &c) { return ComplexRange::add(r, c); }
    friend ComplexRange operator -(const double &r, const ComplexRange &c) { return ComplexRange::sub(r, c); }
    friend ComplexRange operator *(const double &r, const ComplexRange &c) { return ComplexRange::mul(r, c); }
    friend ComplexRange operator /(const double &r, const ComplexRange &c) { return ComplexRange::div(r, c); }

    friend ComplexRange operator +(const Complex &r, const ComplexRange &c) { return ComplexRange::add(r, c); }
    friend ComplexRange operator -(const Complex &r, const ComplexRange &c) { return ComplexRange::sub(r, c); }
    friend ComplexRange operator *(const Complex &r, const ComplexRange &c) { return ComplexRange::mul(r, c); }
    friend ComplexRange operator /(const Complex &r, const ComplexRange &c) { return ComplexRange::div(r, c); }
};
}

#endif // COMPLEX_H
