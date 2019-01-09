#ifndef RANGE_H
#define RANGE_H

#include <QVector>

namespace mp {
class Range : public QVector<double>
{
public:
    Range(int size = 0) : QVector<double>(size){}
    Range(const QList<double> &args);
    Range(const QVector<double> &args);
    Range(std::initializer_list<double> args);
    Range(const double &a, const double &b, const double &step);

    static Range linspace(const double &a, const double &b, const int &count);

    Range add(const double &f) const;
    Range sub(const double &f) const;
    Range mul(const double &f) const;
    Range div(const double &f) const;
    Range add(const Range &r) const;
    Range sub(const Range &r) const;
    Range mul(const Range &r) const;
    Range div(const Range &r) const;

    static Range add(const double &f, const Range &r);
    static Range sub(const double &f, const Range &r);
    static Range mul(const double &f, const Range &r);
    static Range div(const double &f, const Range &r);

    Range neg() const;
    Range pow(const double &f) const;

    Range copy() const;

    inline Range operator -() const { return this->neg(); }

    inline Range operator +(const Range &r) const { return this->add(r); }
    inline Range operator -(const Range &r) const { return this->sub(r); }
    inline Range operator *(const Range &r) const { return this->mul(r); }
    inline Range operator /(const Range &r) const { return this->div(r); }
    inline Range operator +(const double &f) const { return this->add(f); }
    inline Range operator -(const double &f) const { return this->sub(f); }
    inline Range operator *(const double &f) const { return this->mul(f); }
    inline Range operator /(const double &f) const { return this->div(f); }
    inline friend Range operator +(const double &f, const Range &r) { return Range::add(f, r); }
    inline friend Range operator -(const double &f, const Range &r) { return Range::sub(f, r); }
    inline friend Range operator *(const double &f, const Range &r) { return Range::mul(f, r); }
    inline friend Range operator /(const double &f, const Range &r) { return Range::div(f, r); }
};
}

#endif // RANGE_H
