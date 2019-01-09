#include "polynomial.h"

using namespace mp;

Polynomial::Polynomial()
{
    this->_terms.append(0);
}

Polynomial::Polynomial(double d)
{
    this->_terms.append(d);
}

Polynomial::Polynomial(double d, int degree)
{
    this->_terms = Range(degree + 1);
    this->_terms[degree] = d;
}

Polynomial::Polynomial(const Range &coefs)
{
    if(coefs.size() == 0)
        this->_terms.append(0);
    else {
        int size = coefs.size();
        const double *data = coefs.data();
        int degree = size;
        for(int i = size - 1; i >= 0; --i){
            if(data[i] == 0)
                --degree;
            else break;
        }
        for(int i = 0; i < degree; ++i){
            this->_terms.append(data[i]);
        }
        if(degree == 0)
            this->_terms.append(0);
    }
}

Polynomial Polynomial::hermite(int n)
{
    Polynomial hn((n & 1) == 0 ? 1 : -1);
    Polynomial min2x(mp::Range({0, -2}));
    for(int i = 0; i < n; ++i){
        hn = Polynomial::add(hn.derivate(), Polynomial::mul(min2x, hn));
    }
    return hn;
}

int Polynomial::degree() const { return this->_terms.size() - 1; }

double Polynomial::coefficient(int term) const
{
    if (term < 0 && this->_terms.size() <= term)
        return 0;
    return this->_terms[0];
}

Polynomial Polynomial::add(const Polynomial &l, const Polynomial &r)
{
    const double *l_data = l._terms.data();
    const double *r_data = r._terms.data();
    int r_size = r._terms.size();
    int l_size = l._terms.size();

    if(l_size > r_size){
        Range coefs(l_size);
        double *c_data = coefs.data();
        for(int i = 0; i < r_size; ++i)
            c_data[i] = l_data[i] + r_data[i];
        for(int i = r_size; i < l_size; ++i)
            c_data[i] = l_data[i];
        return Polynomial(coefs);

    } else {
        Range coefs(r_size);
        double *c_data = coefs.data();
        for(int i = 0; i < l_size; ++i)
            c_data[i] = l_data[i] + r_data[i];
        for(int i = l_size; i < r_size; ++i)
            c_data[i] = r_data[i];
        return Polynomial(coefs);
    }
}

Polynomial Polynomial::sub(const Polynomial &l, const Polynomial &r)
{
    const double *l_data = l._terms.data();
    const double *r_data = r._terms.data();
    int r_size = r._terms.size();
    int l_size = l._terms.size();

    if(l_size > r_size){
        Range coefs(l_size);
        double *c_data = coefs.data();
        for(int i = 0; i < r_size; ++i)
            c_data[i] = l_data[i] - r_data[i];
        for(int i = r_size; i < l_size; ++i)
            c_data[i] = l_data[i];
        return Polynomial(coefs);

    } else {
        Range coefs(r_size);
        double *c_data = coefs.data();
        for(int i = 0; i < l_size; ++i)
            c_data[i] = l_data[i] - r_data[i];
        for(int i = l_size; i < r_size; ++i)
            c_data[i] = -r_data[i];
        return Polynomial(coefs);
    }
}

Polynomial Polynomial::mul(const Polynomial &l, const Polynomial &r)
{
    Range coefs(l.degree() + r.degree() + 1);

          double *p_data = coefs.data();
    const double *l_data = l._terms.data();
    const double *r_data = r._terms.data();

    int l_size = l._terms.size();
    int r_size = r._terms.size();

    for(int i = 0; i < l_size; ++i){
        for(int j = 0; j < r_size; ++j){
            p_data[i + j] += l_data[i] * r_data[j];
        }
    }
    return Polynomial(coefs);
}

QString Polynomial::toString(char f, int precision) const
{
    if(this->degree() == 0)
        return QString::number(this->coefficient(0));
    QString str;
    int i = 0;
    const double *data = this->_terms.data();
    int size = this->_terms.size();
    for( ; i < size; ++i){
        if(data[i] != 0 && i > 0) {
            str = QString("%1*x^%2").arg(QString::number(data[i], f, precision), QString::number(i));
            break;
        } else if(data[i] != 0 && i == 0) {
            str = QString("%1").arg(QString::number(data[i], f, precision));
            break;
        }
    }
    ++i;
    for( ; i < size; ++i){
        if(data[i] != 0)
            str.append(QString(" + %1*x^%2").arg(QString::number(data[i], f, precision), QString::number(i)));
    }
    return str;
}

double Polynomial::operator ()(double x) const
{
    double pow = 1;
    double y = this->coefficient(0);
    const double *p_data = this->_terms.data();
    int size = this->_terms.size();
    for(int i = 1; i < size; ++i){
        pow *= x;
        y += p_data[i] * pow;
    }
    return y;
}

Range Polynomial::operator ()(const Range &x) const
{
    int size = x.size();
    mp::Range y(size);
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < size; ++i)
        y_data[i] = this->operator ()(x_data[i]);
    return y;
}

Polynomial Polynomial::derivate() const
{
    Range coefs(this->_terms.size() - 1);
    const double *t_data = this->_terms.data();
    double *c_data = coefs.data();
    int size = coefs.size();
    int d = 0;
    for(int i = 0; i < size; ){
        d = i + 1;
        c_data[i] = t_data[d] * d;
        i = d;
    }
    return Polynomial(coefs);
}

Polynomial Polynomial::integrate(double constant) const
{
    Range coefs(this->_terms.size() + 1);
    const double *t_data = this->_terms.data();
    double *c_data = coefs.data();
    int size = coefs.size();
    for(int i = 1; i < size; ++i){
        c_data[i] = t_data[i - 1] / i;
    }
    c_data[0] = constant;
    return Polynomial(coefs);
}
