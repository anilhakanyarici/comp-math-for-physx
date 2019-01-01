#include "complex.h"

#include "func.h"

mp::Complex::Complex(double r) { this->_r = r; this->_i = 0; }

mp::Complex::Complex(double r, double i) : _r(r), _i(i) { }

double mp::Complex::real() const { return this->_r; }

double &mp::Complex::real() { return this->_r; }

double mp::Complex::imaginer() const { return this->_i; }

double &mp::Complex::imaginer() { return this->_i; }

double mp::Complex::magnitude() const { return mp::sqrt(mp::square(this->_r) + mp::square(this->_i)); }

mp::Complex mp::Complex::conjugate() const { return mp::Complex(this->_r, -this->_i); }

mp::Complex mp::Complex::operator +(const double &r) { return Complex(this->_r + r, this->_i); }

mp::Complex mp::Complex::operator -(const double &r) { return Complex(this->_r - r, this->_i); }

mp::Complex mp::Complex::operator *(const double &r) { return Complex(this->_r * r, this->_i * r); }

mp::Complex mp::Complex::operator /(const double &r) { return Complex(this->_r / r, this->_i / r); }

mp::Complex mp::Complex::operator +(const mp::Complex &c) { return Complex(this->_r + c._r, this->_i + c._i); }

mp::Complex mp::Complex::operator -(const mp::Complex &c) { return Complex(this->_r - c._r, this->_i - c._i); }

mp::Complex mp::Complex::operator *(const mp::Complex &c)
{
    double r = this->_r * c._r - this->_i * c._i;
    double i = this->_r * c._i + this->_i * c._r;
    return Complex(r, i);
}

mp::Complex mp::Complex::operator /(const mp::Complex &c)
{
    double r = this->_r * c._r + this->_i * c._i;
    double i = this->_i * c._r - this->_r * c._i;
    double m = c.magnitude();
    return Complex(r / m, i / m);
}

mp::Complex mp::Complex::add(const double &r, const mp::Complex &c) { return Complex(c._r + r, c._i); }

mp::Complex mp::Complex::sub(const double &r, const mp::Complex &c) { return Complex(r - c._r, c._i); }

mp::Complex mp::Complex::mul(const double &r, const mp::Complex &c) { return Complex(r * c._r, r * c._i); }

mp::Complex mp::Complex::div(const double &r, const mp::Complex &c)
{
    double cr = r * c._r;
    double ci = -r * c._i;
    double m = c.magnitude();
    return Complex(cr / m, ci / m);
}

mp::ComplexRange::ComplexRange(const mp::Range &r)
{
    this->_r = mp::Range(r.size());
    double *cr_data = this->_r.data();
    const double *r_data = r.data();
    for(int k = 0; k < r.size(); ++k)
        cr_data[k] = r_data[k];

    this->_i = mp::Range(r.size());
}

mp::ComplexRange::ComplexRange(int size)
{
    this->_r = Range(size);
    this->_i = Range(size);
}

mp::ComplexRange::ComplexRange(const mp::Range &r, const mp::Range &i)
{
    this->_r = mp::Range(r.size());
    double *cr_data = this->_r.data();
    const double *r_data = r.data();
    this->_i = mp::Range(i.size());
    double *ci_data = this->_i.data();
    const double *i_data = i.data();
    for(int k = 0; k < r.size(); ++k){
        cr_data[k] = r_data[k];
        ci_data[k] = i_data[k];
    }
}

void mp::ComplexRange::append(const mp::Complex &c)
{
    this->_r.append(c._r);
    this->_i.append(c._i);
}

void mp::ComplexRange::append(const double &r, const double &i)
{
    this->_r.append(r);
    this->_i.append(i);
}

mp::Range mp::ComplexRange::real() const { return this->_r; }

mp::Range &mp::ComplexRange::real() { return this->_r; }

mp::Range mp::ComplexRange::imaginer() const { return this->_i; }

mp::Range &mp::ComplexRange::imaginer() { return this->_i; }

mp::Range mp::ComplexRange::magnitude() const { return mp::sqrt(mp::square(this->_r) + mp::square(this->_i)); }

mp::ComplexRange mp::ComplexRange::conjugate() const
{
    return mp::ComplexRange(this->_r, -this->_i);
}

mp::ComplexRange mp::ComplexRange::operator +(const double &r)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] + r;
        ci_data[k] = i_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator -(const double &r)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] - r;
        ci_data[k] = i_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator *(const double &r)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] * r;
        ci_data[k] = i_data[k] * r;
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator /(const double &r)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] / r;
        ci_data[k] = i_data[k] / r;
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator +(const mp::Complex &c)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] + c._r;
        ci_data[k] = i_data[k] + c._i;
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator -(const mp::Complex &c)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] - c._r;
        ci_data[k] = i_data[k] - c._i;
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator *(const mp::Complex &c)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] * c._r - i_data[k] * c._i;
        ci_data[k] = r_data[k] * c._i + i_data[k] * c._r;
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator /(const mp::Complex &c)
{
    double m = c.magnitude();
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = (r_data[k] * c._r + i_data[k] * c._i) / m;
        ci_data[k] = (i_data[k] * c._r - r_data[k] * c._i) / m;
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator +(const mp::ComplexRange &c)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    const double *ccr_data = c._r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    const double *cci_data = c._i.data();

    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] + ccr_data[k];
        ci_data[k] = i_data[k] + cci_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator -(const mp::ComplexRange &c)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    const double *ccr_data = c._r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    const double *cci_data = c._i.data();

    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] - ccr_data[k];
        ci_data[k] = i_data[k] - cci_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator *(const mp::ComplexRange &c)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    const double *ccr_data = c._r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    const double *cci_data = c._i.data();

    for(int k = 0; k < this->_r.size(); ++k){
        cr_data[k] = r_data[k] * ccr_data[k] - i_data[k] * cci_data[k];
        ci_data[k] = r_data[k] * cci_data[k] + i_data[k] * ccr_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::operator /(const mp::ComplexRange &c)
{
    ComplexRange range(this->size());
    double *cr_data = range._r.data();
    const double *r_data = this->_r.data();
    const double *ccr_data = c._r.data();
    double *ci_data = range._i.data();
    const double *i_data = this->_i.data();
    const double *cci_data = c._i.data();

    for(int k = 0; k < this->_r.size(); ++k){
        double m = mp::sqrt(ccr_data[k] * ccr_data[k] + cci_data[k] * cci_data[k]);
        cr_data[k] = (r_data[k] * ccr_data[k] + i_data[k] * cci_data[k]) / m;
        ci_data[k] = (i_data[k] * ccr_data[k] - r_data[k] * cci_data[k]) / m;
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::add(const double &r, const ComplexRange &c)
{
    ComplexRange range(c.size());
    double *cr_data = range._r.data();
    const double *r_data = c._r.data();
    double *ci_data = range._i.data();
    const double *i_data = c._i.data();
    for(int k = 0; k < c._r.size(); ++k){
        cr_data[k] = r + r_data[k];
        ci_data[k] = i_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::sub(const double &r, const mp::ComplexRange &c)
{
    ComplexRange range(c.size());
    double *cr_data = range._r.data();
    const double *r_data = c._r.data();
    double *ci_data = range._i.data();
    const double *i_data = c._i.data();
    for(int k = 0; k < c._r.size(); ++k){
        cr_data[k] = r - r_data[k];
        ci_data[k] = i_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::mul(const double &r, const mp::ComplexRange &c)
{
    ComplexRange range(c.size());
    double *cr_data = range._r.data();
    const double *r_data = c._r.data();
    double *ci_data = range._i.data();
    const double *i_data = c._i.data();
    for(int k = 0; k < c._r.size(); ++k){
        cr_data[k] = r * r_data[k];
        ci_data[k] = r * i_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::div(const double &r, const mp::ComplexRange &c)
{
    ComplexRange range(c.size());
    double *cr_data = range._r.data();
    const double *r_data = c._r.data();
    double *ci_data = range._i.data();
    const double *i_data = c._i.data();
    for(int k = 0; k < c._r.size(); ++k){
        double m = mp::sqrt(r_data[k] * r_data[k] + i_data[k] * i_data[k]);
        cr_data[k] = (r * r_data[k]) / m;
        ci_data[k] = (-r * i_data[k]) / m;
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::add(const mp::Complex &c, const mp::ComplexRange &cr)
{
    ComplexRange range(cr.size());
    double *range_r_data = range._r.data();
    const double *crr_data = cr._r.data();
    double *range_i_data = range._i.data();
    const double *cri_data = cr._i.data();
    for(int k = 0; k < cr._r.size(); ++k){
        range_r_data[k] = c._r + crr_data[k];
        range_i_data[k] = c._r + cri_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::sub(const mp::Complex &c, const mp::ComplexRange &cr)
{
    ComplexRange range(cr.size());
    double *range_r_data = range._r.data();
    const double *crr_data = cr._r.data();
    double *range_i_data = range._i.data();
    const double *cri_data = cr._i.data();
    for(int k = 0; k < cr._r.size(); ++k){
        range_r_data[k] = c._r - crr_data[k];
        range_i_data[k] = c._r - cri_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::mul(const mp::Complex &c, const mp::ComplexRange &cr)
{
    ComplexRange range(cr.size());
    double *range_r_data = range._r.data();
    const double *crr_data = cr._r.data();
    double *range_i_data = range._i.data();
    const double *cri_data = cr._i.data();
    for(int k = 0; k < cr._r.size(); ++k){
        range_r_data[k] = c._r + crr_data[k] - c._i * cri_data[k];
        range_i_data[k] = c._r + cri_data[k] + c._i * crr_data[k];
    }
    return range;
}

mp::ComplexRange mp::ComplexRange::div(const mp::Complex &c, const mp::ComplexRange &cr)
{
    ComplexRange range(cr.size());
    double *range_r_data = range._r.data();
    const double *crr_data = cr._r.data();
    double *range_i_data = range._i.data();
    const double *cri_data = cr._i.data();
    for(int k = 0; k < cr._r.size(); ++k){
        double m = mp::sqrt(crr_data[k] * crr_data[k] + cri_data[k] * cri_data[k]);
        range_r_data[k] = (c._r + crr_data[k] + c._i * cri_data[k]) / m;
        range_i_data[k] = (c._r + cri_data[k] - c._i * crr_data[k]) / m;
    }
    return range;
}
