#include "math.h"

#include "range.h"

mp::Range::Range(const QList<double> &args)
{
    this->resize(args.size());
    double *d = this->data();
    QList<double>::ConstIterator it = args.constBegin();
    for(int i = 0; i < args.size(); ++i){
        d[i] = *it;
        ++it;
    }
}

mp::Range::Range(const QVector<double> &args)
{
    this->resize(args.size());
    double *d = this->data();
    const double *b = args.data();
    for(int i = 0; i < args.size(); ++i){
        d[i] = b[i];
    }
}

mp::Range::Range(std::initializer_list<double> args) : QVector<double>(args.size())
{
    double *d = this->data();
    const double *b = args.begin();
    for(unsigned int i = 0; i < args.size(); ++i){
        d[i] = b[i];
    }
}

mp::Range::Range(const double &a, const double &b, const double &step) : QVector<double>(0)
{
    double min = ::fmin(a, b);
    double max = ::fmax(a, b);
    int count = (int)::qRound((max - min) / step);
    this->resize(count);
    double *y_data = this->data();
    for(int i = 0; i < count; ++i)
        y_data[i] = min + i * step;
}

mp::Range mp::Range::linspace(const double &a, const double &b, const int &count)
{
    double min = ::fmin(a, b);
    double max = ::fmax(a, b);
    double step = (max - min) / (count - 1);
    mp::Range y(count);
    double *y_data = y.data();
    for(int i = 0; i < count; ++i)
        y_data[i] = min + i * step;
    return y;
}

mp::Range mp::Range::add(const double &f) const
{
    Range copy(this->size());
    double *c_data = copy.data();
    const double *t_data = this->data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = t_data[i] + f;
    }
    return copy;
}

mp::Range mp::Range::sub(const double &f) const
{
    Range copy(this->size());
    double *c_data = copy.data();
    const double *t_data = this->data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = t_data[i] - f;
    }
    return copy;
}

mp::Range mp::Range::mul(const double &f) const
{
    Range copy(this->size());
    double *c_data = copy.data();
    const double *t_data = this->data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = t_data[i] * f;
    }
    return copy;
}

mp::Range mp::Range::div(const double &f) const
{
    Range copy(this->size());
    double *c_data = copy.data();
    const double *t_data = this->data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = t_data[i] / f;
    }
    return copy;
}

mp::Range mp::Range::add(const mp::Range &r) const
{
    Range copy(this->size());
    double *c_data = copy.data();
    const double *t_data = this->data();
    const double *r_data = r.data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = t_data[i] + r_data[i];
    }
    return copy;
}

mp::Range mp::Range::sub(const mp::Range &r) const
{
    Range copy(this->size());
    double *c_data = copy.data();
    const double *t_data = this->data();
    const double *r_data = r.data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = t_data[i] - r_data[i];
    }
    return copy;
}

mp::Range mp::Range::mul(const mp::Range &r) const
{
    Range copy(this->size());
    double *c_data = copy.data();
    const double *t_data = this->data();
    const double *r_data = r.data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = t_data[i] * r_data[i];
    }
    return copy;
}

mp::Range mp::Range::div(const mp::Range &r) const
{
    Range copy(this->size());
    double *c_data = copy.data();
    const double *t_data = this->data();
    const double *r_data = r.data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = t_data[i] / r_data[i];
    }
    return copy;
}

mp::Range mp::Range::add(const double &f, const mp::Range &r)
{
    Range copy(r.size());
    double *c_data = copy.data();
    const double *r_data = r.data();
    for(int i = 0; i < r.size(); ++i){
        c_data[i] = f + r_data[i];
    }
    return copy;
}

mp::Range mp::Range::sub(const double &f, const mp::Range &r)
{
    Range copy(r.size());
    double *c_data = copy.data();
    const double *r_data = r.data();
    for(int i = 0; i < r.size(); ++i){
        c_data[i] = f - r_data[i];
    }
    return copy;
}

mp::Range mp::Range::mul(const double &f, const mp::Range &r)
{
    Range copy(r.size());
    double *c_data = copy.data();
    const double *r_data = r.data();
    for(int i = 0; i < r.size(); ++i){
        c_data[i] = f * r_data[i];
    }
    return copy;
}

mp::Range mp::Range::div(const double &f, const mp::Range &r)
{
    Range copy(r.size());
    double *c_data = copy.data();
    const double *r_data = r.data();
    for(int i = 0; i < r.size(); ++i){
        c_data[i] = f / r_data[i];
    }
    return copy;
}

mp::Range mp::Range::neg() const
{
    Range copy(this->size());
    double *c_data = copy.data();
    const double *t_data = this->data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = -t_data[i];
    }
    return copy;
}

mp::Range mp::Range::pow(const double &f) const
{
    mp::Range y(this->size());
    double *y_data = y.data();
    const double *x_data = this->data();
    for(int i = 0; i < this->size(); ++i)
        y_data[i] = ::pow(x_data[i], f);
    return y;
}

mp::Range mp::Range::copy() const
{
    mp::Range c(this->size());
    const double *t_data = this->data();
    double *c_data = c.data();
    for(int i = 0; i < this->size(); ++i){
        c_data[i] = t_data[i];
    }
    return c;
}

mp::Range mp::Range::concat(const mp::Range &right) const
{
    mp::Range c(this->size() + right.size());
    double *d = c.data();
    const double *t = this->data();
    const double *r = right.data();
    for(int i = 0; i < this->size(); ++i)
        d[i] = t[i];
    for(int i = 0; i < right.size(); ++i)
        d[i + this->size()] = r[i];
    return c;
}
