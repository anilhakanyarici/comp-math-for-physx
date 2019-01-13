#include "vector.h"

#include <math.h>

#include "vector3.h"
#include "vector4.h"

using namespace mp;

struct Vector::pimpl
{
    float *data = nullptr;
    int size = 0;
    bool copy = false;

    ~pimpl()
    {
        if(this->copy)
            ::free(this->data);
    }
};

Vector::Vector()
{
    this->_pimpl = std::shared_ptr<Vector::pimpl>(new Vector::pimpl());
}

Vector::Vector(const Vector3 &v)
{
    this->_pimpl = std::shared_ptr<Vector::pimpl>(new Vector::pimpl());
    this->_pimpl->data = (float*)::malloc(sizeof(float) * 3);
    this->_pimpl->data[0] = v[0];
    this->_pimpl->data[1] = v[1];
    this->_pimpl->data[2] = v[2];
    this->_pimpl->copy = true;
    this->_pimpl->size = 3;
}

Vector::Vector(const Vector4 &v)
{
    this->_pimpl = std::shared_ptr<Vector::pimpl>(new Vector::pimpl());
    this->_pimpl->data = (float*)::malloc(sizeof(float) * 4);
    this->_pimpl->data[0] = v[0];
    this->_pimpl->data[1] = v[1];
    this->_pimpl->data[2] = v[2];
    this->_pimpl->data[3] = v[3];
    this->_pimpl->copy = true;
    this->_pimpl->size = 4;
}

Vector::Vector(std::initializer_list<float> list)
{
    this->_pimpl = std::shared_ptr<Vector::pimpl>(new Vector::pimpl());
    this->_pimpl->data = (float*)::malloc(sizeof(float) * list.size());
    const float *v = list.begin();
    for(size_t i = 0; i < list.size(); ++i)
        this->_pimpl->data[i] = v[i];
    this->_pimpl->copy = true;
    this->_pimpl->size = list.size();
}

Vector::Vector(int size, bool zeros)
{
    this->_pimpl = std::shared_ptr<Vector::pimpl>(new Vector::pimpl());
    this->_pimpl->data = (float*)::malloc(sizeof(float) * size);
    if(zeros) {
        for(int i = 0; i < size; ++i)
            this->_pimpl->data[i] = 0;
    }
    this->_pimpl->copy = true;
    this->_pimpl->size = size;
}

Vector::Vector(const float *v, int size, bool copy)
{
    this->_pimpl = std::shared_ptr<Vector::pimpl>(new Vector::pimpl());
    this->_pimpl->copy = copy;
    this->_pimpl->size = size;

    if(copy){
        this->_pimpl->data = (float*)::malloc(sizeof(float) * size);
        for(int i = 0; i < size; ++i)
            this->_pimpl->data[i] = v[i];
    } else {
        this->_pimpl->data = const_cast<float*>(v);
    }
}

float *Vector::data()
{
    return this->_pimpl->data;
}

const float *Vector::data() const
{
    return this->_pimpl->data;
}

const float &Vector::operator [](const int &i) const
{
    return this->_pimpl->data[i];
}

float &Vector::operator [](const int &i)
{
    return this->_pimpl->data[i];
}

const float &Vector::at(const int &i) const
{
    return this->_pimpl->data[i];
}

float &Vector::at(const int &i)
{
    return this->_pimpl->data[i];
}

int Vector::size() const
{
    return this->_pimpl->size;
}

float Vector::magnitude() const
{
    return ::sqrt(this->sqrMagnitude());
}

float Vector::sqrMagnitude() const
{
    float m = 0;
    for(int i = 0; i < this->_pimpl->size; ++i){
        m += this->_pimpl->data[i] * this->_pimpl->data[i];
    }
    return m;
}

Vector Vector::normalize() const
{
    return *this / this->magnitude();
}

Vector Vector::proj(const Vector &p)
{
    float l = *this * p;
    float m = p.magnitude();
    return (l / m) * p;
}

Vector Vector::lerp(const Vector &to, float t) const
{
    if(t > 1) t = 1;
    else if(t < 0) t = 0;
    return *this + (*this - to) * t;
}

Vector Vector::slerp(const Vector &to, float t) const
{
    float theta = ::acos(this->normalize() * to.normalize());
    return (::sin((1 - t) * theta) / ::sin(theta)) * *this + (::sin(t * theta) * ::sin(theta)) * to;
}

float Vector::distance(const Vector &p) const
{
    return Vector::sub(*this, p).magnitude();
}

std::string Vector::toString() const
{
    std::string str = "(";
    int iLast = this->_pimpl->size - 1;
    for(int i = 0; i < iLast; ++i){
        str.append(std::to_string(this->_pimpl->data[i])).append("; ");
    }
    if(iLast >= 0)
        str.append(std::to_string(this->_pimpl->data[iLast]));
    str.append(")");
    return str;
}

Vector Vector::copy() const
{
    Vector v(this->_pimpl->size);
    for(int i = 0; i < this->_pimpl->size; ++i)
        v.at(i) = this->at(i);
    return v;
}

Vector Vector::add(const Vector &l, const Vector &r)
{
    int min_size = l.size() < r.size() ? l.size() : r.size();
    int max_size = l.size() > r.size() ? l.size() : r.size();
    Vector max = l.size() < r.size() ? r : l;
    Vector min = l.size() > r.size() ? r : l;

    Vector v(max_size);
    float *d = v.data();
    for(int i = 0; i < min_size; ++i)
        d[i] = max.at(i) + min.at(i);
    for(int i = min_size; i < max_size; ++i)
        d[i] = max.at(i);
    return v;
}

Vector Vector::sub(const Vector &l, const Vector &r)
{
    int min_size = l.size() < r.size() ? l.size() : r.size();
    int max_size = l.size() > r.size() ? l.size() : r.size();
    Vector max = l.size() < r.size() ? r : l;
    Vector min = l.size() > r.size() ? r : l;
    bool min_is_r = l.size() > r.size();

    Vector v(max_size);
    float *d = v.data();
    if(min_is_r){
        for(int i = 0; i < min_size; ++i)
            d[i] = max.at(i) - min.at(i);
        for(int i = min_size; i < max_size; ++i)
            d[i] = max.at(i);
    } else {
        for(int i = 0; i < min_size; ++i)
            d[i] = -max.at(i) + min.at(i);
        for(int i = min_size; i < max_size; ++i)
            d[i] = -max.at(i);
    }
    return v;
}

float Vector::mul(const Vector &l, const Vector &r)
{
    int min_size = l.size() < r.size() ? l.size() : r.size();
    float d = 0;
    for(int i = 0; i < min_size; ++i)
        d += l.at(i) * r.at(i);
    return d;
}

Vector Vector::mul(const Vector &l, const float &r)
{
    Vector c(l.data(), l.size());
    for(int i = 0; i < l.size(); ++i)
        c.at(i) = l.at(i) * r;
    return c;
}

Vector Vector::div(const Vector &l, const float &r)
{
    Vector c(l.data(), l.size());
    for(int i = 0; i < l.size(); ++i)
        c.at(i) = l.at(i) / r;
    return c;
}

bool Vector::equal(const Vector &l, const Vector &r)
{
    static float epsilon = 1e-15;

    Vector maxv = l;
    Vector minv = r;
    int maxl = l.size();
    int minl = r.size();
    if(l.size() < r.size())
    {
        maxv = r;
        minv = l;
        maxl = r.size();
        minl = l.size();
    }

    for(int i = 0; i < minl; ++i){
        if(::fabs(maxv[i] - minv[i]) >= epsilon)
            return false;
    }
    for(int i = minl; i < maxl; ++i){
        if(::fabs(maxv[i]) >= epsilon)
            return false;
    }
    return true;
}
