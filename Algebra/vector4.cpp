#include "vector4.h"

#include "../func.h"
#include "vector3.h"

#include "matrix4x4.h"

#include <math.h>

using namespace mp;

Vector4::Vector4(std::initializer_list<float> list)
{
    assert(list.size() <= 4);
    const float *d = list.begin();
    for(size_t i = 0; i < list.size(); ++i)
        this->_d[i] = d[i];
}

Vector4::Vector4(const Vector3 &cartesian, float w)
{
    this->_d[0] = cartesian.x();
    this->_d[1] = cartesian.y();
    this->_d[2] = cartesian.z();
    this->_d[3] = w;
}

Vector4::Vector4(float x, float y, float z, float w)
{
    this->_d[0] = x;
    this->_d[1] = y;
    this->_d[2] = z;
    this->_d[3] = w;
}

float *Vector4::data()
{
    return this->_d;
}

const float *Vector4::data() const
{
    return this->_d;
}

Vector4 Vector4::normalize() const
{
    float l = this->magnitude();
    return Vector4(this->_d[0] / l, this->_d[1] / l, this->_d[2] / l, this->_d[3] / l);
}

float Vector4::magnitude() const
{
    return mp::sqrt(this->sqrMagnitude());
}

float Vector4::sqrMagnitude() const
{
    return this->_d[0] * this->_d[0] + this->_d[1] * this->_d[1] + this->_d[2] * this->_d[2] + this->_d[3] * this->_d[3];
}

Vector4 Vector4::proj(const Vector4 &p) const
{
    float d = Vector4::dot(*this, p);
    return Vector4::mul(p.normalize(), d);
}

Vector4 Vector4::lerp(const Vector4 &to, float t) const
{
    if(t > 1) t = 1;
    else if(t < 0) t = 0;
    return *this + (*this - to) * t;
}

Vector4 Vector4::slerp(const Vector4 &to, float t) const
{
    float theta = ::acos(this->normalize() * to.normalize());
    return (::sin((1 - t) * theta) / ::sin(theta)) * *this + (::sin(t * theta) * ::sin(theta)) * to;
}

float Vector4::distance(const Vector4 &p) const
{
    return Vector4::sub(*this, p).magnitude();
}

Vector3 Vector4::toCartesian() const
{
    return Vector3(this->_d[0] / this->_d[3], this->_d[1] / this->_d[3], this->_d[2] / this->_d[3]);
}

std::string Vector4::toString() const
{
    std::string str = "(";
    str.append(std::to_string(this->_d[0])).append("; ");
    str.append(std::to_string(this->_d[1])).append("; ");
    str.append(std::to_string(this->_d[2])).append("; ");
    str.append(std::to_string(this->_d[3]));
    str.append(")");
    return str;
}

Vector4 Vector4::copy() const
{
    return Vector4(this->_d[0], this->_d[1], this->_d[2], this->_d[3]);
}

Vector4 Vector4::add(const Vector4 &left, const Vector4 &right)
{
    return Vector4(left._d[0] + right._d[0], left._d[1] + right._d[1], left._d[2] + right._d[2], left._d[3] + right._d[3]);
}

Vector4 Vector4::sub(const Vector4 &left, const Vector4 &right)
{
    return Vector4(left._d[0] - right._d[0], left._d[1] - right._d[1], left._d[2] - right._d[2], left._d[3] - right._d[3]);
}

float Vector4::dot(const Vector4 &left, const Vector4 &right)
{
    return left._d[0] * right._d[0] + left._d[1] * right._d[1] + left._d[2] * right._d[2] + left._d[3] * right._d[3];
}

Vector4 Vector4::mul(const Vector4 &left, const float &right)
{
    return Vector4(left._d[0] * right, left._d[1] * right, left._d[2] * right, left._d[3] * right);
}

Vector4 Vector4::mul(const Matrix4x4 &left, const Vector4 &right)
{
    Vector4 r;
    r[0] = left.at(0, 0) * right[0] +  left.at(1, 0) * right[1] +  left.at(2, 0) * right[2] + left.at(3, 0) * right[3];
    r[1] = left.at(0, 1) * right[0] +  left.at(1, 1) * right[1] +  left.at(2, 1) * right[2] + left.at(3, 1) * right[3];
    r[2] = left.at(0, 2) * right[0] +  left.at(1, 2) * right[1] +  left.at(2, 2) * right[2] + left.at(3, 2) * right[3];
    r[3] = left.at(0, 3) * right[0] +  left.at(1, 3) * right[1] +  left.at(2, 3) * right[2] + left.at(3, 2) * right[3];
    return r;
}

Vector4 Vector4::div(const Vector4 &left, const float &right)
{
    return Vector4(left._d[0] / right, left._d[1] / right, left._d[2] / right, left._d[3] / right);
}

Vector4 Vector4::neg(const Vector4 &v)
{
    return Vector4(-v._d[0], -v._d[1], -v._d[2], -v._d[3]);
}

bool Vector4::equal(const Vector4 &left, const Vector4 &right)
{
    return ::fabs(left[0] - right[0]) < 1e-15 && ::fabs(left[1] - right[1]) < 1e-15 && ::fabs(left[2] - right[2]) < 1e-15 && ::fabs(left[3] - right[3]) < 1e-15;
}
