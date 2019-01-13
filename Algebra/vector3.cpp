#include "vector3.h"
#include "matrix3x3.h"
#include "quaternion.h"

#include "../func.h"

#include <math.h>

using namespace mp;

Vector3::Vector3(std::initializer_list<float> list)
{
    assert(list.size() <= 3);
    const float *d = list.begin();
    for(size_t i = 0; i < list.size(); ++i)
        this->_d[i] = d[i];
}

Vector3::Vector3(float x, float y, float z)
{
    this->_d[0] = x; this->_d[1] = y; this->_d[2] = z;
}

float *Vector3::data()
{
    return this->_d;
}

const float *Vector3::data() const
{
    return this->_d;
}

const float &Vector3::operator [](const int &i) const
{
    return this->_d[i];
}

float &Vector3::operator [](const int &i)
{
    return this->_d[i];
}

const float &Vector3::at(const int &i) const
{
    return this->_d[i];
}

float &Vector3::at(const int &i)
{
    return this->_d[i];
}

Vector3 Vector3::normalize() const
{
    float l = this->magnitude();
    return Vector3(this->_d[0] / l, this->_d[1] / l, this->_d[2] / l);
}

float Vector3::magnitude() const
{
    return mp::sqrt(this->sqrMagnitude());
}

float Vector3::sqrMagnitude() const
{
    return this->_d[0] * this->_d[0] + this->_d[1] * this->_d[1] + this->_d[2] * this->_d[2];
}

Vector3 Vector3::proj(const Vector3 &p) const
{
    float d = Vector3::dot(*this, p);
    return Vector3::mul(p.normalize(), d);
}

Vector3 Vector3::lerp(const Vector3 &to, float t) const
{
    if(t > 1) t = 1;
    else if(t < 0) t = 0;
    return *this + (*this - to) * t;
}

Vector3 Vector3::slerp(const Vector3 &to, float t) const
{
    float theta = ::acos(this->normalize() * to.normalize());
    return (::sin((1 - t) * theta) / ::sin(theta)) * *this + (::sin(t * theta) * ::sin(theta)) * to;
}

float Vector3::distance(const Vector3 &p) const
{
    return Vector3::sub(*this, p).magnitude();
}

Vector3 Vector3::rotate(const Quaternion &q) const
{
    Quaternion p = Quaternion::mul(q, *this);
    p = Quaternion::mul(p, q.reciprocal());
    return p.vector();
}

Vector3 Vector3::rotate(float angle, const Vector3 &rotation_axis) const
{
    Quaternion q(angle / 2, rotation_axis);
    return Vector3::rotate(q);
}

std::string Vector3::toString() const
{
    std::string str = "(";
    str.append(std::to_string(this->_d[0])).append("; ");
    str.append(std::to_string(this->_d[1])).append("; ");
    str.append(std::to_string(this->_d[2]));
    str.append(")");
    return str;
}

Vector3 Vector3::copy() const
{
    return Vector3(this->_d[0], this->_d[1], this->_d[2]);
}

Vector3 Vector3::add(const Vector3 &left, const Vector3 &right)
{
    return Vector3(left._d[0] + right._d[0], left._d[1] + right._d[1], left._d[2] + right._d[2]);
}

Vector3 Vector3::sub(const Vector3 &left, const Vector3 &right)
{
    return Vector3(left._d[0] - right._d[0], left._d[1] - right._d[1], left._d[2] - right._d[2]);
}

float Vector3::dot(const Vector3 &left, const Vector3 &right)
{
    return left._d[0] * right._d[0] + left._d[1] * right._d[1] + left._d[2] * right._d[2];
}

Vector3 Vector3::cross(const Vector3 &left, const Vector3 &right)
{
    float x = left._d[1] * right._d[2] - left._d[2] * right._d[1];
    float y = left._d[0] * right._d[2] - left._d[2] * right._d[0];
    float z = left._d[0] * right._d[1] - left._d[1] * right._d[0];
    return Vector3(x, -y, z);
}

Vector3 Vector3::mul(const Vector3 &left, const float &right)
{
    return Vector3(left._d[0] * right, left._d[1] * right, left._d[2] * right);
}

Vector3 Vector3::mul(const Matrix3x3 &m, const Vector3 &v)
{
    float x = m[0][0] * v.x() + m[0][1] * v.y() + m[0][2] * v.z();
    float y = m[1][0] * v.x() + m[1][1] * v.y() + m[1][2] * v.z();
    float z = m[2][0] * v.x() + m[2][1] * v.y() + m[2][2] * v.z();
    return Vector3(x, y, z);
}

Vector3 Vector3::div(const Vector3 &left, const float &right)
{
    return Vector3(left._d[0] / right, left._d[1] / right, left._d[2] / right);
}

bool Vector3::equal(const Vector3 &left, const Vector3 &right)
{
    return ::fabs(left[0] - right[0]) < 1e-15 && ::fabs(left[1] - right[1]) < 1e-15 && ::fabs(left[2] - right[2]) < 1e-15;
}
