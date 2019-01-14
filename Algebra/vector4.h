#ifndef VECTOR4_H
#define VECTOR4_H

#include <string>
#include <initializer_list>

#include "quaternion.h"

namespace mp
{

class Vector3;
class Matrix4x4;

class Vector4
{
    float _d[4] { 0, 0, 0, 0 };
public:
    Vector4(){}
    Vector4(std::initializer_list<float> list);
    explicit Vector4(const Vector3 &cartesian, float w = 1);
    explicit Vector4(float x, float y = 0, float z = 0, float w = 1);

    inline static Vector4 zero() { return Vector4(); }
    inline static Vector4 i() { return Vector4(1, 0, 0, 0); }
    inline static Vector4 j() { return Vector4(0, 1, 0, 0); }
    inline static Vector4 k() { return Vector4(0, 0, 1, 0); }
    inline static Vector4 h() { return Vector4(0, 0, 0, 1); }

    float *data();
    const float *data() const;

    const float &operator [](const int &i) const { return this->_d[i]; }
    float &operator [](const int &i) { return this->_d[i]; }

    const float &at(const int &i) const { return this->_d[i]; }
    float &at(const int &i) { return this->_d[i]; }

    inline const float &get(const int &i) const { return this->at(i); }
    inline float &get(const int &i) { return this->at(i); }

    const float &x() const { return this->_d[0]; }
    float &x() { return this->_d[0]; }
    const float &y() const { return this->_d[1]; }
    float &y() { return this->_d[1]; }
    const float &z() const { return this->_d[2]; }
    float &z() { return this->_d[2]; }
    const float &w() const { return this->_d[3]; }
    float &w() { return this->_d[3]; }

    Vector4 normalize() const;
    float magnitude() const;
    float sqrMagnitude() const;
    Vector4 proj(const Vector4 &p) const;
    Vector4 lerp(const Vector4 &to, float t) const;
    Vector4 slerp(const Vector4 &to, float t) const;
    float distance(const Vector4 &p) const;
    Vector3 toCartesian() const;

    std::string toString() const;
    Vector4 copy() const;

    static Vector4 add(const Vector4 &left, const Vector4 &right);
    static Vector4 sub(const Vector4 &left, const Vector4 &right);
    static float dot(const Vector4 &left, const Vector4 &right);
    static Vector4 mul(const Vector4 &left, const float &right);
    static Vector4 mul(const Matrix4x4 &left, const Vector4 &right);
    static Vector4 div(const Vector4 &left, const float &right);
    static Vector4 neg(const Vector4 &v);
    static bool equal(const Vector4 &left, const Vector4 &right);

    friend inline Vector4 operator -(const Vector4 &v) { return Vector4::neg(v); }

    friend inline Vector4 operator +(const Vector4 &v, const Vector4 &r) { return Vector4::add(v, r); }
    friend inline Vector4 operator -(const Vector4 &v, const Vector4 &r) { return Vector4::sub(v, r); }
    friend inline Vector4 operator *(const Vector4 &v, float f) { return Vector4::mul(v, f); }
    friend inline Vector4 operator *(float f, const Vector4 &v) { return Vector4::mul(v, f); }
    friend inline float operator *(const Vector4 &v, const Vector4 &w) { return Vector4::dot(v, w); }
    friend inline Vector4 operator *(const Matrix4x4 &v, const Vector4 &w) { return Vector4::mul(v, w); }
    friend inline Vector4 operator /(const Vector4 &v, float d) { return Vector4::div(v, d); }

    friend inline bool operator ==(const Vector4 &v, Vector4 d) { return Vector4::equal(v, d); }
    friend inline bool operator !=(const Vector4 &v, Vector4 d) { return !Vector4::equal(v, d); }
};
}

#endif // VECTOR4_H
