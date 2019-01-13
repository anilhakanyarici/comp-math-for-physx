#ifndef VECTOR3_H
#define VECTOR3_H

#include <string>
#include <initializer_list>

namespace mp
{

class Matrix3x3;
class Quaternion;

class Vector3
{
    float _d[3] { 0, 0, 0 };
public:
    Vector3(){}
    Vector3(std::initializer_list<float> list);
    explicit Vector3(float x, float y = 0, float z = 0);

    inline static Vector3 zero() { return Vector3(); }
    inline static Vector3 i() { return Vector3(1, 0, 0); }
    inline static Vector3 j() { return Vector3(0, 1, 0); }
    inline static Vector3 k() { return Vector3(0, 0, 1); }

    float *data();
    const float *data() const;

    const float &operator [](const int &i) const;
    float &operator [](const int &i);

    const float &at(const int &i) const;
    float &at(const int &i);

    inline const float &get(const int &i) const { return this->at(i); }
    inline float &get(const int &i) { return this->at(i); }

    const float &x() const { return this->_d[0]; }
    float &x() { return this->_d[0]; }
    const float &y() const { return this->_d[1]; }
    float &y() { return this->_d[1]; }
    const float &z() const { return this->_d[2]; }
    float &z() { return this->_d[2]; }

    Vector3 normalize() const;
    float magnitude() const;
    float sqrMagnitude() const;
    Vector3 proj(const Vector3 &p) const;
    Vector3 lerp(const Vector3 &to, float t) const;
    Vector3 slerp(const Vector3 &to, float t) const;
    float distance(const Vector3 &p) const;
    Vector3 rotate(const Quaternion &q) const;
    Vector3 rotate(float angle, const Vector3 &rotation_axis) const;

    std::string toString() const;
    Vector3 copy() const;

    static Vector3 add(const Vector3 &left, const Vector3 &right);
    static Vector3 sub(const Vector3 &left, const Vector3 &right);
    static float dot(const Vector3 &left, const Vector3 &right);
    static Vector3 cross(const Vector3 &left, const Vector3 &right);
    static Vector3 mul(const Vector3 &left, const float &right);
    static Vector3 mul(const Matrix3x3 &m, const Vector3 &v);
    static Vector3 div(const Vector3 &left, const float &right);
    static bool equal(const Vector3 &left, const Vector3 &right);

    friend inline Vector3 operator +(const Vector3 &v, const Vector3 &r) { return Vector3::add(v, r); }
    friend inline Vector3 operator -(const Vector3 &v, const Vector3 &r) { return Vector3::sub(v, r); }
    friend inline Vector3 operator *(const Vector3 &v, float f) { return Vector3::mul(v, f); }
    friend inline Vector3 operator *(float f, const Vector3 &v) { return Vector3::mul(v, f); }
    friend inline float operator *(const Vector3 &v, const Vector3 &w) { return Vector3::dot(v, w); }
    friend inline Vector3 operator *(const Matrix3x3 &m, const Vector3 &v) { return Vector3::mul(m, v); }
    friend inline Vector3 operator %(const Vector3 &v, const Vector3 &w) { return Vector3::cross(v, w); }
    friend inline Vector3 operator /(const Vector3 &v, float d) { return Vector3::div(v, d); }

    friend inline bool operator ==(const Vector3 &v, Vector3 d) { return Vector3::equal(v, d); }
    friend inline bool operator !=(const Vector3 &v, Vector3 d) { return !Vector3::equal(v, d); }
};
}

#endif // VECTOR3_H
