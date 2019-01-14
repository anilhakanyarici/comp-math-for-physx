#ifndef QUATERNION_H
#define QUATERNION_H

#include <QString>

namespace mp
{

class Vector3;
class Matrix3x3;

class Quaternion
{
    float _d[4];
public:
    Quaternion(){}
    Quaternion(const float &w);
    Quaternion(const Vector3 &v);
    Quaternion(float w, float x, float y, float z);
    Quaternion(float angle, const Vector3 &axis);

    const float &w() const { return this->_d[0]; }
    float &w() { return this->_d[0]; }
    const float &x() const { return this->_d[1]; }
    float &x() { return this->_d[1]; }
    const float &y() const { return this->_d[2]; }
    float &y() { return this->_d[2]; }
    const float &z() const { return this->_d[3]; }
    float &z() { return this->_d[3]; }

    const float &operator [](const int &i) const;
    float &operator [](const int &i);

    const float &at(const int &i) const;
    float &at(const int &i);

    inline const float &get(const int &i) const { return this->at(i); }
    inline float &get(const int &i) { return this->at(i); }

    float scalar() const;
    Vector3 vector() const;
    void setScalar(float w);
    void setVector(const Vector3 &v);

    float angle() const;
    Vector3 axis() const;
    Quaternion normalize() const;
    float magnitude() const;
    float sqrMagnitude() const;
    Quaternion conjugate() const;
    Quaternion reciprocal() const; //multiplicative inverse

    std::string toString() const;
    Quaternion copy() const;

    Vector3 toEuler() const;
    Matrix3x3 toRotationMatrix() const;

    static Quaternion euler(float pitch, float roll, float yaw);
    static Quaternion euler(const Vector3 &euler_angles);
    static Quaternion identity();

    static Quaternion add(const Quaternion &p, const Quaternion &q);
    static Quaternion sub(const Quaternion &p, const Quaternion &q);
    static Quaternion mul(const Quaternion &p, const Quaternion &q);
    static Quaternion mul(const Vector3 &v, const Quaternion &q);
    static Quaternion mul(const Quaternion &p, const Vector3 &v);
    static Quaternion mul(const Quaternion &p, float quot);
    static Quaternion div(const Quaternion &p, const Quaternion &q);
    static Quaternion div(const Quaternion &p, float den);
    static Quaternion neg(const Quaternion &q);
    static bool equal(const Quaternion &p, const Quaternion &q);

    friend inline Quaternion operator -(const Quaternion &q) { return Quaternion::neg(q); }

    friend inline Quaternion operator +(const Quaternion &q, const Quaternion &p) { return Quaternion::add(q, p); }
    friend inline Quaternion operator +(const Quaternion &q, const Vector3 &v) { return Quaternion::add(q, v); }
    friend inline Quaternion operator +(const Vector3 &v, const Quaternion &q) { return Quaternion::add(v, q); }
    friend inline Quaternion operator +(const Quaternion &q, const float &f) { return Quaternion::add(q, f); }
    friend inline Quaternion operator +(const float &f, const Quaternion &q) { return Quaternion::add(f, q); }
    friend inline Quaternion operator +(const Vector3 &v, const float &f) { return Quaternion::add(v, f); }
    friend inline Quaternion operator +(const float &f, const Vector3 &v) { return Quaternion::add(f, v); }

    friend inline Quaternion operator -(const Quaternion &q, const Quaternion &p) { return Quaternion::sub(q, p); }
    friend inline Quaternion operator -(const Quaternion &q, const Vector3 &v) { return Quaternion::sub(q, v); }
    friend inline Quaternion operator -(const Vector3 &v, const Quaternion &q) { return Quaternion::sub(v, q); }
    friend inline Quaternion operator -(const Quaternion &q, const float &f) { return Quaternion::sub(q, f); }
    friend inline Quaternion operator -(const float &f, const Quaternion &q) { return Quaternion::sub(f, q); }
    friend inline Quaternion operator -(const Vector3 &v, const float &f) { return Quaternion::sub(v, f); }
    friend inline Quaternion operator -(const float &f, const Vector3 &v) { return Quaternion::sub(f, v); }

    friend inline Quaternion operator *(const float &f, const Quaternion &q) { return Quaternion::mul(q, f); }
    friend inline Quaternion operator *(const Quaternion &q, const float &f) { return Quaternion::mul(q, f); }
    friend inline Quaternion operator *(const Quaternion &q, const Vector3 &v) { return Quaternion::mul(q, v); }
    friend inline Quaternion operator *(const Vector3 &v, const Quaternion &q) { return Quaternion::mul(v, q); }
    friend inline Quaternion operator *(const Quaternion &q, const Quaternion &p) { return Quaternion::mul(q, p); }

    friend inline Quaternion operator /(const Quaternion &q, const float &w) { return Quaternion::mul(q, w); }
    friend inline Quaternion operator /(const float &w, const Quaternion &q) { return Quaternion::mul(q.reciprocal(), w); }
    friend inline Quaternion operator /(const Vector3 &v, const Quaternion &q) { return Quaternion::mul(v, q.reciprocal()); }
    friend inline Quaternion operator /(const Quaternion &p, const Quaternion &q) { return Quaternion::mul(p, q.reciprocal()); }

    friend inline Quaternion operator ==(const Quaternion &q, const Quaternion &p) { return Quaternion::equal(q, p); }
    friend inline Quaternion operator ==(const Vector3 &q, const Quaternion &p) { return Quaternion::equal(q, p); }
    friend inline Quaternion operator ==(const float &q, const Quaternion &p) { return Quaternion::equal(q, p); }
    friend inline Quaternion operator !=(const Quaternion &q, const Quaternion &p) { return !Quaternion::equal(q, p); }
    friend inline Quaternion operator !=(const Vector3 &q, const Quaternion &p) { return !Quaternion::equal(q, p); }
    friend inline Quaternion operator !=(const float &q, const Quaternion &p) { return !Quaternion::equal(q, p); }
};
}

#endif // QUATERNION_H
