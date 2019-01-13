#ifndef VECTOR_H
#define VECTOR_H

#include <memory>

#include <QString>

namespace mp
{

class Vector3;
class Vector4;

class Vector
{
    struct pimpl;
    std::shared_ptr<pimpl> _pimpl;

public:
    Vector();
    Vector(const Vector3 &v);
    Vector(const Vector4 &v);
    Vector(std::initializer_list<float> list);
    explicit Vector(int size, bool zeros = false);
    Vector(const float *v, int size, bool copy = true);

    float *data();
    const float *data() const;

    const float &operator [](const int &i) const;
    float &operator [](const int &i);

    const float &at(const int &i) const;
    float &at(const int &i);

    inline const float &get(const int &i) const { return this->at(i); }
    inline float &get(const int &i) { return this->at(i); }

    int size() const;

    float magnitude() const;
    float sqrMagnitude() const;
    Vector normalize() const;
    Vector proj(const Vector &p);
    Vector lerp(const Vector &to, float t) const;
    Vector slerp(const Vector &to, float t) const;
    float distance(const Vector &p) const;

    std::string toString() const;
    Vector copy() const;

    static Vector add(const Vector &l, const Vector &r);
    static Vector sub(const Vector &l, const Vector &r);
    static float mul(const Vector &l, const Vector &r);
    static Vector mul(const Vector &l, const float &r);
    static Vector div(const Vector &l, const float &r);
    static bool equal(const Vector &l, const Vector &r);

    friend inline Vector operator +(const Vector &l, const Vector &r) { return Vector::add(l, r); }
    friend inline Vector operator -(const Vector &l, const Vector &r) { return Vector::sub(l, r); }
    friend inline float operator *(const Vector &l, const Vector &r) { return Vector::mul(l, r); }
    friend inline Vector operator *(const Vector &l, const float &r) { return Vector::mul(l, r); }
    friend inline Vector operator *(const float &l, const Vector &r) { return Vector::mul(r, l); }
    friend inline Vector operator /(const Vector &l, const float &r) { return Vector::div(l, r); }

    friend inline bool operator ==(const Vector &l, const Vector &r) { return Vector::equal(l, r); }
    friend inline bool operator !=(const Vector &l, const Vector &r) { return !Vector::equal(l, r); }
};
}

#endif // VECTOR_H
