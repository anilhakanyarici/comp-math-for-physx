#ifndef MATRIX4X4_H
#define MATRIX4X4_H

#include <string>
#include <vector>

#include "vector4.h"

namespace mp
{

class Vector4;

class Matrix4x4
{
    std::vector<Vector4> _d;
public:
    Matrix4x4(float diags = 0);
    Matrix4x4(std::initializer_list<Vector4> list);

    static inline Matrix4x4 identity() { return Matrix4x4(1); }

    Vector4 *data();
    const Vector4 *data() const;

    const Vector4 &operator [](const int &col) const;
    Vector4 &operator [](const int &col);

    const float &at(const int &col, const int &row) const;
    float &at(const int &col, const int &row);

    inline const float &get(const int &row, const int &col) const { return this->at(col, row); }
    inline float &get(const int &row, const int &col) { return this->at(col, row); }

    float determinant() const;
    Matrix4x4 transpose() const;
    Matrix4x4 inverse() const;

    static Matrix4x4 add(const Matrix4x4 &l, const Matrix4x4 &r);
    static Matrix4x4 sub(const Matrix4x4 &l, const Matrix4x4 &r);
    static Matrix4x4 mul(const Matrix4x4 &l, const Matrix4x4 &r);
    static Vector4 mul(const Matrix4x4 &m, const Vector4 &v);
    static Matrix4x4 mul(const Matrix4x4 &m, const float &f);
    static Matrix4x4 div(const Matrix4x4 &m, const float &f);
    static Matrix4x4 neg(const Matrix4x4 &m);
    static bool equal(const Matrix4x4 &l, const Matrix4x4 &r);

    std::string toString() const;
    Matrix4x4 copy() const;

    friend inline Matrix4x4 operator -(const Matrix4x4 &v) { return Matrix4x4::neg(v); }

    friend Matrix4x4 operator +(const Matrix4x4 &m, const Matrix4x4 &r) { return Matrix4x4::add(m, r); }
    friend Matrix4x4 operator -(const Matrix4x4 &m, const Matrix4x4 &r) { return Matrix4x4::sub(m, r); }
    friend Matrix4x4 operator *(const Matrix4x4 &m, const Matrix4x4 &r) { return Matrix4x4::mul(m, r); }
    friend Matrix4x4 operator *(const Matrix4x4 &m, const float &f) { return Matrix4x4::mul(m, f); }
    friend Matrix4x4 operator /(const Matrix4x4 &m, const float &f) { return Matrix4x4::div(m, f); }

    friend Matrix4x4 operator ==(const Matrix4x4 &m, const Matrix4x4 &r) { return Matrix4x4::equal(m, r); }
    friend Matrix4x4 operator !=(const Matrix4x4 &m, const Matrix4x4 &r) { return !Matrix4x4::equal(m, r); }
};
}

#endif // MATRIX4X4_H
