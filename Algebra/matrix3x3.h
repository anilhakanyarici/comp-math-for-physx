#ifndef MATRIX3X3_H
#define MATRIX3X3_H

#include <string>
#include <vector>

namespace mp {

class Vector3;

class Matrix3x3
{
    std::vector<Vector3> _d;
public:
    Matrix3x3(float diags = 0);
    Matrix3x3(std::initializer_list<Vector3> list);
    Matrix3x3(const float &m00, const float &m10, const float &m20, const float &m01, const float &m11, const float &m21, const float &m02, const float &m12, const float &m22);

    static Matrix3x3 identity() { return Matrix3x3(1); }
    static Matrix3x3 random();

    Vector3 *data();
    const Vector3 *data() const;

    const Vector3 &operator [](const int &col) const;
    Vector3 &operator [](const int &col);

    const float &at(const int &col, const int &row) const;
    float &at(const int &col, const int &row);

    inline const float &get(const int &row, const int &col) const { return this->at(col, row); }
    inline float &get(const int &row, const int &col) { return this->at(col, row); }

    float determinant() const;
    Matrix3x3 transpose() const;
    Matrix3x3 inverse() const;
    Matrix3x3 minorMatrix() const;
    Matrix3x3 cofactorMatrix() const;
    Matrix3x3 adjoint() const;

    static Matrix3x3 add(const Matrix3x3 &l, const Matrix3x3 &r);
    static Matrix3x3 sub(const Matrix3x3 &l, const Matrix3x3 &r);
    static Matrix3x3 mul(const Matrix3x3 &l, const Matrix3x3 &r);
    static Vector3 mul(const Matrix3x3 &m, const Vector3 &v);
    static Matrix3x3 mul(const Matrix3x3 &m, const float &f);
    static Matrix3x3 div(const Matrix3x3 &m, const float &f);
    static Matrix3x3 neg(const Matrix3x3 &m);
    static bool equal(const Matrix3x3 &l, const Matrix3x3 &r);

    std::string toString() const;
    Matrix3x3 copy() const;

    friend inline Matrix3x3 operator -(const Matrix3x3 &v) { return Matrix3x3::neg(v); }

    friend Matrix3x3 operator +(const Matrix3x3 &m, const Matrix3x3 &r) { return Matrix3x3::add(m, r); }
    friend Matrix3x3 operator -(const Matrix3x3 &m, const Matrix3x3 &r) { return Matrix3x3::sub(m, r); }
    friend Matrix3x3 operator *(const Matrix3x3 &m, const Matrix3x3 &r) { return Matrix3x3::mul(m, r); }
    friend Matrix3x3 operator *(const Matrix3x3 &m, const float &f) { return Matrix3x3::mul(m, f); }
    friend Matrix3x3 operator /(const Matrix3x3 &m, const float &f) { return Matrix3x3::div(m, f); }

    friend Matrix3x3 operator ==(const Matrix3x3 &m, const Matrix3x3 &r) { return Matrix3x3::equal(m, r); }
    friend Matrix3x3 operator !=(const Matrix3x3 &m, const Matrix3x3 &r) { return !Matrix3x3::equal(m, r); }
};
}

#endif // MATRIX3X3_H
