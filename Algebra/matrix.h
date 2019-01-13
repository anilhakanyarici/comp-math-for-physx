#ifndef MATRIX_H
#define MATRIX_H

#include "Algebra/vector.h"

namespace mp
{

class Matrix3x3;
class Matrix4x4;

class Matrix
{
    struct pimpl;
    std::shared_ptr<pimpl> _pimpl;

public:
    Matrix();
    Matrix(const Matrix3x3 &m);
    Matrix(const Matrix4x4 &m);
    Matrix(float v, int dimension = 1);
    Matrix(std::initializer_list<Vector> vs);
    Matrix(int column_size, int row_size, bool zeros = false);

    static Matrix identity(int dimension = 1);
    static Matrix random(int column_size, int row_size);
    static Matrix diagonal(float v, int dimension = 1);

    Vector *data();
    const Vector *data() const;

    const Vector &operator [](const int &col) const;
    Vector &operator [](const int &col);

    const float &at(const int &col, const int &row) const;
    float &at(const int &col, const int &row);

    inline const float &get(const int &row, const int &col) const { return this->at(col, row); }
    inline float &get(const int &row, const int &col) { return this->at(col, row); }

    int rowSize() const;
    int columnSize() const;
    bool isSquare() const;

    float determinant() const;
    float minorValue(const int &col, const int &row) const;
    Matrix minorMatrix() const;
    float cofactorValue(const int &col, const int &row) const;
    Matrix cofactorMatrix() const;
    Matrix transpose() const;
    Matrix adjoint() const;
    Matrix inverse() const;
    float trace() const;

    Matrix deleteColumn(const int &col);
    Matrix deleteRow(const int &row);
    Matrix deleteColumnAndRow(const int &col, const int &row) const;

    static Matrix add(const Matrix &l, const Matrix &r);
    static Matrix sub(const Matrix &l, const Matrix &r);
    static Vector mul(const Matrix &m, const Vector &v);
    static Matrix mul(const Matrix &l, const float &r);
    static Matrix mul(const Matrix &l, const Matrix &r);
    static Matrix div(const Matrix &l, const float &r);
    static bool equal(const Matrix &l, const Matrix &r);

    std::string toString() const;
    Matrix copy() const;

    friend inline Matrix operator +(const Matrix &l, const Matrix &r) { return Matrix::add(l, r); }
    friend inline Matrix operator -(const Matrix &l, const Matrix &r) { return Matrix::sub(l, r); }
    friend inline Matrix operator *(const Matrix &l, const Matrix &r) { return Matrix::mul(l, r); }
    friend inline Vector operator *(const Matrix &l, const Vector &r) { return Matrix::mul(l, r); }
    friend inline Vector operator *(const Vector &l, const Matrix &r) { return Matrix::mul(r, l); }
    friend inline Matrix operator *(const Matrix &l, const float &r) { return Matrix::mul(l, r); }
    friend inline Matrix operator *(const float &l, const Matrix &r) { return Matrix::mul(r, l); }
    friend inline Matrix operator /(const Matrix &l, const float &r) { return Matrix::div(l, r); }
    friend inline Matrix operator /(const Matrix &l, const Matrix &r) { return Matrix::mul(l, r.inverse()); }

    friend inline bool operator ==(const Matrix &l, const Matrix &r) { return Matrix::equal(l, r); }
    friend inline bool operator !=(const Matrix &l, const Matrix &r) { return !Matrix::equal(l, r); }
};
}

#endif // MATRIX_H
