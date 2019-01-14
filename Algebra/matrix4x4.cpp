#include "matrix4x4.h"

#include <math.h>

#include "vector4.h"
#include "matrix.h"

using namespace mp;

Matrix4x4::Matrix4x4(float diags)
{
    this->_d = std::vector<Vector4>(4);
    this->_d[0][0] = diags;
    this->_d[0][1] = 0;
    this->_d[0][2] = 0;
    this->_d[0][3] = 0;
    this->_d[1][0] = 0;
    this->_d[1][1] = diags;
    this->_d[1][2] = 0;
    this->_d[1][3] = 0;
    this->_d[2][2] = 0;
    this->_d[2][1] = 0;
    this->_d[2][2] = diags;
    this->_d[2][3] = 0;
    this->_d[3][0] = 0;
    this->_d[3][1] = 0;
    this->_d[3][2] = 0;
    this->_d[3][3] = diags;
}

Matrix4x4::Matrix4x4(std::initializer_list<Vector4> list)
{
    assert(list.size() <= 4);
    const Vector4 *d = list.begin();
    this->_d = std::vector<Vector4>(4);
    for(size_t i = 0; i < list.size(); ++i){
        this->_d[i][0] = d[0][i];
        this->_d[i][1] = d[1][i];
        this->_d[i][2] = d[2][i];
        this->_d[i][3] = d[3][i];
    }
}

Vector4 *Matrix4x4::data()
{
    return this->_d.data();
}

const Vector4 *Matrix4x4::data() const
{
    return this->_d.data();
}

const Vector4 &Matrix4x4::operator [](const int &col) const
{
    return this->_d[col];
}

Vector4 &Matrix4x4::operator [](const int &col)
{
    return this->_d[col];
}

const float &Matrix4x4::at(const int &col, const int &row) const
{
    return this->_d[col][row];
}

float &Matrix4x4::at(const int &col, const int &row)
{
    return this->_d[col][row];
}

float Matrix4x4::determinant() const
{
    return Matrix(*this).determinant();
}

Matrix4x4 Matrix4x4::transpose() const
{
    Matrix4x4 m;
    m._d[0][0] = this->_d[0][0];
    m._d[0][1] = this->_d[1][0];
    m._d[0][2] = this->_d[2][0];
    m._d[0][3] = this->_d[3][0];
    m._d[1][0] = this->_d[0][1];
    m._d[1][1] = this->_d[1][1];
    m._d[1][2] = this->_d[2][1];
    m._d[1][3] = this->_d[3][1];
    m._d[2][0] = this->_d[0][2];
    m._d[2][1] = this->_d[1][2];
    m._d[2][2] = this->_d[2][2];
    m._d[2][3] = this->_d[3][2];
    m._d[3][0] = this->_d[0][3];
    m._d[3][1] = this->_d[1][3];
    m._d[3][2] = this->_d[2][3];
    m._d[3][3] = this->_d[3][3];
    return m;
}

Matrix4x4 Matrix4x4::inverse() const
{
    Matrix i = Matrix(*this).inverse();
    Matrix4x4 m;
    m._d[0][0] = i[0][0];
    m._d[0][1] = i[0][1];
    m._d[0][2] = i[0][2];
    m._d[0][3] = i[0][3];
    m._d[1][0] = i[1][0];
    m._d[1][1] = i[1][1];
    m._d[1][2] = i[1][2];
    m._d[1][3] = i[1][3];
    m._d[2][0] = i[2][0];
    m._d[2][1] = i[2][1];
    m._d[2][2] = i[2][2];
    m._d[2][3] = i[2][3];
    m._d[3][0] = i[3][0];
    m._d[3][1] = i[3][1];
    m._d[3][2] = i[3][2];
    m._d[3][3] = i[3][3];
    return m;
}

Matrix4x4 Matrix4x4::add(const Matrix4x4 &l, const Matrix4x4 &r)
{
    return Matrix4x4({ { l._d[0][0] + r._d[0][0], l._d[0][1] + r._d[0][1], l._d[0][2] + r._d[0][2], l._d[0][3] + r._d[0][3] },
                       { l._d[1][0] + r._d[1][0], l._d[1][1] + r._d[1][1], l._d[1][2] + r._d[1][2], l._d[1][3] + r._d[1][3] },
                       { l._d[2][0] + r._d[2][0], l._d[2][1] + r._d[2][1], l._d[2][2] + r._d[2][2], l._d[2][3] + r._d[2][3] },
                       { l._d[3][0] + r._d[3][0], l._d[3][1] + r._d[3][1], l._d[3][2] + r._d[3][2], l._d[3][3] + r._d[3][3] } });
}

Matrix4x4 Matrix4x4::sub(const Matrix4x4 &l, const Matrix4x4 &r)
{
    return Matrix4x4({ { l._d[0][0] - r._d[0][0], l._d[0][1] - r._d[0][1], l._d[0][2] - r._d[0][2], l._d[0][3] - r._d[0][3] },
                       { l._d[1][0] - r._d[1][0], l._d[1][1] - r._d[1][1], l._d[1][2] - r._d[1][2], l._d[1][3] - r._d[1][3] },
                       { l._d[2][0] - r._d[2][0], l._d[2][1] - r._d[2][1], l._d[2][2] - r._d[2][2], l._d[2][3] - r._d[2][3] },
                       { l._d[3][0] - r._d[3][0], l._d[3][1] - r._d[3][1], l._d[3][2] - r._d[3][2], l._d[3][3] - r._d[3][3] } });
}

Matrix4x4 Matrix4x4::mul(const Matrix4x4 &l, const Matrix4x4 &r)
{
    Matrix4x4 m;
    m.at(0, 0) = l.at(0, 0) * r.at(0, 0) + l.at(1, 0) * r.at(0, 1) + l.at(2, 0) * r.at(0, 2) + l.at(3, 0) * r.at(0, 3);
    m.at(1, 0) = l.at(0, 0) * r.at(1, 0) + l.at(1, 0) * r.at(1, 1) + l.at(2, 0) * r.at(1, 2) + l.at(3, 0) * r.at(1, 3);
    m.at(2, 0) = l.at(0, 0) * r.at(2, 0) + l.at(1, 0) * r.at(2, 1) + l.at(2, 0) * r.at(2, 2) + l.at(3, 0) * r.at(2, 3);
    m.at(3, 0) = l.at(0, 0) * r.at(3, 0) + l.at(1, 0) * r.at(3, 1) + l.at(2, 0) * r.at(3, 2) + l.at(3, 0) * r.at(3, 3);
    m.at(0, 1) = l.at(0, 1) * r.at(0, 0) + l.at(1, 1) * r.at(0, 1) + l.at(2, 1) * r.at(0, 2) + l.at(3, 1) * r.at(0, 3);
    m.at(1, 1) = l.at(0, 1) * r.at(1, 0) + l.at(1, 1) * r.at(1, 1) + l.at(2, 1) * r.at(1, 2) + l.at(3, 1) * r.at(1, 3);
    m.at(2, 1) = l.at(0, 1) * r.at(2, 0) + l.at(1, 1) * r.at(2, 1) + l.at(2, 1) * r.at(2, 2) + l.at(3, 1) * r.at(2, 3);
    m.at(3, 1) = l.at(0, 1) * r.at(3, 0) + l.at(1, 1) * r.at(3, 1) + l.at(2, 1) * r.at(3, 2) + l.at(3, 1) * r.at(3, 3);
    m.at(0, 2) = l.at(0, 2) * r.at(0, 0) + l.at(1, 2) * r.at(0, 1) + l.at(2, 2) * r.at(0, 2) + l.at(3, 2) * r.at(0, 3);
    m.at(1, 2) = l.at(0, 2) * r.at(1, 0) + l.at(1, 2) * r.at(1, 1) + l.at(2, 2) * r.at(1, 2) + l.at(3, 2) * r.at(1, 3);
    m.at(2, 2) = l.at(0, 2) * r.at(2, 0) + l.at(1, 2) * r.at(2, 1) + l.at(2, 2) * r.at(2, 2) + l.at(3, 2) * r.at(2, 3);
    m.at(3, 2) = l.at(0, 2) * r.at(3, 0) + l.at(1, 2) * r.at(3, 1) + l.at(2, 2) * r.at(3, 2) + l.at(3, 2) * r.at(3, 3);
    m.at(0, 3) = l.at(0, 3) * r.at(0, 0) + l.at(1, 3) * r.at(0, 1) + l.at(2, 3) * r.at(0, 2) + l.at(3, 3) * r.at(0, 3);
    m.at(1, 3) = l.at(0, 3) * r.at(1, 0) + l.at(1, 3) * r.at(1, 1) + l.at(2, 3) * r.at(1, 2) + l.at(3, 3) * r.at(1, 3);
    m.at(2, 3) = l.at(0, 3) * r.at(2, 0) + l.at(1, 3) * r.at(2, 1) + l.at(2, 3) * r.at(2, 2) + l.at(3, 3) * r.at(2, 3);
    m.at(3, 3) = l.at(0, 3) * r.at(3, 0) + l.at(1, 3) * r.at(3, 1) + l.at(2, 3) * r.at(3, 2) + l.at(3, 3) * r.at(3, 3);
    return m;
}

Vector4 Matrix4x4::mul(const Matrix4x4 &m, const Vector4 &v)
{
    Vector4 r;
    r[0] = m.at(0, 0) * v[0] +  m.at(1, 0) * v[1] +  m.at(2, 0) * v[2] + m.at(3, 0) * v[3];
    r[1] = m.at(0, 1) * v[0] +  m.at(1, 1) * v[1] +  m.at(2, 1) * v[2] + m.at(3, 1) * v[3];
    r[2] = m.at(0, 2) * v[0] +  m.at(1, 2) * v[1] +  m.at(2, 2) * v[2] + m.at(3, 2) * v[3];
    r[3] = m.at(0, 3) * v[0] +  m.at(1, 3) * v[1] +  m.at(2, 3) * v[2] + m.at(3, 2) * v[3];
    return r;
}

Matrix4x4 Matrix4x4::mul(const Matrix4x4 &m, const float &f)
{
    Matrix4x4 n;
    n._d[0][0] = m._d[0][0] * f;
    n._d[0][1] = m._d[0][1] * f;
    n._d[0][2] = m._d[0][2] * f;
    n._d[0][3] = m._d[0][3] * f;
    n._d[1][0] = m._d[1][0] * f;
    n._d[1][1] = m._d[1][1] * f;
    n._d[1][2] = m._d[1][2] * f;
    n._d[1][3] = m._d[1][3] * f;
    n._d[2][0] = m._d[2][0] * f;
    n._d[2][1] = m._d[2][1] * f;
    n._d[2][2] = m._d[2][2] * f;
    n._d[2][3] = m._d[2][3] * f;
    n._d[3][0] = m._d[3][0] * f;
    n._d[3][1] = m._d[3][1] * f;
    n._d[3][2] = m._d[3][2] * f;
    n._d[3][3] = m._d[3][3] * f;
    return n;
}

Matrix4x4 Matrix4x4::div(const Matrix4x4 &m, const float &f)
{
    Matrix4x4 n;
    n._d[0][0] = m._d[0][0] / f;
    n._d[0][1] = m._d[0][1] / f;
    n._d[0][2] = m._d[0][2] / f;
    n._d[0][3] = m._d[0][3] / f;
    n._d[1][0] = m._d[1][0] / f;
    n._d[1][1] = m._d[1][1] / f;
    n._d[1][2] = m._d[1][2] / f;
    n._d[1][3] = m._d[1][3] / f;
    n._d[2][0] = m._d[2][0] / f;
    n._d[2][1] = m._d[2][1] / f;
    n._d[2][2] = m._d[2][2] / f;
    n._d[2][3] = m._d[2][3] / f;
    n._d[3][0] = m._d[3][0] / f;
    n._d[3][1] = m._d[3][1] / f;
    n._d[3][2] = m._d[3][2] / f;
    n._d[3][3] = m._d[3][3] / f;
    return n;
}

Matrix4x4 Matrix4x4::neg(const Matrix4x4 &m)
{
    Matrix4x4 n;
    n._d[0][0] = -m._d[0][0];
    n._d[0][1] = -m._d[0][1];
    n._d[0][2] = -m._d[0][2];
    n._d[0][3] = -m._d[0][3];
    n._d[1][0] = -m._d[1][0];
    n._d[1][1] = -m._d[1][1];
    n._d[1][2] = -m._d[1][2];
    n._d[1][3] = -m._d[1][3];
    n._d[2][0] = -m._d[2][0];
    n._d[2][1] = -m._d[2][1];
    n._d[2][2] = -m._d[2][2];
    n._d[2][3] = -m._d[2][3];
    n._d[3][0] = -m._d[3][0];
    n._d[3][1] = -m._d[3][1];
    n._d[3][2] = -m._d[3][2];
    n._d[3][3] = -m._d[3][3];
    return n;
}

bool Matrix4x4::equal(const Matrix4x4 &l, const Matrix4x4 &r)
{
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            if(::fabs(l.at(i, j) - r.at(i, j)) >= 1e-15)
                return false;
        }
    }
    return true;
}

std::string Matrix4x4::toString() const
{
    std::string str;
    for(int i = 0; i < 4; ++i){
        std::string row_str = "";
        for(int j = 0; j < 4; ++j){
            row_str.append(std::to_string(this->at(j, i))).append("\t");
        }
        str.append(row_str).append("\n");
    }
    return str;
}

Matrix4x4 Matrix4x4::copy() const
{
    Matrix4x4 n;
    n._d[0][0] = this->_d[0][0];
    n._d[0][1] = this->_d[0][1];
    n._d[0][2] = this->_d[0][2];
    n._d[0][3] = this->_d[0][3];
    n._d[1][0] = this->_d[1][0];
    n._d[1][1] = this->_d[1][1];
    n._d[1][2] = this->_d[1][2];
    n._d[1][3] = this->_d[1][3];
    n._d[2][0] = this->_d[2][0];
    n._d[2][1] = this->_d[2][1];
    n._d[2][2] = this->_d[2][2];
    n._d[2][3] = this->_d[2][3];
    n._d[3][0] = this->_d[3][0];
    n._d[3][1] = this->_d[3][1];
    n._d[3][2] = this->_d[3][2];
    n._d[3][3] = this->_d[3][3];
    return n;
}
