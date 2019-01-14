#include "matrix3x3.h"

#include "../func.h"
#include "vector3.h"

#include <math.h>

using namespace mp;

Matrix3x3::Matrix3x3(float diags)
{
    this->_d = std::vector<Vector3>(3);
    this->_d[0][0] = diags;
    this->_d[0][1] = 0;
    this->_d[0][2] = 0;
    this->_d[1][0] = 0;
    this->_d[1][1] = diags;
    this->_d[1][2] = 0;
    this->_d[2][0] = 0;
    this->_d[2][1] = 0;
    this->_d[2][2] = diags;

}

Matrix3x3::Matrix3x3(std::initializer_list<Vector3> list)
{
    assert(list.size() <= 3);
    const Vector3 *d = list.begin();
    this->_d = std::vector<Vector3>(3);
    for(size_t i = 0; i < list.size(); ++i){
        this->_d[i][0] = d[0][i];
        this->_d[i][1] = d[1][i];
        this->_d[i][2] = d[2][i];
    }
}

Matrix3x3::Matrix3x3(const float &m00, const float &m10, const float &m20, const float &m01, const float &m11, const float &m21, const float &m02, const float &m12, const float &m22)
{
    this->_d = std::vector<Vector3>(3);
    this->_d[0][0] = m00;
    this->_d[0][1] = m01;
    this->_d[0][2] = m02;
    this->_d[1][0] = m10;
    this->_d[1][1] = m11;
    this->_d[1][2] = m12;
    this->_d[2][0] = m20;
    this->_d[2][1] = m21;
    this->_d[2][2] = m22;
}

Matrix3x3 Matrix3x3::random()
{
    Matrix3x3 m;
    m._d[0][0] = mp::random(-1.0, 1.0);
    m._d[0][1] = mp::random(-1.0, 1.0);
    m._d[0][2] = mp::random(-1.0, 1.0);
    m._d[1][0] = mp::random(-1.0, 1.0);
    m._d[1][1] = mp::random(-1.0, 1.0);
    m._d[1][2] = mp::random(-1.0, 1.0);
    m._d[2][0] = mp::random(-1.0, 1.0);
    m._d[2][1] = mp::random(-1.0, 1.0);
    m._d[2][2] = mp::random(-1.0, 1.0);
    return m;
}

Vector3 *Matrix3x3::data()
{
    return this->_d.data();
}

const Vector3 *Matrix3x3::data() const
{
    return this->_d.data();
}

const Vector3 &Matrix3x3::operator [](const int &col) const
{
    return this->_d[col];
}

Vector3 &Matrix3x3::operator [](const int &col)
{
    return this->_d[col];
}

const float &Matrix3x3::at(const int &col, const int &row) const
{
    return this->_d[col][row];
}

float &Matrix3x3::at(const int &col, const int &row)
{
    return this->_d[col][row];
}

float Matrix3x3::determinant() const
{
    float e00 = this->at(0, 0);
    float e01 = this->at(0, 1);
    float e02 = this->at(0, 2);

    float m00 = this->at(1, 1) * this->at(2, 2) - this->at(1, 2) * this->at(2, 1);
    float m01 = this->at(1, 2) * this->at(2, 0) - this->at(1, 0) * this->at(2, 2);
    float m02 = this->at(1, 0) * this->at(2, 1) - this->at(1, 1) * this->at(2, 0);

    return e00 * m00 + e01 * m01 + e02 * m02;
}

Matrix3x3 Matrix3x3::transpose() const
{
    Matrix3x3 m;
    m._d[0][0] = this->_d[0][0];
    m._d[0][1] = this->_d[1][0];
    m._d[0][2] = this->_d[2][0];
    m._d[1][0] = this->_d[0][1];
    m._d[1][1] = this->_d[1][1];
    m._d[1][2] = this->_d[2][1];
    m._d[2][0] = this->_d[0][2];
    m._d[2][1] = this->_d[1][2];
    m._d[2][2] = this->_d[2][2];
    return m;
}

Matrix3x3 Matrix3x3::inverse() const
{
    float m00 = this->_d[1][1] * this->_d[2][2] - this->_d[1][2] * this->_d[2][1];
    float m10 = this->_d[1][0] * this->_d[2][2] - this->_d[1][2] * this->_d[2][0];
    float m20 = this->_d[1][0] * this->_d[2][1] - this->_d[1][1] * this->_d[2][0];

    float m01 = this->_d[0][1] * this->_d[2][2] - this->_d[0][2] * this->_d[2][1];
    float m11 = this->_d[0][0] * this->_d[2][2] - this->_d[0][2] * this->_d[2][0];
    float m21 = this->_d[0][0] * this->_d[2][1] - this->_d[0][1] * this->_d[2][0];

    float m02 = this->_d[0][1] * this->_d[1][2] - this->_d[0][2] * this->_d[1][1];
    float m12 = this->_d[0][0] * this->_d[1][2] - this->_d[0][2] * this->_d[1][0];
    float m22 = this->_d[0][0] * this->_d[1][1] - this->_d[0][1] * this->_d[1][0];

    float d = this->_d[0][0] * m00 - this->_d[0][1] * m10 + this->_d[0][2] * m20;

    m01 = -m01 / d; m10 = -m10 / d; m12 = -m12 / d; m21 = -m21 / d;
    m00 = m00 / d; m02 = m02 / d; m11 = m11 / d; m20 = m20 / d; m22 = m22 / d;

    return Matrix3x3(m00, m10, m20, m01, m11, m21, m02, m12, m22);
}

Matrix3x3 Matrix3x3::minorMatrix() const
{
    Matrix3x3 mm;
    mm._d[0][0] = this->_d[1][1] * this->_d[2][2] - this->_d[2][1] * this->_d[1][2];
    mm._d[1][0] = this->_d[0][1] * this->_d[2][2] - this->_d[2][1] * this->_d[0][2];
    mm._d[2][0] = this->_d[0][1] * this->_d[1][2] - this->_d[1][1] * this->_d[0][2];
    mm._d[0][1] = this->_d[1][0] * this->_d[2][2] - this->_d[2][0] * this->_d[1][2];
    mm._d[1][1] = this->_d[0][0] * this->_d[2][2] - this->_d[2][0] * this->_d[0][2];
    mm._d[2][1] = this->_d[0][0] * this->_d[1][2] - this->_d[1][0] * this->_d[0][2];
    mm._d[0][2] = this->_d[1][0] * this->_d[2][1] - this->_d[2][0] * this->_d[1][1];
    mm._d[1][2] = this->_d[0][0] * this->_d[2][1] - this->_d[2][0] * this->_d[0][1];
    mm._d[2][2] = this->_d[0][0] * this->_d[1][1] - this->_d[1][0] * this->_d[0][1];
    return mm;
}

Matrix3x3 Matrix3x3::cofactorMatrix() const
{
    Matrix3x3 mm;
    mm._d[0][0] = this->_d[1][1] * this->_d[2][2] - this->_d[2][1] * this->_d[1][2];
    mm._d[1][0] = -(this->_d[0][1] * this->_d[2][2] - this->_d[2][1] * this->_d[0][2]);
    mm._d[2][0] = this->_d[0][1] * this->_d[1][2] - this->_d[1][1] * this->_d[0][2];
    mm._d[0][1] = -(this->_d[1][0] * this->_d[2][2] - this->_d[2][0] * this->_d[1][2]);
    mm._d[1][1] = this->_d[0][0] * this->_d[2][2] - this->_d[2][0] * this->_d[0][2];
    mm._d[2][1] = -(this->_d[0][0] * this->_d[1][2] - this->_d[1][0] * this->_d[0][2]);
    mm._d[0][2] = this->_d[1][0] * this->_d[2][1] - this->_d[2][0] * this->_d[1][1];
    mm._d[1][2] = -(this->_d[0][0] * this->_d[2][1] - this->_d[2][0] * this->_d[0][1]);
    mm._d[2][2] = this->_d[0][0] * this->_d[1][1] - this->_d[1][0] * this->_d[0][1];
    return mm;
}

Matrix3x3 Matrix3x3::adjoint() const
{
    Matrix3x3 mm;
    mm._d[0][0] = this->_d[1][1] * this->_d[2][2] - this->_d[2][1] * this->_d[1][2];
    mm._d[0][1] = -(this->_d[0][1] * this->_d[2][2] - this->_d[2][1] * this->_d[0][2]);
    mm._d[0][2] = this->_d[0][1] * this->_d[1][2] - this->_d[1][1] * this->_d[0][2];
    mm._d[1][0] = -(this->_d[1][0] * this->_d[2][2] - this->_d[2][0] * this->_d[1][2]);
    mm._d[1][1] = this->_d[0][0] * this->_d[2][2] - this->_d[2][0] * this->_d[0][2];
    mm._d[1][2] = -(this->_d[0][0] * this->_d[1][2] - this->_d[1][0] * this->_d[0][2]);
    mm._d[2][0] = this->_d[1][0] * this->_d[2][1] - this->_d[2][0] * this->_d[1][1];
    mm._d[2][1] = -(this->_d[0][0] * this->_d[2][1] - this->_d[2][0] * this->_d[0][1]);
    mm._d[2][2] = this->_d[0][0] * this->_d[1][1] - this->_d[1][0] * this->_d[0][1];
    return mm;
}

Matrix3x3 Matrix3x3::add(const Matrix3x3 &l, const Matrix3x3 &r)
{
    return Matrix3x3(
            l._d[0][0] + r._d[0][0],
            l._d[1][0] + r._d[1][0],
            l._d[2][0] + r._d[2][0],
            l._d[0][1] + r._d[0][1],
            l._d[1][1] + r._d[1][1],
            l._d[2][1] + r._d[2][1],
            l._d[0][2] + r._d[0][2],
            l._d[1][2] + r._d[1][2],
            l._d[2][2] + r._d[2][2]);
}

Matrix3x3 Matrix3x3::sub(const Matrix3x3 &l, const Matrix3x3 &r)
{
    return Matrix3x3(
            l._d[0][0] - r._d[0][0],
            l._d[1][0] - r._d[1][0],
            l._d[2][0] - r._d[2][0],
            l._d[0][1] - r._d[0][1],
            l._d[1][1] - r._d[1][1],
            l._d[2][1] - r._d[2][1],
            l._d[0][2] - r._d[0][2],
            l._d[1][2] - r._d[1][2],
            l._d[2][2] - r._d[2][2]);
}

Matrix3x3 Matrix3x3::mul(const Matrix3x3 &l, const Matrix3x3 &r)
{
    Matrix3x3 m;
    m.at(0, 0) = l.at(0, 0) * r.at(0, 0) + l.at(1, 0) * r.at(0, 1) + l.at(2, 0) * r.at(0, 2);
    m.at(1, 0) = l.at(0, 0) * r.at(1, 0) + l.at(1, 0) * r.at(1, 1) + l.at(2, 0) * r.at(1, 2);
    m.at(2, 0) = l.at(0, 0) * r.at(2, 0) + l.at(1, 0) * r.at(2, 1) + l.at(2, 0) * r.at(2, 2);
    m.at(0, 1) = l.at(0, 1) * r.at(0, 0) + l.at(1, 1) * r.at(0, 1) + l.at(2, 1) * r.at(0, 2);
    m.at(1, 1) = l.at(0, 1) * r.at(1, 0) + l.at(1, 1) * r.at(1, 1) + l.at(2, 1) * r.at(1, 2);
    m.at(2, 1) = l.at(0, 1) * r.at(2, 0) + l.at(1, 1) * r.at(2, 1) + l.at(2, 1) * r.at(2, 2);
    m.at(0, 2) = l.at(0, 2) * r.at(0, 0) + l.at(1, 2) * r.at(0, 1) + l.at(2, 2) * r.at(0, 2);
    m.at(1, 2) = l.at(0, 2) * r.at(1, 0) + l.at(1, 2) * r.at(1, 1) + l.at(2, 2) * r.at(1, 2);
    m.at(2, 2) = l.at(0, 2) * r.at(2, 0) + l.at(1, 2) * r.at(2, 1) + l.at(2, 2) * r.at(2, 2);
    return m;
}

Vector3 Matrix3x3::mul(const Matrix3x3 &m, const Vector3 &v)
{
    Vector3 r;
    r[0] = m.at(0, 0) * v[0] +  m.at(1, 0) * v[1] +  m.at(2, 0) * v[2];
    r[1] = m.at(0, 1) * v[0] +  m.at(1, 1) * v[1] +  m.at(2, 1) * v[2];
    r[2] = m.at(0, 2) * v[0] +  m.at(1, 2) * v[1] +  m.at(2, 2) * v[2];
    return r;
}

Matrix3x3 Matrix3x3::mul(const Matrix3x3 &m, const float &f)
{
    Matrix3x3 n;
    n._d[0][0] = m._d[0][0] * f;
    n._d[0][1] = m._d[0][1] * f;
    n._d[0][2] = m._d[0][2] * f;
    n._d[1][0] = m._d[1][0] * f;
    n._d[1][1] = m._d[1][1] * f;
    n._d[1][2] = m._d[1][2] * f;
    n._d[2][0] = m._d[2][0] * f;
    n._d[2][1] = m._d[2][1] * f;
    n._d[2][2] = m._d[2][2] * f;
    return n;
}

Matrix3x3 Matrix3x3::div(const Matrix3x3 &m, const float &f)
{
    Matrix3x3 n;
    n._d[0][0] = m._d[0][0] / f;
    n._d[0][1] = m._d[0][1] / f;
    n._d[0][2] = m._d[0][2] / f;
    n._d[1][0] = m._d[1][0] / f;
    n._d[1][1] = m._d[1][1] / f;
    n._d[1][2] = m._d[1][2] / f;
    n._d[2][0] = m._d[2][0] / f;
    n._d[2][1] = m._d[2][1] / f;
    n._d[2][2] = m._d[2][2] / f;
    return n;
}

Matrix3x3 Matrix3x3::neg(const Matrix3x3 &m)
{
    Matrix3x3 n;
    n._d[0][0] = -m._d[0][0];
    n._d[0][1] = -m._d[0][1];
    n._d[0][2] = -m._d[0][2];
    n._d[1][0] = -m._d[1][0];
    n._d[1][1] = -m._d[1][1];
    n._d[1][2] = -m._d[1][2];
    n._d[2][0] = -m._d[2][0];
    n._d[2][1] = -m._d[2][1];
    n._d[2][2] = -m._d[2][2];
    return n;
}

bool Matrix3x3::equal(const Matrix3x3 &l, const Matrix3x3 &r)
{
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            if(::fabs(l.at(i, j) - r.at(i, j)) >= 1e-15)
                return false;
        }
    }
    return true;
}

std::string Matrix3x3::toString() const
{
    std::string str;
    for(int i = 0; i < 3; ++i){
        std::string row_str = "";
        for(int j = 0; j < 3; ++j){
            row_str.append(std::to_string(this->at(j, i))).append("\t");
        }
        str.append(row_str).append("\n");
    }
    return str;
}

Matrix3x3 Matrix3x3::copy() const
{
    Matrix3x3 n;
    n._d[0][0] = this->_d[0][0];
    n._d[0][1] = this->_d[0][1];
    n._d[0][2] = this->_d[0][2];
    n._d[1][0] = this->_d[1][0];
    n._d[1][1] = this->_d[1][1];
    n._d[1][2] = this->_d[1][2];
    n._d[2][0] = this->_d[2][0];
    n._d[2][1] = this->_d[2][1];
    n._d[2][2] = this->_d[2][2];
    return n;
}
