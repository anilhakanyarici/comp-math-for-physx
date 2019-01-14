#include "matrix.h"

#include <math.h>
#include <vector>

#include "../func.h"
#include "Algebra/vector3.h"
#include "Algebra/matrix3x3.h"
#include "Algebra/matrix4x4.h"

using namespace mp;

struct Matrix::pimpl
{
    std::vector<Vector> data;
    int row_size = 0;
    int column_size = 0;
};

Matrix::Matrix()
{
    this->_pimpl = std::shared_ptr<pimpl>(new Matrix::pimpl());
}

Matrix::Matrix(const Matrix3x3 &m)
{
    this->_pimpl = std::shared_ptr<pimpl>(new Matrix::pimpl());
    this->_pimpl->row_size = 3;
    this->_pimpl->column_size = 3;

    this->_pimpl->data = std::vector<Vector>(3);
    for(int i = 0; i < 3; ++i)
        this->_pimpl->data[i] = Vector(3, false);

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            this->_pimpl->data[i][j] = m.at(i, j);
        }
    }
}

Matrix::Matrix(const Matrix4x4 &m)
{
    this->_pimpl = std::shared_ptr<pimpl>(new Matrix::pimpl());
    this->_pimpl->row_size = 4;
    this->_pimpl->column_size = 4;

    this->_pimpl->data = std::vector<Vector>(4);
    for(int i = 0; i < 4; ++i)
        this->_pimpl->data[i] = Vector(4, false);

    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            this->_pimpl->data[i][j] = m.at(i, j);
        }
    }
}

Matrix::Matrix(float v, int dimension)
{
    this->_pimpl = std::shared_ptr<pimpl>(new Matrix::pimpl());
    if(dimension == 0)
        return;
    this->_pimpl->row_size = dimension;
    this->_pimpl->column_size = dimension;

    this->_pimpl->data = std::vector<Vector>(dimension);
    for(int i = 0; i < dimension; ++i){
        this->_pimpl->data[i] = Vector(dimension, true);
        this->at(i, i) = v;
    }
}

Matrix::Matrix(std::initializer_list<Vector> vs)
{
    const Vector *cols = vs.begin();
    int column_size = cols[0].size(), row_size = vs.size();
    this->_pimpl = std::shared_ptr<pimpl>(new Matrix::pimpl());
    if(column_size == 0 || row_size == 0)
        return;
    this->_pimpl->row_size = row_size;
    this->_pimpl->column_size = column_size;

    this->_pimpl->data = std::vector<Vector>(column_size);
    for(int i = 0; i < row_size; ++i){
        assert(cols[i].size() == column_size);
        this->_pimpl->data[i] = Vector(column_size);
        for(int j = 0; j < column_size; ++j){
            this->_pimpl->data[i][j] = cols[j][i];
        }
    }
}

Matrix::Matrix(int column_size, int row_size, bool zeros)
{
    this->_pimpl = std::shared_ptr<pimpl>(new Matrix::pimpl());
    if(column_size == 0 || row_size == 0)
        return;
    this->_pimpl->row_size = row_size;
    this->_pimpl->column_size = column_size;

    this->_pimpl->data = std::vector<Vector>(column_size);
    for(int i = 0; i < column_size; ++i)
        this->_pimpl->data[i] = Vector(row_size, zeros);
}

Matrix Matrix::identity(int dimension)
{
    return Matrix(1.0f, dimension);
}

Matrix Matrix::random(int column_size, int row_size)
{
    Matrix m(column_size, row_size);
    for(int i = 0; i < column_size; ++i){
        for(int j = 0; j < row_size; ++j){
            m.at(i, j) = mp::random(-1.0, 1.0);
        }
    }
    return m;
}

Matrix Matrix::diagonal(float v, int dimension)
{
    return Matrix(v, dimension);
}

Vector *Matrix::data()
{
    return this->_pimpl->data.data();
}

const Vector *Matrix::data() const
{
    return this->_pimpl->data.data();
}

const Vector &Matrix::operator [](const int &col) const
{
    return this->_pimpl->data[col];
}

Vector &Matrix::operator [](const int &col)
{
    return this->_pimpl->data[col];
}

const float &Matrix::at(const int &col, const int &row) const
{
    return this->_pimpl->data[col][row];
}

float &Matrix::at(const int &col, const int &row)
{
    return this->_pimpl->data[col][row];
}

int Matrix::rowSize() const
{
    return this->_pimpl->row_size;
}

int Matrix::columnSize() const
{
    return this->_pimpl->column_size;
}

bool Matrix::isSquare() const
{
    return this->_pimpl->row_size == this->_pimpl->column_size;
}

float Matrix::determinant() const
{
    if(this->_pimpl->column_size != this->_pimpl->row_size)
        return ::sqrt(Matrix::mul(*this, this->transpose()).determinant());
    if(this->_pimpl->column_size == 1)
        return this->at(0, 0);
    if(this->_pimpl->column_size == 2)
        return this->at(0, 0) * this->at(1, 1) - this->at(0, 1) * this->at(1, 0);

    float d = 0;
    for(int i = 0; i < this->_pimpl->column_size; ++i){
        d += this->at(i, 0) * this->minorValue(i, 0) * ::pow(-1, i);
    }
    return d;
}

float Matrix::minorValue(const int &col, const int &row) const
{
    return this->deleteColumnAndRow(col, row).determinant();
}

Matrix Matrix::minorMatrix() const
{
    Matrix m(this->_pimpl->column_size, this->_pimpl->row_size);
    for(int i = 0; i < this->_pimpl->column_size; ++i){
        for(int j = 0; j < this->_pimpl->row_size; ++j){
            m.at(i, j) = this->minorValue(i, j);
        }
    }
    return m;
}

float Matrix::cofactorValue(const int &col, const int &row) const
{
    return this->minorValue(col, row) * ::pow(-1, col + row);
}

Matrix Matrix::cofactorMatrix() const
{
    Matrix m(this->_pimpl->column_size, this->_pimpl->row_size);
    for(int i = 0; i < this->_pimpl->column_size; ++i){
        for(int j = 0; j < this->_pimpl->row_size; ++j){
            m.at(i, j) = this->minorValue(i, j) * ::pow(-1, i + j);
        }
    }
    return m;
}

Matrix Matrix::transpose() const
{
    Matrix m(this->_pimpl->row_size, this->_pimpl->column_size);
    for(int i = 0; i < this->_pimpl->column_size; ++i){
        for(int j = 0; j < this->_pimpl->row_size; ++j){
            m.at(j, i) = this->at(i, j);
        }
    }
    return m;
}

Matrix Matrix::adjoint() const
{
    Matrix m(this->_pimpl->row_size, this->_pimpl->column_size);
    for(int i = 0; i < this->_pimpl->column_size; ++i){
        for(int j = 0; j < this->_pimpl->row_size; ++j){
            m.at(j, i) = this->minorValue(i, j) * ::pow(-1, i + j);
        }
    }
    return m;
}

Matrix Matrix::inverse() const
{
    if(this->_pimpl->column_size != this->_pimpl->row_size)
    {
        Matrix t = this->transpose();
        return Matrix::mul(t, Matrix::mul(*this, t).inverse());
    }
    float d = this->determinant();
    Matrix m(this->_pimpl->row_size, this->_pimpl->column_size);
    for(int i = 0; i < this->_pimpl->column_size; ++i){
        for(int j = 0; j < this->_pimpl->row_size; ++j){
            m.at(j, i) = this->minorValue(i, j) * ::pow(-1, i + j) / d;
        }
    }
    return m;
}

float Matrix::trace() const
{
    assert(this->isSquare());

    float t = 0;
    for(int i = 0; i < this->_pimpl->column_size; ++i)
        t += this->at(i, i);
    return t;
}

Matrix Matrix::deleteColumn(const int &col)
{
    Matrix m(this->_pimpl->column_size - 1, this->_pimpl->row_size);
    for(int j = 0; j < col; ++j)
        for(int i = 0; i < this->_pimpl->row_size; ++i)
            m.at(j, i) = this->at(j, i);

    for(int j = this->_pimpl->column_size - 1; j > col; --j)
        for(int i = 0; i < this->_pimpl->row_size; ++i)
            m.at(j - 1, i) = this->at(j, i);

    return m;
}

Matrix Matrix::deleteRow(const int &row)
{
    Matrix m(this->_pimpl->column_size, this->_pimpl->row_size - 1);

    for(int j = 0; j < this->_pimpl->column_size; ++j)
    {
        for(int i = 0; i < row; ++i)
            m.at(j, i) = this->at(j, i);
        for(int i = this->_pimpl->row_size - 1; i > row; --i)
            m.at(j, i - 1) = this->at(j, i);
    }
    return m;
}

Matrix Matrix::deleteColumnAndRow(const int &col, const int &row) const
{
    Matrix m(this->_pimpl->column_size - 1, this->_pimpl->row_size - 1);

    for(int j = 0; j < col; ++j)
    {
        for(int i = 0; i < row; ++i)
            m.at(j, i) = this->at(j, i);
        for(int i = this->_pimpl->row_size - 1; i > row; --i)
            m.at(j, i - 1) = this->at(j, i);
    }
    for(int j = this->_pimpl->column_size - 1; j > col; --j)
    {
        for(int i = 0; i < row; ++i)
            m.at(j - 1, i) = this->at(j, i);
        for(int i = this->_pimpl->row_size - 1; i > row; --i)
            m.at(j - 1, i - 1)= this->at(j, i);
    }
    return m;
}

Matrix Matrix::add(const Matrix &l, const Matrix &r)
{
    assert(l.columnSize() == r.columnSize() && l.rowSize() == r.rowSize());
    Matrix m(l._pimpl->column_size, l._pimpl->row_size);
    for(int i = 0; i < l._pimpl->column_size; ++i){
        for(int j = 0; j < l._pimpl->row_size; ++j){
            m.at(i, j) = l.at(i, j) + r.at(i, j);
        }
    }
    return m;
}

Matrix Matrix::sub(const Matrix &l, const Matrix &r)
{
    assert(l.columnSize() == r.columnSize() && l.rowSize() == r.rowSize());
    Matrix m(l._pimpl->column_size, l._pimpl->row_size);
    for(int i = 0; i < l._pimpl->column_size; ++i){
        for(int j = 0; j < l._pimpl->row_size; ++j){
            m.at(i, j) = l.at(i, j) - r.at(i, j);
        }
    }
    return m;
}

Vector Matrix::mul(const Matrix &m, const Vector &v)
{
    Vector r(m._pimpl->row_size, true);
    int cols = v.size() < m._pimpl->column_size ? v.size() : m._pimpl->column_size;
    for(int i = 0; i < m._pimpl->row_size; ++i){
        for(int j = 0; j < cols; ++j)
            r[i] += m.at(j, i) * v[j];
    }
    return r;
}

Matrix Matrix::mul(const Matrix &l, const float &r)
{
    Matrix m(l._pimpl->column_size, l._pimpl->row_size, true);
    for(int i = 0; i < l._pimpl->row_size; ++i){
        for(int j = 0; j < l._pimpl->row_size; ++j){
            m.at(i, j) = l.at(i, j) * r;
        }
    }
    return m;
}

Matrix Matrix::mul(const Matrix &l, const Matrix &r)
{
    assert(l.columnSize() == r.rowSize());
    Matrix m(r._pimpl->column_size, l._pimpl->row_size, true);
    for(int i = 0; i < l._pimpl->row_size; ++i){
        for(int j = 0; j < r._pimpl->column_size; ++j){
            for(int k = 0; k < l._pimpl->column_size; ++k){
                m.at(j, i) += l.at(k, i) * r.at(j, k);
            }
        }
    }
    return m;
}

Matrix Matrix::div(const Matrix &l, const float &r)
{
    Matrix m(l._pimpl->column_size, l._pimpl->row_size, true);
    for(int i = 0; i < l._pimpl->row_size; ++i){
        for(int j = 0; j < l._pimpl->row_size; ++j){
            m.at(i, j) = l.at(i, j) / r;
        }
    }
    return m;
}

Matrix Matrix::neg(const Matrix &l)
{
    Matrix m(l._pimpl->column_size, l._pimpl->row_size, true);
    for(int i = 0; i < l._pimpl->row_size; ++i){
        for(int j = 0; j < l._pimpl->row_size; ++j){
            m.at(i, j) = -l.at(i, j);
        }
    }
    return m;
}

bool Matrix::equal(const Matrix &l, const Matrix &r)
{
    static float epsilon = 1e-15;
    if(l._pimpl->row_size != r._pimpl->row_size || l._pimpl->column_size != r._pimpl->column_size)
        return false;

    for(int i = 0; i < l._pimpl->row_size; ++i){
        for(int j = 0; j < l._pimpl->row_size; ++j){
            if(::fabs(l.at(i, j) - r.at(i, j)) >= epsilon)
                return false;
        }
    }
    return true;
}

std::string Matrix::toString() const
{
    std::string str;
    for(int i = 0; i < this->_pimpl->row_size; ++i){
        std::string row_str = "";
        for(int j = 0; j < this->_pimpl->column_size; ++j){
            row_str.append(std::to_string(this->at(j, i))).append("\t");
        }
        str.append(row_str).append("\n");
    }
    return str;
}

Matrix Matrix::copy() const
{
    Matrix m(this->_pimpl->column_size, this->_pimpl->row_size);
    for(int i = 0; i < this->_pimpl->column_size; ++i){
        for(int j = 0; j < this->_pimpl->column_size; ++j){
            m.at(i, j) = this->at(i, j);
        }
    }
    return m;
}

