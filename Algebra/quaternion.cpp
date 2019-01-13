#include "quaternion.h"
#include "math.h"
#include "vector3.h"
#include "matrix3x3.h"
#include "../func.h"

using namespace mp;

Quaternion::Quaternion(const float &w) { this->_d[0] = w; }

Quaternion::Quaternion(const Vector3 &v) { this->_d[1] = v.x(); this->_d[2] = v.y(); this->_d[3] = v.z(); this->_d[0] = 0; }

Quaternion::Quaternion(float w, float x, float y, float z) { this->_d[0] = w; this->_d[1] = x; this->_d[2] = y; this->_d[3] = z; }

Quaternion::Quaternion(float angle, const Vector3 &axis)
{
    Vector3 v = axis.normalize();
    float sa = sin(angle);
    float ca = cos(angle);
    v = Vector3::mul(v, sa);
    float len = sqrt(ca * ca + v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
    this->_d[0] = ca / len;
    this->_d[1] = v.x() / len;
    this->_d[2] = v.y() / len;
    this->_d[3] = v.z() / len;
}

const float &Quaternion::operator [](const int &i) const
{
    return this->_d[i];
}

float &Quaternion::operator [](const int &i)
{
    return this->_d[i];
}

const float &Quaternion::at(const int &i) const
{
    return this->_d[i];
}

float &Quaternion::at(const int &i)
{
    return this->_d[i];
}

float Quaternion::scalar() const { return this->_d[0]; }
Vector3 Quaternion::vector() const { return Vector3(this->_d[1], this->_d[2], this->_d[3]); }
void Quaternion::setScalar(float w) { this->_d[0] = w; }
void Quaternion::setVector(const Vector3 &v) { this->_d[1] = v.x(); this->_d[2] = v.y(); this->_d[3] = v.z(); }

float Quaternion::angle() const
{
    return ::atan(this->_d[0] / this->vector().magnitude());
}

Vector3 Quaternion::axis() const
{
    return this->vector().normalize();
}

Quaternion Quaternion::unit() const
{
    float mag = sqrt(this->_d[0] * this->_d[0] + this->_d[1] * this->_d[1] + this->_d[2] * this->_d[2] + this->_d[3] * this->_d[3]);
    return Quaternion(this->_d[0] / mag, this->_d[1] / mag, this->_d[2] / mag, this->_d[3] / mag);
}

float Quaternion::magnitude() const
{
    return sqrt(this->sqrMagnitude());
}

float Quaternion::sqrMagnitude() const
{
    return this->_d[0] * this->_d[0] + this->_d[1] * this->_d[1] + this->_d[2] * this->_d[2] + this->_d[3] * this->_d[3];
}

Quaternion Quaternion::conjugate() const
{
    return Quaternion(this->_d[0], -this->_d[1], -this->_d[2], - this->_d[3]);
}

Quaternion Quaternion::reciprocal() const
{
    float sLen = this->_d[0] * this->_d[0] + this->_d[1] * this->_d[1] + this->_d[2] * this->_d[2] + this->_d[3] * this->_d[3];
    return Quaternion(this->_d[0] / sLen, -this->_d[1] / sLen, -this->_d[2] / sLen, - this->_d[3] / sLen);
}

Vector3 Quaternion::toEuler() const
{
    float roll, pitch, yaw;
    float sinr = 2.0 * (this->_d[0] * this->_d[1] + this->_d[2] * this->_d[3]);
    float cosr = 1.0 - 2.0 * (this->_d[1] * this->_d[1] + this->_d[2] * this->_d[2]);
    roll = ::atan2(sinr, cosr);

    float sinp = 2.0 * (this->_d[0] * this->_d[2] - this->_d[3] * this->_d[1]);
    if (::fabs(sinp) >= 1)
        pitch = ::copysign(M_PI / 2, sinp); //copysign(x, y) = abs(x) * sign(y)
    else
        pitch = ::asin(sinp);

    float siny = 2.0 * (this->_d[0] * this->_d[3] + this->_d[1] * this->_d[2]);
    float cosy = 1.0 - 2.0 * (this->_d[2] * this->_d[2] + this->_d[3] * this->_d[3]);
    yaw = ::atan2(siny, cosy);
    return Vector3(pitch, roll, yaw);
}

Matrix3x3 Quaternion::toRotationMatrix() const
{
    float s = 1.0 / this->sqrMagnitude();
    Matrix3x3 rm;
    rm.at(0, 0) = 1 - 2 * s * (mp::square(this->y()) + mp::square(this->z()));
    rm.at(0, 1) = 2 * s * (this->x() * this->y() - this->z() * this->w());
    rm.at(0, 2) = 2 * s * (this->x() * this->z() + this->y() * this->w());
    rm.at(1, 0) = 2 * s * (this->x() * this->y() + this->z() * this->w());
    rm.at(1, 1) = 1 - 2 * s * (mp::square(this->x()) + mp::square(this->z()));
    rm.at(1, 2) = 2 * s * (this->y() * this->z() - this->x() * this->w());
    rm.at(2, 0) = 2 * s * (this->x() * this->z() - this->y() * this->w());
    rm.at(2, 1) = 2 * s * (this->y() * this->z() + this->x() * this->w());
    rm.at(2, 2) = 1 - 2 * s * (mp::square(this->x()) + mp::square(this->y()));
    return rm;
}

Quaternion Quaternion::euler(float pitch, float roll, float yaw)
{
    Quaternion q;
    float cy = cos(yaw * 0.5);
    float sy = sin(yaw * 0.5);
    float cr = cos(roll * 0.5);
    float sr = sin(roll * 0.5);
    float cp = cos(pitch * 0.5);
    float sp = sin(pitch * 0.5);

    q._d[0] = cy * cr * cp + sy * sr * sp;
    q._d[1] = cy * sr * cp - sy * cr * sp;
    q._d[2] = cy * cr * sp + sy * sr * cp;
    q._d[3] = sy * cr * cp - cy * sr * sp;
    return q;
}

Quaternion Quaternion::euler(const Vector3 &euler_angles)
{
    return Quaternion::euler(euler_angles.x(), euler_angles.y(), euler_angles.z());
}

Quaternion Quaternion::identity()
{
    Quaternion q;
    q._d[0] = 1;
    return q;
}

Quaternion Quaternion::add(const Quaternion &p, const Quaternion &q)
{
    Quaternion r;
    r._d[0] = p._d[0] + q._d[0];
    r._d[1] = p._d[1] + q._d[1];
    r._d[2] = p._d[2] + q._d[2];
    r._d[3] = p._d[3] + q._d[3];
    return r;
}

Quaternion Quaternion::sub(const Quaternion &p, const Quaternion &q)
{
    Quaternion r;
    r._d[0] = p._d[0] - q._d[0];
    r._d[1] = p._d[1] - q._d[1];
    r._d[2] = p._d[2] - q._d[2];
    r._d[3] = p._d[3] - q._d[3];
    return r;
}

Quaternion Quaternion::mul(const Quaternion &p, const Quaternion &q)
{
    Quaternion r;
    r._d[0] = p._d[0] * q._d[0] - (p.x() * q.x() + p.y() * q.y() + p.z() * q.z());
    Vector3 v1 = Vector3::mul(q.vector(), p._d[0]);
    Vector3 v2 = Vector3::mul(p.vector(), q._d[0]);
    Vector3 pvxqv = Vector3::cross(p.vector(), q.vector());
    Vector3 v = Vector3::add(v1, v2);
    v = Vector3::add(v, pvxqv);
    r._d[1] = v.x();
    r._d[2] = v.y();
    r._d[3] = v.z();
    return r;
}

Quaternion Quaternion::mul(const Vector3 &v, const Quaternion &q)
{
    Quaternion r;
    r._d[0] = -(v.x() * q.x() + v.y() * q.y() + v.z() * q.z());
    Vector3 vector = Vector3::mul(v, q._d[0]);
    Vector3 pvxqv = Vector3::cross(v, q.vector());
    vector = Vector3::add(vector, pvxqv);
    r._d[1] = vector.x();
    r._d[2] = vector.y();
    r._d[3] = vector.z();
    return r;
}

Quaternion Quaternion::mul(const Quaternion &p, const Vector3 &v)
{
    Quaternion r;
    r._d[0] = -(p.x() * v.x() + p.y() * v.y() + p.z() * v.z());
    Vector3 vector = Vector3::mul(v, p._d[0]);
    Vector3 pvxqv = Vector3::cross(p.vector(), v);
    vector = Vector3::add(vector, pvxqv);
    r._d[1] = vector.x();
    r._d[2] = vector.y();
    r._d[3] = vector.z();
    return r;
}

Quaternion Quaternion::mul(const Quaternion &p, float quot)
{
    return Quaternion(p._d[0] * quot, p._d[1] * quot, p._d[2] * quot, p._d[3] * quot);
}

Quaternion Quaternion::div(const Quaternion &p, const Quaternion &q)
{
    return Quaternion::mul(p, q.reciprocal());
}

Quaternion Quaternion::div(const Quaternion &p, float den)
{
    return Quaternion(p._d[0] / den, p._d[1] / den, p._d[2] / den, p._d[3] / den);
}

bool Quaternion::equal(const Quaternion &p, const Quaternion &q)
{
    return ::fabs(p[0] - q[0]) < 1e-15 && ::fabs(p[1] - q[1]) < 1e-15 && ::fabs(p[2] - q[2]) < 1e-15 && ::fabs(p[3] - q[3]) < 1e-15;
}
