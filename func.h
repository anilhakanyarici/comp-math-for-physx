#ifndef FUNC_H
#define FUNC_H

#include "range.h"
#include "complex.h"
#include "polynomial.h"

namespace mp {

const double pi         = 3.1415926535897932384;
const double twopi      = 6.2831853071795864769;
const double fourpi     = 12.566370614359172953;
const double sqrtpi     = 3.1415926535897932384;
const double e          = 2.7182818284590452354;
const double ln2        = 0.6931471805599453094;
const double ln3        = 1.0986122886681096913;
const double ln4        = 1.3862943611198906188;
const double ln5        = 1.6094379124341003746;
const double ln6        = 1.7917594692280550008;
const double ln7        = 1.9459101490553133051;
const double ln8        = 2.0794415416798359282;
const double ln9        = 2.1972245773362193827;
const double ln10       = 2.3025850929940456840;
const double invsqrtpi  = 1.1283791670955125739; // [2 / sqrt(pi)]

Range range(const double &a, const double &b, const double &step);
Range rangeCount(const double &a, const double &b, const int &count);
Range linspace(const double &a, const double &b, const int &count);

int sgn(const double &v);
double erf(const double &z);
Range erf(const Range &z);
double abs(const double &v);
Range abs(const Range &x);
double max(const Range &r);
double min(const Range &r);
Range neg(const Range &x);
double factorial(const double &v);
double combination(const int &n, const int &r);

Range linear(const Range &x, double m, double n); //m * x + n
Range poly(const Range &x, const Range &coefs); //poly(x) = sum(coefs[i] * x^i)
Range polyexp(const Range &x, const Range &coefs, double s); //polyexp(x) = poly(x) * exp(-sx)

constexpr double cube(const double &x) { return x * x * x; }
Range cube(const Range &x);
constexpr double square(const double &x) { return x * x; }
Range square(const Range &x);
double sqrt(const double &x);
Range sqrt(const Range &x);
double exp(const double &x);
Range exp(const Range &x);
double exp2(const double &x);
Range exp2(const Range &x);
double exp10(const double &x);
Range exp10(const Range &x);
double pow(const double &x, const double &f);
Range pow(const Range &x, double f);

double sin(const double &x);
Range sin(const Range &x);
double cos(const double &x);
Range cos(const Range &x);
double tan(const double &x);
Range tan(const Range &x);
double cot(const double &x);
Range cot(const Range &x);
double sec(const double &x);
Range sec(const Range &x);
double csc(const double &x);
Range csc(const Range &x);

double sinh(const double &x);
Range sinh(const Range &x);
double cosh(const double &x);
Range cosh(const Range &x);
double tanh(const double &x);
Range tanh(const Range &x);
double coth(const double &x);
Range coth(const Range &x);

double lb(const double &x);
Range lb(const Range &x);
double ln(const double &x);
Range ln(const Range &x);
double ld(const double &x);
Range ld(const Range &x);
double log(const double &base, const double &x);
Range log(double base, const Range &x);

double gamma(const double &z);
double beta(const double &m, const double &n);
double zeta(const double &x);
const Range gamma(const Range &x);
const Range zeta(const Range &x);

double sum(const Range &x);
double mean(const Range &x);
double dev(const Range &x); //deviation
double var(const Range &x); //variance
double mse(const Range &x2, const Range &x1); //mean-square error
double cov(const Range &x, const Range &y); //covariance
double cor(const Range &x, const Range &y); //correlation

Range hamming(const int &N, const double &a = 0.53836);
Range hamming(const Range &r, const double &a = 0.53836);
Range smooth(const mp::Range &v, const int &radius, const int &times = 1);
Range normalize(const mp::Range &v);

void linreg(const Range &x, const Range &y, double &m, double &n); //linear-regression parameters
Range linreg(const Range &x, const Range &y, double *m = nullptr, double *n = nullptr); //linear-regression parameters with linear function.

double integral(const Range &x, const Range &y);
QMap<double, double> frequency(const Range &values, const int &bins = 10);
mp::Range random(int count = 1, double min = 0, double max = 1);

QVector<Range> loadFromCSVFile(const QString &filename, QStringList *column_names = nullptr);
}

#endif // FUNC_H
