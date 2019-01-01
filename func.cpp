#include "func.h"

#include <QFile>
#include <QtMath>
#include <QTextStream>

mp::Range mp::range(const double &a, const double &b, const double &step)
{
    double min = ::fmin(a, b);
    double max = ::fmax(a, b);
    int count = (int)::qRound((max - min) / step);
    mp::Range y(count);
    double *y_data = y.data();
    for(int i = 0; i < count; ++i)
        y_data[i] = min + i * step;
    return y;
}

mp::Range mp::rangeCount(const double &a, const double &b, const int &count)
{
    double min = ::fmin(a, b);
    double max = ::fmax(a, b);
    double step = (max - min) / count;
    mp::Range y(count);
    double *y_data = y.data();
    for(int i = 0; i < count; ++i)
        y_data[i] = min + i * step;
    return y;
}

mp::Range mp::linspace(const double &a, const double &b, const int &count)
{
    double min = ::fmin(a, b);
    double max = ::fmax(a, b);
    double step = (max - min) / (count - 1);
    mp::Range y(count);
    double *y_data = y.data();
    for(int i = 0; i < count; ++i)
        y_data[i] = min + i * step;
    return y;
}

int mp::sgn(const double &v)
{
    return (0 < v) - (v < 0);
}

double mp::erf(const double &z)
{
    int sign = mp::sgn(z);
    double to = mp::abs(z);
    double dt = to / 1000;
    double t = 0;
    double integral = 0;
    while(t < to) {
        double g1 = ::exp(-(t * t));
        t += dt;
        double g2 = ::exp(-(t * t));
        double g = (g1 + g2) / 2;
        integral += g * dt;
    }
    return sign * integral * mp::invsqrtpi;
}

mp::Range mp::erf(const mp::Range &z)
{
    mp::Range y(z.size());
    double *y_data = y.data();
    const double *z_data = z.data();
    for(int i = 0; i < z.size(); ++i)
        y_data[i] = mp::erf(z_data[i]);
    return y;
}

double mp::abs(const double &v)
{
    return v < 0 ? -v : v;
}

mp::Range mp::abs(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::fabs(x_data[i]);
    return y;
}

double mp::max(const mp::Range &r)
{
    const double *r_data = r.data();
    double max = r_data[0];
    for(int i = 1; i < r.size(); ++i)
        if(max < r_data[i])
            max = r_data[i];
    return max;
}

double mp::min(const mp::Range &r)
{
    const double *r_data = r.data();
    double min = r_data[0];
    for(int i = 1; i < r.size(); ++i)
        if(min > r_data[i])
            min = r_data[i];
    return min;
}

mp::Range mp::neg(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = -x_data[i];
    return y;
}

mp::Range mp::linear(const mp::Range &x, double m, double n)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = x_data[i] * m + n;
    return y;
}

mp::Range mp::poly(const mp::Range &x, const mp::Range &coefs)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i) {
        double y_term = 0;
        for(int j = 0; j < coefs.size(); ++j) {
            y_term += ::pow(x_data[i], j) * coefs[j];
        }
        y_data[i] = y_term;
    }
    return y;
}

mp::Range mp::polyexp(const mp::Range &x, const mp::Range &coefs, double s)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i) {
        double y_term = 0;
        for(int j = 0; j < coefs.size(); ++j) {
            y_term += ::pow(x_data[i], j) * coefs[j];
        }
        y_data[i] = y_term * ::exp(-s * x_data[i]);
    }
    return y;
}

mp::Range mp::cube(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = x_data[i] * x_data[i] * x_data[i];
    return y;
}

mp::Range mp::square(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = x_data[i] * x_data[i];
    return y;
}

double mp::sqrt(const double &x) { return ::sqrt(x); }

mp::Range mp::sqrt(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::sqrt(x_data[i]);
    return y;
}

double mp::exp(const double &x) { return ::exp(x); }

mp::Range mp::exp(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::exp(x_data[i]);
    return y;
}

double mp::exp2(const double &x) { return ::exp2(x); }

mp::Range mp::exp2(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::exp2(x_data[i]);
    return y;
}

double mp::exp10(const double &x) { return exp10(x); }

mp::Range mp::exp10(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::exp10(x_data[i]);
    return y;
}

double mp::pow(const double &x, const double &f) { return ::pow(x, f); }

mp::Range mp::pow(const mp::Range &x, double f)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::pow(x_data[i], f);
    return y;
}

double mp::sin(const double &x) { return ::sin(x); }

mp::Range mp::sin(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::qSin(x_data[i]);
    return y;
}

double mp::cos(const double &x) { return ::cos(x); }

mp::Range mp::cos(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::qCos(x_data[i]);
    return y;
}

double mp::tan(const double &x) { return ::tan(x); }

mp::Range mp::tan(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::qTan(x_data[i]);
    return y;
}

double mp::cot(const double &x) { return 1.0 / ::tan(x); }

mp::Range mp::cot(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = 1.0 / ::qTan(x_data[i]);
    return y;
}

double mp::sec(const double &x) { return 1.0 / ::cos(x); }

mp::Range mp::sec(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = 1.0 / ::cos(x_data[i]);
    return y;
}

double mp::csc(const double &x) { return 1.0 / ::sin(x); }

mp::Range mp::csc(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = 1.0 / ::sin(x_data[i]);
    return y;
}

double mp::sinh(const double &x) { return ::sinh(x); }

mp::Range mp::sinh(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::sinh(x_data[i]);
    return y;
}

double mp::cosh(const double &x) { return ::cosh(x); }

mp::Range mp::cosh(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::cosh(x_data[i]);
    return y;
}

double mp::tanh(const double &x) { return ::tanh(x); }

mp::Range mp::tanh(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::tanh(x_data[i]);
    return y;
}

double mp::coth(const double &x) { return 1.0 / ::tanh(x); }

mp::Range mp::coth(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = 1.0 / ::tanh(x_data[i]);
    return y;
}

double mp::lb(const double &x) { return ::qLn(x) / mp::ln2; }

mp::Range mp::lb(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::qLn(x_data[i]) / mp::ln2;
    return y;
}

double mp::ln(const double &x) { return ::qLn(x); }

mp::Range mp::ln(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::qLn(x_data[i]);
    return y;
}

double mp::ld(const double &x) { return ::qLn(x) / mp::ln10; }

mp::Range mp::ld(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::qLn(x_data[i]) / mp::ln10;
    return y;
}

double mp::log(const double &base, const double &x) { return ::qLn(x) / ::qLn(base); }

mp::Range mp::log(double base, const mp::Range &x)
{
    double log_base = ::qLn(base);
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = ::qLn(x_data[i]) / log_base;
    return y;
}

double mp::gamma(const double &z)
{
    double dt = z / 100;
    double to = z * 10;
    double t = dt;
    double integral = dt;
    while(t < to) {
        double g1 = ::exp(-t) * ::pow(t, z - 1);
        t += dt;
        double g2 = ::exp(-t) * ::pow(t, z - 1);
        double g = (g1 + g2);
        integral += g * dt;
    }
    return integral / 2;
}

double mp::beta(const double &m, const double &n)
{
    double du = 1 / 1000;
    double u = 0;
    double integral = 0;
    while(u < 1) {
        double g1 = ::pow(u, m - 1) * ::pow(1 - u, n - 1);
        u += du;
        double g2 = ::pow(u, m - 1) * ::pow(1 - u, n - 1);
        double g = (g1 + g2) / 2;
        integral += g * du;
    }
    return integral;
}

double mp::zeta(const double &x)
{
    if(x <= 1)
        return INFINITY;
    double sum = 0;
    for(int i = 1; i <= 100; ++i) {
        sum += ::pow(i, -x);
    }
    return sum;
}

const mp::Range mp::gamma(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = mp::gamma(x_data[i]);
    return y;
}

const mp::Range mp::zeta(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = mp::zeta(x_data[i]);
    return y;
}

double mp::sum(const mp::Range &x)
{
    double sum = 0;
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        sum += x_data[i];
    return sum;
}

double mp::mean(const mp::Range &x)
{
    double sum = 0;
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        sum += x_data[i];
    return sum / x.length();
}

double mp::dev(const mp::Range &x)
{
    return mp::sqrt(mp::var(x));
}

double mp::var(const mp::Range &x)
{
    double mean = mp::mean(x);
    double var = 0;
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i) {
        double dif = x_data[i] - mean;
        var += dif * dif;
    }
    return var / x.length();
}

double mp::mse(const mp::Range &x2, const mp::Range &x1)
{
    double mse = 0;
    const double *x1_data = x1.data();
    const double *x2_data = x2.data();
    for(int i = 0; i < x1.size(); ++i) {
        double dif = x1_data[i] - x2_data[i];
        mse += dif * dif;
    }
    return mse / x1.length();
}

double mp::cov(const mp::Range &x, const mp::Range &y)
{
    double x_mean = mp::mean(x);
    double y_mean = mp::mean(y);

    double cov = 0;
    const double *x_data = x.data();
    const double *y_data = y.data();
    for(int i = 0; i < x.size(); ++i) {
        double dif_x = x_data[i] - x_mean;
        double dif_y = y_data[i] - y_mean;
        cov += dif_x * dif_y;
    }
    return cov / x.length();
}

mp::Range mp::hamming(const int &N, const double &a)
{
    assert(a <= 1);
    double a1 = 1 - a;
    mp::Range h(N);
    double *h_data = h.data();
    double k = 2 * mp::pi / (N - 1);
    for(int i = 0; i < N; ++i) {
        h_data[i] = a - a1 * mp::cos(k * i);
    }
    return h;
}

mp::Range mp::hamming(const Range &r, const double &a)
{
    assert(a <= 1);
    int N = r.size();
    double a1 = 1 - a;
    mp::Range h(N);
    double *h_data = h.data();
    const double *r_data = r.data();
    double k = 2 * mp::pi / (N - 1);
    for(int i = 0; i < N; ++i) {
        h_data[i] = a - a1 * mp::cos(k * r_data[i]);
    }
    return h;
}

mp::Range mp::smooth(const mp::Range &v, const int &point, const int &times)
{
    mp::Range smoother = mp::hamming(2 * point + 1);
    smoother = smoother / mp::sum(smoother);

    mp::Range s(v.size());
    mp::Range c = v.copy();
    for(int i = 0; i < point; ++i){
        s[i] = v[i];
        s[v.size() - 1 - i] = v[v.size() - 1 - i];
    }
    int middles = s.size() - point;
    for(int k = 0; k < times; ++k){
        double *s_data = s.data();
        double *c_data = c.data();
        for(int i = point; i < middles; ++i) {
            s_data[i] = 0;
            for(int j = -point; j <= point; ++j) {
                s_data[i] += c_data[i + j] * smoother[j + point];
            }
        }
        c = s;
    }
    return c;
}

void mp::linreg(const mp::Range &x, const mp::Range &y, double &m, double &n)
{
    m = mp::cov(x, y) / mp::var(x);
    n = mp::mean(y) - m * mp::mean(x);
}

mp::Range mp::linreg(const mp::Range &x, const mp::Range &y, double *m, double *n)
{
    double m0, n0;
    if(m == nullptr) m = &m0;
    if(n == nullptr) n = &n0;
    *m = mp::cov(x, y) / mp::var(x);
    *n = mp::mean(y) - *m * mp::mean(x);
    return mp::linear(x, *m, *n);
}

double mp::integral(const mp::Range &x, const mp::Range &y)
{
    int i_last = x.size() - 1;
    const double *x_data = x.data();
    const double *y_data = y.data();
    double dx, integral = 0;
    for(int i = 0; i < i_last; ++i) {
        dx = x_data[i + 1] - x_data[i];
        integral += dx * (y_data[i + 1] + y_data[i]);
    }
    return integral / 2;
}

QMap<double, double> mp::frequency(const Range &values, const int &bins)
{
    assert(bins > 0);
    QMap<double, double> freq;
    double min = mp::min(values);
    double max = mp::max(values);
    double diff = ::ceil(max - min) / bins;

    for(int i = 0; i < values.size(); ++i){
        double b = (int)(values[i] / diff) * diff;
        ++freq[b];
    }
    return freq;
}

mp::Range mp::random(int count, double min, double max)
{
    assert(min < max);
    assert(count > 0);
    mp::Range rnds(count);
    double *rnds_data = rnds.data();
    for(int i = 0; i < count; ++i){
        double r = (double)qrand() / RAND_MAX;
        rnds_data[i] = min + r * (max - min);
    }
    return rnds;
}

QVector<mp::Range> mp::loadFromCSVFile(const QString &filename, QStringList *column_names)
{
    QVector<mp::Range> matrix;
    QFile file(filename);
    if(!file.open(QFile::ReadOnly))
        return matrix;
    if(column_names != nullptr)
        column_names->clear();

    QString line = QString(file.readLine()).trimmed();
    while(line.isEmpty())
        line = QString(file.readLine()).trimmed();

    bool is_first_line = true;
    while(!line.isEmpty()){
        int i = 0;
        int column = 0;
        while(i < line.size()){
            QString value;
            for(; i < line.size(); ++i) {
                QChar c = line[i];
                if(c != QChar('\r') && c != QChar(' ') && c != QChar('\t')) {
                    break;
                }
            }
            for(; i < line.size(); ++i){
                QChar c = line[i];
                if(c == QChar('\r') || c == QChar(' ') || c == QChar('\t'))  {
                    break;
                }
                value.append(c);
            }
            while(column >= matrix.size())
                matrix.append(mp::Range());
            bool ok = false;
            float f_value = value.toDouble(&ok);
            if(column_names != nullptr && is_first_line) {
                if(ok)
                    column_names->append("");
                else {
                    column_names->append(value);
                    continue;
                }
            }
            matrix[column].append(f_value);
            assert(ok || is_first_line);
            ++column;
        }
        is_first_line = false;
        line = QString(file.readLine()).trimmed();
    }
    return matrix;
}
