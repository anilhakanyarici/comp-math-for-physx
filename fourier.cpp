#include "fourier.h"

#include <QtMath>

mp::Fourier::Fourier(const mp::Range &x, const mp::Range &y)
{
    this->_x = x;
    this->_y = y;

    double *x_data = this->_x.data();
    double *y_data = this->_y.data();
    int integral_size = this->_x.size();
    double p = x_data[integral_size - 1] - x_data[0];

    double a_integral = 0;
    for(int i = 1; i < integral_size; ++i){
        a_integral += y_data[i] * (x_data[i] - x_data[i - 1]);
    }
    this->_a.append(2 * a_integral / p);
    this->_b.append(0);
    this->_it = 1;
}

int mp::Fourier::iteration() const { return this->_it; }

void mp::Fourier::iterate(const int &N)
{
    mp::Range a(N);
    mp::Range b(N);

    double *x_data = this->_x.data();
    double *y_data = this->_y.data();
    int integral_size = this->_x.size();
    double p = x_data[integral_size - 1] - x_data[0];
    double twopip = 2 * 3.14159265358979323846 / p;

    for(int n = 0; n < N; ++n){
        double n_twopip = n * twopip;
        double a_integral = 0, b_integral = 0, diff = 0, phase = 0, del = 0;
        for(int i = 1; i < integral_size; ++i){
            diff = x_data[i] - x_data[i - 1];
            phase = n_twopip * x_data[i];
            del = y_data[i] * diff;
            a_integral += del * ::cos(phase);
            b_integral += del * ::sin(phase);
        }
        a[n] = (2 * a_integral / p);
        b[n] = (2 * b_integral / p);
    }

    this->_it = N;
    this->_a = a;
    this->_b = b;
}

void mp::Fourier::iterate()
{
    double *x_data = this->_x.data();
    double *y_data = this->_y.data();
    int integral_size = this->_x.size();
    double p = x_data[integral_size - 1] - x_data[0];
    double twopip = 2 * 3.14159265358979323846 / p;

    double it_twopip = this->_it * twopip;
    double a_integral = 0, b_integral = 0, diff = 0, phase = 0, del = 0;
    for(int i = 1; i < integral_size; ++i){
        diff = x_data[i] - x_data[i - 1];
        phase = it_twopip * x_data[i];
        del = y_data[i] * diff;
        a_integral += del * ::cos(phase);
        b_integral += del * ::sin(phase);
    }
    this->_a.append(2 * a_integral / p);
    this->_b.append(2 * b_integral / p);
    ++this->_it;
}

const mp::Range &mp::Fourier::a()
{
    return this->_a;
}

const mp::Range &mp::Fourier::b()
{
    return this->_b;
}

double mp::Fourier::operator()(const double &x)
{
    double s = this->_a[0] / 2;
    double twopi = 2 * 3.14159265358979323846;
    double p = this->_x[this->_x.size() - 1] - this->_x[0];
    for(int i = 1; i < this->_a.size(); ++i){
        s += this->_a[i] * ::cos(twopi * x * i / p) + this->_b[i] * ::sin(twopi * x * i / p);
    }
    return s;
}

mp::Range mp::Fourier::operator()(const mp::Range &x)
{
    mp::Range y(x.size());
    double *y_data = y.data();
    const double *x_data = x.data();
    for(int i = 0; i < x.size(); ++i)
        y_data[i] = (*this)(x_data[i]);
    return y;
}

mp::Range mp::Fourier::errors(const int &iteration_count)
{
    mp::Fourier fourier(this->_x, this->_y);
    mp::Range errs;
    for(int i = 1; i < iteration_count; ++i) {
        mp::Range f = fourier(this->_x);
        errs.append(mp::mse(f, this->_y));
        fourier.iterate();
    }
    return errs;
}

mp::Complex mp::Fourier::transform(const double &k)
{
    double *x_data = this->_x.data();
    double *y_data = this->_y.data();
    int integral_size = this->_x.size();

    double r = 0, i = 0;
    double twopik = mp::twopi * k;
    double phase = twopik * x_data[0];
    double r_prev_integrand =  y_data[0] * ::cos(phase);
    double i_prev_integrand =  y_data[0] * ::sin(phase);
    for(int j = 1; j < integral_size; ++j) {
        phase = twopik * x_data[j];
        double r_curr_integrand = y_data[j] * ::cos(phase);
        double i_curr_integrand = y_data[j] * ::sin(phase);
        double diff = x_data[j] - x_data[j - 1];
        r += (r_curr_integrand + r_prev_integrand) * diff;
        i += (i_curr_integrand + i_prev_integrand) * diff;
        r_prev_integrand = r_curr_integrand;
        i_prev_integrand = i_curr_integrand;
    }
    return Complex(r / 2, i / 2);
}

double mp::Fourier::transformReal(const double &k)
{
    double *x_data = this->_x.data();
    double *y_data = this->_y.data();
    int integral_size = this->_x.size();

    double twopik = mp::twopi * k, integral = 0;
    double prev_integrand =  y_data[0] * ::cos(twopik * x_data[0]);
    for(int j = 1; j < integral_size; ++j) {
        double curr_integrand = y_data[j] * ::cos(twopik * x_data[j]);
        integral += (curr_integrand + prev_integrand) * (x_data[j] - x_data[j - 1]);
        prev_integrand = curr_integrand;
    }
    return integral / 2;
}

double mp::Fourier::transformImaginer(const double &k)
{
    double *x_data = this->_x.data();
    double *y_data = this->_y.data();
    int integral_size = this->_x.size();

    double twopik = mp::twopi * k, integral = 0;
    double prev_integrand =  y_data[0] * ::sin(twopik * x_data[0]);
    for(int j = 1; j < integral_size; ++j) {
        double curr_integrand = y_data[j] * ::sin(twopik * x_data[j]);
        integral += (curr_integrand + prev_integrand) * (x_data[j] - x_data[j - 1]);
        prev_integrand = curr_integrand;
    }
    return integral / 2;
}

mp::ComplexRange mp::Fourier::transform(const mp::Range &k)
{
    mp::ComplexRange f(k.size());
    double *fr_data = f.real().data();
    double *fi_data = f.imaginer().data();
    const double *k_data = k.data();
    for(int n = 0; n < k.size(); ++n){
        Complex t = this->transform(k_data[n]);
        fr_data[n] = t.real();
        fi_data[n] = t.imaginer();
    }
    return f;
}

mp::Range mp::Fourier::transformReal(const mp::Range &k)
{
    mp::Range f(k.size());
    double *f_data = f.data();
    const double *k_data = k.data();
    for(int n = 0; n < k.size(); ++n)
        f_data[n] = this->transformReal(k_data[n]);
    return f;
}

mp::Range mp::Fourier::transformImaginer(const mp::Range &k)
{
    mp::Range f(k.size());
    double *f_data = f.data();
    const double *k_data = k.data();
    for(int n = 0; n < k.size(); ++n)
        f_data[n] = this->transformImaginer(k_data[n]);
    return f;
}
