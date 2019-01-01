#ifndef FOURIER_H
#define FOURIER_H

#include "func.h"

namespace mp {
class Fourier {

    Range _x, _y;
    Range _a, _b;
    int _it = 0;

public:
    Fourier(const Range &x, const Range &y);

    int iteration() const;

    void iterate(const int &N);
    void iterate();

    const Range &a();
    const Range &b();

    double operator()(const double &x);//f(x) = a0/2 + a(i)cos(ikx) + b(i)sin(ikx) (1 <= i < N)
    Range operator()(const Range &x);

    mp::Range errors(const int &iteration_count);

    Complex transform(const double &k);
    double transformReal(const double &k);
    double transformImaginer(const double &k);
    ComplexRange transform(const Range &k);
    Range transformReal(const Range &k);
    Range transformImaginer(const Range &k);
};
}

#endif // FOURIER_H
